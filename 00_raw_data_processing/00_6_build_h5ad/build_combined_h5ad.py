#!/usr/bin/env python3

# build_combined_h5ad.py
# Load filtered_feature_bc_matrix.h5 for each sample in unique_samples.txt, attach the per-sample metadata (batch, cohort, treatment, fraction, dataset_origin, response) from sample_metadata_final.tsv, and combine everything into a single AnnData.
#
# Object contract (for Methods): 
# - input filtered_feature_bc_matrix.h5 (Cell Ranger 4.0.0 cell-calling already applied), so counts differ from the paper BY DESIGN
# - var_names = gene symbols made unique with var_names_make_unique() (CellTypist Cells_Adult_Breast.pkl works on symbols)
# - Ensembl IDs kept in .var['gene_ids'] via merge='same'
# - raw counts stored in .layers['counts'] before any downstream normalisation
# - barcodes made unique per sample as f"{sample}_{bc}", so concat uses index_unique=None.
#
# Usage:
#     python3 build_combined_h5ad.py \
#         --cellranger-dir cellranger_out \
#         --samples-file unique_samples.txt \
#         --metadata-tsv sample_metadata_final.tsv \
#         --out all_samples_combined.h5ad

import argparse
import scanpy as sc
import anndata as ad
import pandas as pd
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cellranger-dir", default="cellranger_out")
    ap.add_argument("--samples-file", default="unique_samples.txt")
    ap.add_argument("--metadata-tsv", default="sample_metadata_final.tsv")
    ap.add_argument("--out", default="all_samples_combined.h5ad")
    args = ap.parse_args()
    # Parse input/output paths

    print(f"anndata {ad.__version__} | scanpy {sc.__version__}")
    # Print library versions to record them for Methods (concat/merge behaviour is version-dependent)

    cellranger_dir = Path(args.cellranger_dir)
    meta_df = pd.read_csv(args.metadata_tsv, sep="\t").set_index("sample_id")
    # Load the metadata table, indexed by sample_id for per-sample lookup

    with open(args.samples_file) as f:
        samples = [line.strip() for line in f if line.strip()]
    print(f"Samples to process: {len(samples)}")
    # Read the list of samples to combine

    adatas = []
    missing, unmatched_meta = [], []
    # Use a list: passing a dict to ad.concat makes anndata treat the
    # keys as labels, and the exact behaviour is version-dependent; a list is
    # explicit and version-robust. Trackers record skipped / unmatched samples

    for sample in samples:
        h5_path = cellranger_dir / sample / "outs" / "filtered_feature_bc_matrix.h5"
        if not h5_path.exists():
            print(f"[SKIP] Missing output: {sample}")
            missing.append(sample)
            continue
        # Locate each sample's filtered matrix; skip and record samples without output.

        adata = sc.read_10x_h5(h5_path)
        adata.var_names_make_unique()
        # Read the 10x matrix and make gene symbols unique (kept as var_names for CellTypist).

        adata.layers["counts"] = adata.X.copy()
        # Preserve the raw integer counts in a dedicated layer BEFORE any downstream
        # normalisation overwrites .X. Required by scVI/scANVI/DRVI (need raw counts)
        # and by several scib metrics; recovering them later would mean re-reading all .h5.

        adata.obs_names = [f"{sample}_{bc}" for bc in adata.obs_names]
        # Prefix barcodes with the sample name so they are unique across samples before the merge.

        if sample in meta_df.index:
            row = meta_df.loc[sample]
            adata.obs["sample"] = sample
            adata.obs["batch"] = row["batch"]
            adata.obs["cohort"] = row["cohort"]
            adata.obs["treatment"] = row["treatment"]
            adata.obs["fraction"] = row["fraction"]
            adata.obs["dataset_origin"] = row["dataset_origin"]
            adata.obs["response"] = row["response"]
        else:
            print(f"[WARN] No metadata found for: {sample}")
            unmatched_meta.append(sample)
            adata.obs["sample"] = sample
        # Inject the per-sample metadata into .obs; for samples missing from the table,
        # still set .obs['sample'] so the column exists after concat.

        adatas.append(adata)
        print(f"[OK] {sample}: {adata.n_obs} cells, {adata.n_vars} genes")
        # Keep per-sample logging for traceability of the run.

    print("\n--- Summary ---")
    print(f"Samples loaded successfully: {len(adatas)}")
    if missing:
        print(f"Missing samples (h5 not found): {missing}")
    if unmatched_meta:
        print(f"Samples without clinical metadata: {unmatched_meta}")
    # Report load outcome and any problems before merging.

    print("\nMerging into a single AnnData...")
    combined = ad.concat(adatas, join="outer", index_unique=None, merge="same")
    # Concatenate all samples: outer join on genes; index_unique=None because barcodes
    # are already unique (prefixed above); merge='same' propagates .var columns shared
    # across samples (notably .var['gene_ids'], the Ensembl IDs) into the combined object.

    combined.write_h5ad(args.out, compression="gzip")
    print(f"\nFinal total: {combined.n_obs} cells, {combined.n_vars} genes")
    print(f"Saved to: {args.out}")
    # Write the combined object and report its final dimensions.


if __name__ == "__main__":
    main()