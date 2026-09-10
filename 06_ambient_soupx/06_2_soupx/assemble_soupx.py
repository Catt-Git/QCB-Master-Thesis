#!/usr/bin/env python3
"""06_2 step 3: the 148 per-channel corrections put back into one object.

Reads `shiao.h5ad`, subtracts what `run_soupx.R` removed in each channel, and writes the
result as a full-size object with the same cells, the same genes and the same `.obs` as
phase 01 produced - the only difference being that the counts have had their ambient
component taken out, and that every cell now carries how much of it that was.

  $DATA_DIR/06_amb/shiao_soupx_all_cells.h5ad     619,693 cells, nothing dropped yet

**Nothing is dropped here.** Which cells are too far gone to keep is a threshold read off a
distribution, so it is taken in `06_3_ambient_qc/ambient_qc.ipynb` and applied by 06_4 -
the same split as 05_1, where `run_infercnv.R` produces the residuals and
`call_malignant.ipynb` decides what they mean.

**What `.X` is in the output, and why it is not what it was.** In `shiao.h5ad`, `.X` is the
scran log-normalised matrix and `obs['size_factors']` are the factors that produced it. Both
describe counts that no longer exist: correcting the counts invalidates the size factors that
were estimated from them, and re-running scran on 619,693 cells is the 470 GB job of
`01_pre_processing/submit_preprocessing_all.slurm`. Rather than carry a stale normalisation
that would look valid, this object holds the **corrected raw counts in both `.X` and
`layers['counts']`**, `size_factors` is dropped, and `uns['ambient_soupx']` records it. That
costs nothing downstream: `05_2/subset_and_qc.ipynb` copies `layers['counts']` into `.X` and
drops the size factors as its first act anyway, and it re-runs scran on the subset it builds.

**What is kept although it was computed before the correction.** The phase-01 PCA, UMAP,
neighbour graph, leiden sweep, cell-cycle scores and scrublet columns all stay, unchanged and
under their original names. They are landmarks - the same status `cell_type_01_4` has in
phase 05 - and the notebook uses the UMAP to show *where* the soup sits. The two things that
would be read as current and are not, `size_factors` and `var['highly_variable']`, are
removed instead. The QC columns are the exception: `total_counts`, `n_genes_by_counts`,
`pct_counts_mt` and the rest are RECOMPUTED on the corrected counts, because 05_2 filters on
them, and the phase-01 values are kept beside them with a `_precorrection` suffix so the
notebook can plot one against the other.

New `.obs` columns, all prefixed `soupx_`:
    soupx_rho            the channel's contamination fraction (constant within a sample)
    soupx_method         autoEstCont / fixed / fallback_prior, from run_soupx.R
    soupx_umis_before    UMIs in this cell before the correction
    soupx_umis_after     UMIs after
    soupx_frac_removed   1 - after/before. THE per-cell quantity 06_3 thresholds on:
                         rho is a property of the channel, this is a property of the cell,
                         and a cell that is mostly soup loses most of its counts here.

Input : $DATA_DIR/shiao.h5ad
        $DATA_DIR/06_amb/adjusted/<sample>/{removed.mtx.gz,rho.csv,cells.csv}
Output: $DATA_DIR/06_amb/shiao_soupx_all_cells.h5ad
        $DATA_DIR/06_amb/soupx_rho.csv      one row per channel
        $DATA_DIR/06_amb/soupx_cells.csv    one row per cell (what the notebook reads)

Local usage (benchmark-py-r), normally through soupx_all.sh:
    export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
    python3 assemble_soupx.py
    python3 assemble_soupx.py --allow-missing   # assemble what is there, report the rest
"""

from __future__ import annotations

import argparse
import gzip
import os
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io as sio
import scipy.sparse as sp

SAMPLE_KEY = "sample"
# Recomputed on the corrected counts; the phase-01 value is kept as <name>_precorrection.
QC_OBS = ["n_genes_by_counts", "total_counts", "log1p_n_genes_by_counts",
          "log1p_total_counts", "total_counts_mt", "log1p_total_counts_mt",
          "pct_counts_mt", "total_counts_ribo", "log1p_total_counts_ribo",
          "pct_counts_ribo"]
# Dropped: both describe the uncorrected counts and would be read as current.
DROP_OBS = ["size_factors"]
DROP_VAR = ["highly_variable"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--allow-missing", action="store_true",
                   help="assemble even if some channels have no removals (they are carried "
                        "through UNCORRECTED and flagged soupx_method='not_run')")
    p.add_argument("--out-name", default="shiao_soupx_all_cells.h5ad",
                   help="output file name under 06_amb/ [%(default)s]")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    data_dir = Path(os.environ["DATA_DIR"]).expanduser().resolve()
    in_path = data_dir / "shiao.h5ad"
    amb_dir = data_dir / "06_amb"
    adj_dir = amb_dir / "adjusted"
    out_path = amb_dir / args.out_name

    assert in_path.exists(), f"missing {in_path} (produced by phase 01)"
    assert adj_dir.is_dir(), f"missing {adj_dir}; run the soupx stage first"

    print(f"DATA_DIR : {data_dir}")
    print(f"Input    : {in_path}")
    print(f"Adjusted : {adj_dir}")
    print(f"Output   : {out_path}")
    print()

    print("Loading data...", flush=True)
    adata = sc.read_h5ad(in_path)
    assert "counts" in adata.layers, "expected raw counts in layers['counts']"
    counts = adata.layers["counts"]
    counts = counts.tocsr() if not sp.isspmatrix_csr(counts) else counts
    print(f"{adata.n_obs:,} cells x {adata.n_vars:,} genes, "
          f"{counts.nnz:,} non-zero counts", flush=True)

    var_names = pd.Index(adata.var_names)
    sample_all = adata.obs[SAMPLE_KEY].astype(str)
    samples = list(adata.obs[SAMPLE_KEY].cat.categories)

    have = [s for s in samples if (adj_dir / s / "removed.mtx.gz").exists()]
    missing = [s for s in samples if s not in set(have)]
    if missing and not args.allow_missing:
        print(f"{len(missing)} channel(s) have no removals in {adj_dir}:", file=sys.stderr)
        print("  " + ", ".join(missing[:20]) + (" ..." if len(missing) > 20 else ""),
              file=sys.stderr)
        print("Run the soupx stage for them, or pass --allow-missing to carry them through "
              "uncorrected.", file=sys.stderr)
        return 1
    print(f"channels: {len(have)} corrected, {len(missing)} without removals\n")

    # -----------------------------------------------------------------------------------
    # Apply the deltas. One channel at a time, in place on the counts matrix's own row
    # blocks: the cells of a sample are contiguous in this object, so each block is a
    # slice assignment rather than a fancy-index scatter, and the 916M-entry matrix is
    # never duplicated.
    # -----------------------------------------------------------------------------------
    blocks = []
    order = []
    for i, sample in enumerate(samples, 1):
        pos = np.flatnonzero((sample_all == sample).to_numpy())
        if len(pos) == 0:
            continue
        block = counts[pos, :]
        if sample in set(have):
            with gzip.open(adj_dir / sample / "removed.mtx.gz", "rb") as fh:
                removed = sio.mmread(fh)
            # run_soupx.R writes genes x cells, as SoupX and 10x do; the object is the
            # transpose of that.
            removed = sp.csr_matrix(removed).T.tocsr()
            assert removed.shape == block.shape, (
                f"{sample}: removals are {removed.shape}, the object's block is "
                f"{block.shape}; the input was prepared from a different object"
            )
            block = block - removed
            block.eliminate_zeros()
            neg = int((block.data < 0).sum()) if block.nnz else 0
            assert neg == 0, f"{sample}: {neg:,} negative counts after subtraction"
            n_removed = int(removed.sum())
        else:
            n_removed = 0
        blocks.append(sp.csr_matrix(block, dtype=np.float32))
        order.append(pos)
        if i % 20 == 0 or i == len(samples):
            print(f"  applied {i}/{len(samples)} channels", flush=True)
        del block

    order = np.concatenate(order)
    corrected = sp.vstack(blocks, format="csr")
    del blocks, counts
    if not np.all(np.diff(order) > 0):
        # The blocks were built in category order; if that is not the object's cell order,
        # put them back. On this object it is (each sample is one contiguous run), so this
        # branch normally does not fire - it is here so the script is not silently wrong on
        # an object where it does not hold.
        corrected = corrected[np.argsort(order), :]
    corrected.eliminate_zeros()

    removed_total = float(adata.layers["counts"].sum() - corrected.sum())
    before_total = float(adata.layers["counts"].sum())
    print(f"\nremoved {removed_total:,.0f} of {before_total:,.0f} UMIs "
          f"({100 * removed_total / before_total:.2f}%), "
          f"{corrected.nnz:,} non-zero left ({100 * corrected.nnz / adata.layers['counts'].nnz:.1f}% "
          "of the entries)\n", flush=True)

    # -----------------------------------------------------------------------------------
    # The per-cell and per-channel tables written by the R stage, merged and attached.
    # -----------------------------------------------------------------------------------
    rho = pd.concat([pd.read_csv(adj_dir / s / "rho.csv") for s in have], ignore_index=True)
    rho.to_csv(amb_dir / "soupx_rho.csv", index=False)
    print(f"rho: median {rho['rho'].median():.3f}, range "
          f"{rho['rho'].min():.3f}-{rho['rho'].max():.3f}")
    print(rho["method"].value_counts().to_string())
    print()

    cells = pd.concat([pd.read_csv(adj_dir / s / "cells.csv").assign(sample=s) for s in have],
                      ignore_index=True).set_index("cell")
    cells = cells.reindex(adata.obs_names)
    # Channels without removals: unchanged counts, and said so rather than left as NaN.
    untouched = cells["umis_before"].isna().to_numpy()
    if untouched.any():
        raw_totals = np.asarray(adata.layers["counts"].sum(axis=1)).ravel()
        cells.loc[untouched, "umis_before"] = raw_totals[untouched]
        cells.loc[untouched, "umis_after"] = raw_totals[untouched]

    rho_by_sample = rho.set_index("sample")
    adata.obs["soupx_rho"] = sample_all.map(rho_by_sample["rho"]).astype(float).to_numpy()
    adata.obs["soupx_method"] = pd.Categorical(
        sample_all.map(rho_by_sample["method"]).fillna("not_run").to_numpy())
    adata.obs["soupx_umis_before"] = cells["umis_before"].to_numpy()
    adata.obs["soupx_umis_after"] = cells["umis_after"].to_numpy()
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = 1.0 - (adata.obs["soupx_umis_after"].to_numpy()
                      / adata.obs["soupx_umis_before"].to_numpy())
    # A cell with zero counts before the correction cannot have lost a fraction of them.
    adata.obs["soupx_frac_removed"] = np.where(
        adata.obs["soupx_umis_before"].to_numpy() > 0, frac, 0.0)

    cells_out = adata.obs[["sample", "cohort", "cell_type", "soupx_rho", "soupx_method",
                           "soupx_umis_before", "soupx_umis_after", "soupx_frac_removed"]]
    cells_out.to_csv(amb_dir / "soupx_cells.csv")
    print(f"per-cell table -> {amb_dir / 'soupx_cells.csv'}")
    print(adata.obs["soupx_frac_removed"].describe().to_string())
    print()

    # -----------------------------------------------------------------------------------
    # The object. See the docstring for why .X becomes the corrected counts and why the
    # phase-01 QC columns are recomputed rather than kept.
    # -----------------------------------------------------------------------------------
    for col in QC_OBS:
        if col in adata.obs:
            adata.obs[f"{col}_precorrection"] = adata.obs[col].to_numpy()
    for col in DROP_OBS:
        if col in adata.obs:
            adata.obs.drop(columns=[col], inplace=True)
            print(f"obs dropped: {col} (estimated from the uncorrected counts)")
    for col in DROP_VAR:
        if col in adata.var:
            adata.var.drop(columns=[col], inplace=True)
            print(f"var dropped: {col} (selected on the uncorrected counts)")

    adata.X = corrected
    adata.layers["counts"] = corrected
    del corrected

    print("\nRecomputing the QC metrics on the corrected counts...", flush=True)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], percent_top=None,
                               log1p=True, inplace=True)

    adata.uns["ambient_soupx"] = {
        "tool": "SoupX",
        "soupx_version": str(rho["soupx_version"].iloc[0]) if len(rho) else "",
        "n_channels_corrected": len(have),
        "n_channels_not_run": len(missing),
        "umis_before": before_total,
        "umis_removed": removed_total,
        # Spelled out in the object itself, so anyone who opens this file without the README
        # still learns that .X is not what it is in shiao.h5ad.
        "X_is": "corrected raw counts (identical to layers['counts'])",
        "size_factors": "dropped: the phase-01 scran factors describe the uncorrected counts; "
                        "05_2 re-runs scran on the subset it builds",
        "qc_metrics": "recomputed on the corrected counts; the phase-01 values are kept "
                      "under <name>_precorrection",
        "stale_but_kept": "obsm['X_pca'], obsm['X_umap'], the neighbour graph, the leiden "
                          "sweep, the cell-cycle scores and the scrublet columns were all "
                          "computed before the correction and are kept as landmarks",
    }

    x_values = adata.X.data
    assert np.allclose(x_values, np.round(x_values)), \
        ".X must hold integer counts; run_soupx.R was run without roundToInt"
    assert x_values.min() >= 0, "negative counts in the corrected matrix"

    print(f"\nWriting {out_path} ...", flush=True)
    adata.write_h5ad(out_path, compression="gzip")
    print(f"done: {out_path} ({out_path.stat().st_size / 1024**3:.2f} GB)")
    print("\nnext: 06_3_ambient_qc/ambient_qc.ipynb")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
