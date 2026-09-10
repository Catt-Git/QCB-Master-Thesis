#!/usr/bin/env python3
"""06_2 step 1: one SoupX input directory per sample, out of shiao.h5ad.

SoupX works per **channel**, i.e. per 10x run, because the soup is a property of the
emulsion: the ambient RNA in `P01_A_P` came from the cells lysed in that tube and nowhere
else. `shiao.h5ad` holds 148 of those channels in one object, so this script takes it apart
along `obs['sample']` and hands each piece to R with everything the correction needs:

  input/genes.tsv                 the object's 30,869 gene symbols, ONE copy, shared
  input/<sample>/counts.mtx.gz    genes x cells, integer, the raw counts of that channel
  input/<sample>/barcodes.tsv     the cell ids of that channel, in the object's order
  input/<sample>/clusters.tsv     one cluster label per cell (see below)
  input/<sample>/soup.tsv         the empty-droplet profile of 06_1, aligned to genes.tsv

**Why gzip, and why per sample.** The counts layer of `shiao.h5ad` has 915,882,723 non-zero
entries. Written out as plain Matrix Market that is ~16 GB of text, and the adjusted matrix
would be another 16 GB, on a laptop with ~48 GB free. Two decisions follow from that number
and they are the only reason this step looks more elaborate than 05_1's equivalent:
gzipping the .mtx (~4x, and `Matrix::readMM(gzfile(...))` reads it directly), and having
`run_soupx.R` write back only the *removals* rather than the whole adjusted matrix - see its
docstring. Peak extra disk is then ~4 GB, and `soupx_all.sh --clean` drops even that once
the assembly is done.

**The clusters.** `autoEstCont()` estimates the contamination fraction by looking, cluster by
cluster, at cells that should not express a gene and do. It therefore needs a grouping at
roughly cell-type resolution, and the object already has one: `optscib_unintegrated_leiden`,
the 32-cluster partition phase 01_5 computed on the unintegrated PCA. It is used here exactly
as `05_2/clustering_tum.py` uses `cell_type_01_4` for its NMI target - as one grouping to
condition an estimate on, never as a biological claim - and it is preferred to `cell_type`
for a reason worth stating: the CellTypist labels are one of the things ambient RNA is
suspected of having distorted, and 06_4 re-runs them on the corrected counts. Estimating the
correction from the labels it is meant to fix would close a loop. `--cluster-key` overrides
it for anyone who wants to check the estimate is not driven by that choice.

Input : $DATA_DIR/shiao.h5ad              (read backed; layers['counts'] = raw counts)
        $DATA_DIR/06_amb/soup/<sample>.csv.gz   from 06_1, one per sample, all required
Output: $DATA_DIR/06_amb/input/...        as above
        $DATA_DIR/06_amb/sample_census.csv      cells, genes and soup UMIs per channel

Nothing under $DATA_DIR/05_tum or $DATA_DIR/04_epi is read or written, and `shiao.h5ad`
is opened read-only: phase 06 is a third branch off the same object and owns 06_amb/ alone.

Local usage (benchmark-py-r), normally through soupx_all.sh:
    export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
    python3 prepare_soupx_input.py
    python3 prepare_soupx_input.py --samples P01_A_P P02_B_P    # a subset
    python3 prepare_soupx_input.py --force                      # rewrite existing
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
import scipy.io as sio
import scipy.sparse as sp

SAMPLE_KEY = "sample"
CLUSTER_KEY = "optscib_unintegrated_leiden"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--samples", nargs="+", default=None,
                   help="only these samples [default: every sample in the object]")
    p.add_argument("--cluster-key", default=CLUSTER_KEY,
                   help=f"obs column used as SoupX clusters [{CLUSTER_KEY}]")
    p.add_argument("--force", action="store_true",
                   help="rewrite input directories that already exist")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    data_dir = Path(os.environ["DATA_DIR"]).expanduser().resolve()
    in_path = data_dir / "shiao.h5ad"
    amb_dir = data_dir / "06_amb"
    soup_dir = amb_dir / "soup"
    input_dir = amb_dir / "input"
    input_dir.mkdir(parents=True, exist_ok=True)

    assert in_path.exists(), f"missing {in_path} (produced by phase 01)"
    if not soup_dir.is_dir():
        print(f"no soup profiles at {soup_dir}; run 06_1 on the cluster and copy them down "
              "(see 06_1_soup_profile/README.md)", file=sys.stderr)
        return 1

    print(f"DATA_DIR : {data_dir}")
    print(f"Input    : {in_path} ({in_path.stat().st_size / 1024**3:.2f} GB)")
    print(f"Soup     : {soup_dir}")
    print(f"Output   : {input_dir}")
    print(f"Clusters : {args.cluster_key}")
    print()

    # Backed: only .obs and one channel's slice of the counts are ever needed at once.
    adata = ad.read_h5ad(in_path, backed="r")
    assert "counts" in adata.layers, "expected raw counts in layers['counts']"
    assert args.cluster_key in adata.obs, \
        f"obs has no '{args.cluster_key}'; pass --cluster-key with a column that exists"

    var_names = pd.Index(adata.var_names)
    genes_path = input_dir / "genes.tsv"
    if args.force or not genes_path.exists():
        pd.Series(var_names).to_csv(genes_path, index=False, header=False)
    else:
        # A shared genes.tsv is only safe if it still describes this object. It is the
        # row space of every counts.mtx.gz and of every soup.tsv, so a mismatch would be
        # a silent gene-shift in the correction rather than an error.
        on_disk = pd.read_csv(genes_path, header=None)[0].to_numpy()
        assert np.array_equal(on_disk, var_names.to_numpy()), (
            f"{genes_path} does not match the genes of {in_path.name}; "
            "re-run with --force"
        )
    print(f"genes: {len(var_names):,} -> {genes_path}")

    samples = (list(adata.obs[SAMPLE_KEY].cat.categories) if args.samples is None
               else list(args.samples))
    clusters_all = adata.obs[args.cluster_key].astype(str)
    sample_all = adata.obs[SAMPLE_KEY].astype(str)

    # Every channel needs its own soup: a missing profile cannot be substituted by another
    # sample's, so it is a hard stop before any work is done rather than 100 samples in.
    missing = [s for s in samples if not (soup_dir / f"{s}.csv.gz").exists()]
    if missing:
        print(f"{len(missing)} sample(s) without a soup profile in {soup_dir}:",
              file=sys.stderr)
        print("  " + ", ".join(missing[:20]) + (" ..." if len(missing) > 20 else ""),
              file=sys.stderr)
        print("Re-run 06_1 for them, or pass --samples to prepare only the ones you have.",
              file=sys.stderr)
        return 1

    rows = []
    for i, sample in enumerate(samples, 1):
        out_dir = input_dir / sample
        done_marker = out_dir / "counts.mtx.gz"
        if done_marker.exists() and not args.force:
            n_cells = sum(1 for _ in open(out_dir / "barcodes.tsv"))
            print(f"[have] ({i}/{len(samples)}) {sample}: input already written, skipping")
            rows.append(dict(sample=sample, status="prepared", n_cells=n_cells))
            continue

        in_sample = (sample_all == sample).to_numpy()
        n_cells = int(in_sample.sum())
        if n_cells == 0:
            print(f"[skip] ({i}/{len(samples)}) {sample}: no cell in the object",
                  file=sys.stderr)
            rows.append(dict(sample=sample, status="no_cells", n_cells=0))
            continue

        # Increasing positions: `.layers['counts']` on a backed object requires it, and the
        # object's own cell order is what barcodes.tsv records, so no reordering happens.
        pos = np.flatnonzero(in_sample)
        counts = sp.csr_matrix(adata.layers["counts"][pos, :])
        cell_ids = adata.obs_names[pos]
        cell_clusters = clusters_all.iloc[pos].to_numpy()

        # The soup profile of 06_1 lives on the Cell Ranger reference (33k genes); the object
        # is the 30,869 that survived phase 01. Reindex onto the object's genes - a gene the
        # object dropped simply leaves the soup - and renormalise `est`, because SoupX reads
        # it as a distribution over the genes it is given.
        soup = pd.read_csv(soup_dir / f"{sample}.csv.gz")
        soup = soup.drop_duplicates("gene").set_index("gene").reindex(var_names)
        soup_counts = soup["counts"].fillna(0.0).to_numpy()
        soup_total = float(soup_counts.sum())
        assert soup_total > 0, f"{sample}: the soup profile is empty on this gene set"
        n_missing_genes = int(soup["counts"].isna().sum())
        soup_out = pd.DataFrame({
            "gene": var_names, "counts": soup_counts, "est": soup_counts / soup_total,
        })

        out_dir.mkdir(parents=True, exist_ok=True)
        # SoupX wants genes as ROWS and cells as COLUMNS, like every 10x matrix.
        with gzip.open(out_dir / "counts.mtx.gz", "wb") as fh:
            sio.mmwrite(fh, counts.T.astype(np.int32), field="integer")
        pd.Series(cell_ids).to_csv(out_dir / "barcodes.tsv", index=False, header=False)
        pd.Series(cell_clusters).to_csv(out_dir / "clusters.tsv", index=False, header=False)
        soup_out.to_csv(out_dir / "soup.tsv", sep="\t", index=False)

        mtx_mb = (out_dir / "counts.mtx.gz").stat().st_size / 1024**2
        n_clusters = len(set(cell_clusters))
        print(f"[ok]   ({i}/{len(samples)}) {sample}: {n_cells:,} cells, "
              f"{counts.nnz:,} non-zero ({mtx_mb:.0f} MB gz), {n_clusters} clusters, "
              f"soup {soup_total:,.0f} UMIs", flush=True)

        rows.append(dict(
            sample=sample, status="prepared", n_cells=n_cells, n_nonzero=int(counts.nnz),
            n_clusters=n_clusters, cell_umis=int(counts.sum()), soup_umis=int(soup_total),
            soup_genes_not_in_object=n_missing_genes,
        ))
        del counts

    census = pd.DataFrame(rows)
    census_path = amb_dir / "sample_census.csv"
    if census_path.exists() and args.samples is not None:
        # A partial run must not truncate the census the R stage reads to decide what to run.
        old = pd.read_csv(census_path)
        census = pd.concat([old[~old["sample"].isin(census["sample"])], census],
                           ignore_index=True).sort_values("sample")
    census.to_csv(census_path, index=False)

    print()
    print(census["status"].value_counts().to_string())
    print(f"\ncensus -> {census_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
