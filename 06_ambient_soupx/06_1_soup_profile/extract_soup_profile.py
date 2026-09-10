#!/usr/bin/env python3
"""06_1: the ambient expression profile of each sample, measured in its empty droplets.

SoupX corrects a cell's counts against a *profile of the soup* - the relative abundance of
each gene in the ambient RNA floating in the emulsion. That profile is not a quantity you
can read off the cells: it has to be measured where there are no cells, i.e. in the droplets
Cell Ranger called EMPTY. This is the one step of the phase that cannot run from
`shiao.h5ad`, because `shiao.h5ad` was built from `filtered_feature_bc_matrix.h5` (see
`00_6_build_h5ad/build_combined_h5ad.py`) and the empty droplets were dropped there and then.

So this script goes back to the Cell Ranger output, reads the UNFILTERED matrix, and reduces
each sample's several million droplets to one table of 30-33k numbers. That reduction is the
whole point of putting the step on the cluster: the raw matrices are tens of GB and live next
to the FASTQ, the profiles are a few hundred KB each and are what comes back down.

  $CELLRANGER_DIR/<sample>/outs/raw_feature_bc_matrix.h5   ~0.5-2 GB, on the cluster
                              |
                              v
  $DATA_DIR/06_amb/soup/<sample>.csv.gz                    ~300 KB, copied to the laptop

Which droplets count as empty. SoupX's own `estimateSoup()` sums the droplets whose total
UMI count falls in `soupRange = c(0, 100)`, and this reproduces that rule rather than
inventing one: droplets with 1-`MAX_UMI` counts. The lower bound is exclusive because a
droplet with zero counts contributes nothing but would inflate the droplet count in the
census. The upper bound is the assumption of the method - that nothing above ~100 UMIs is
reliably cell-free - and it is the one parameter here worth varying if a sample looks odd.

Note what is NOT done: the empty droplets are not intersected with "the barcodes Cell Ranger
did not call". They do not need to be. Any barcode Cell Ranger called as a cell has far more
than 100 UMIs, so the UMI ceiling excludes the called cells by construction, and the rule
stays the same whether or not the filtered matrix is at hand.

Gene names are made unique exactly as `build_combined_h5ad.py` did (`var_names_make_unique()`
on the same GRCh38-3.0.0 reference), so the symbols in the profile and the symbols in
`shiao.h5ad` are the same vocabulary. They are NOT the same length - phase 01 filtered genes -
and the alignment of the two is done downstream, in `06_2/prepare_soupx_input.py`, where the
object's gene order is known.

Input : $CELLRANGER_DIR/<sample>/outs/raw_feature_bc_matrix.h5
Output: $DATA_DIR/06_amb/soup/<sample>.csv.gz    gene, counts, est   (est sums to 1)
        $DATA_DIR/06_amb/soup_census.csv         one row per sample, what went into it

Cluster usage (the way this was run; catalano_env, ~64 GB for the biggest raw matrix):
    cd 06_ambient_soupx/06_1_soup_profile && mkdir -p logs
    export DATA_DIR=/users/genomics/albertoc/Tesi/hopes_and_dreams/datasets
    sbatch --export=ALL,DATA_DIR=$DATA_DIR submit_soup_profile.slurm

Local usage (only if the raw matrices are on this machine, which they are not by default):
    export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
    python3 extract_soup_profile.py --cellranger-dir /path/to/cellranger_out
    python3 extract_soup_profile.py --samples P01_A_P P02_B_P     # a subset
    python3 extract_soup_profile.py --force                       # rewrite existing
"""

from __future__ import annotations

import argparse
import gzip
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

# SoupX's `soupRange` upper bound, verbatim. A droplet with more counts than this is not
# assumed to be cell-free, so it does not contribute to the ambient profile.
MAX_UMI = 100

RAW_H5 = "outs/raw_feature_bc_matrix.h5"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--cellranger-dir", default=None,
                   help="directory holding <sample>/outs/raw_feature_bc_matrix.h5 "
                        "[default: $CELLRANGER_DIR, else $DATA_DIR/cellranger_out - which on "
                        "the cluster used here is WRONG, the batch landed beside datasets/ "
                        "rather than inside it]")
    p.add_argument("--out-dir", default=None,
                   help="where the profiles go [default: $DATA_DIR/06_amb/soup]")
    p.add_argument("--samples", nargs="+", default=None,
                   help="only these samples [default: every sample with a raw matrix]")
    p.add_argument("--samples-file", default=None,
                   help="a file with one sample id per line (e.g. renamed_fastq/unique_samples.txt)")
    p.add_argument("--max-umi", type=int, default=MAX_UMI,
                   help=f"a droplet with at most this many UMIs is soup [{MAX_UMI}]")
    p.add_argument("--force", action="store_true",
                   help="recompute profiles that are already on disk")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    data_dir = Path(os.environ["DATA_DIR"]).expanduser().resolve()
    cr_dir = Path(args.cellranger_dir or os.environ.get("CELLRANGER_DIR")
                  or data_dir / "cellranger_out").expanduser().resolve()
    out_dir = Path(args.out_dir or data_dir / "06_amb" / "soup").expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    if not cr_dir.is_dir():
        print(f"no Cell Ranger output at {cr_dir}\n"
              "This step needs the UNFILTERED matrices, which are not in the repo's\n"
              "datasets/. Pass --cellranger-dir or export CELLRANGER_DIR.\n"
              "\n"
              "The default guessed here is $DATA_DIR/cellranger_out, which is what\n"
              "00_5/run_cellranger_batch.sh documents as its OUTDIR. On the machine this\n"
              "was run on, the batch actually landed one level up, as a SIBLING of\n"
              f"datasets/ - try {cr_dir.parent.parent / 'cellranger_out'} before searching\n"
              "further. See the phase README.", file=sys.stderr)
        return 1

    # Which samples. The listing is driven by what actually has a raw matrix on disk, so a
    # half-finished Cell Ranger batch is reported as missing rather than crashing later.
    if args.samples:
        samples = list(args.samples)
    elif args.samples_file:
        samples = [ln.strip() for ln in open(args.samples_file) if ln.strip()]
    else:
        samples = sorted(p.name for p in cr_dir.iterdir()
                         if (p / RAW_H5).exists())

    print(f"DATA_DIR       : {data_dir}")
    print(f"cellranger dir : {cr_dir}")
    print(f"output dir     : {out_dir}")
    print(f"samples        : {len(samples)}")
    print(f"soup range     : 1 - {args.max_umi} UMIs per droplet")
    print(flush=True)

    rows = []
    for i, sample in enumerate(samples, 1):
        out_path = out_dir / f"{sample}.csv.gz"
        h5_path = cr_dir / sample / RAW_H5

        if out_path.exists() and not args.force:
            # Still censused, so a resumed run does not produce a census with holes in it.
            prof = pd.read_csv(out_path)
            print(f"[have] ({i}/{len(samples)}) {sample}: {out_path.name} exists, skipping")
            rows.append(dict(sample=sample, status="have", n_genes=len(prof),
                             soup_umis=int(prof["counts"].sum())))
            continue

        if not h5_path.exists():
            print(f"[MISS] ({i}/{len(samples)}) {sample}: no {RAW_H5}", file=sys.stderr)
            rows.append(dict(sample=sample, status="missing_raw_matrix"))
            continue

        raw = sc.read_10x_h5(h5_path)
        raw.var_names_make_unique()

        # Per-droplet totals over the whole unfiltered matrix. `.X` here is droplets x genes.
        totals = np.asarray(raw.X.sum(axis=1)).ravel()
        is_soup = (totals >= 1) & (totals <= args.max_umi)
        n_soup = int(is_soup.sum())

        if n_soup == 0:
            # Possible in principle (a pre-filtered matrix), and it would make `est` all-NaN.
            print(f"[FAIL] ({i}/{len(samples)}) {sample}: no droplet in the soup range",
                  file=sys.stderr)
            rows.append(dict(sample=sample, status="no_empty_droplets",
                             n_droplets=int(raw.n_obs)))
            del raw
            continue

        counts = np.asarray(raw.X[is_soup, :].sum(axis=0)).ravel()
        total_umis = float(counts.sum())
        prof = pd.DataFrame({
            "gene": raw.var_names.to_numpy(),
            "counts": counts,
            # `est` is SoupX's own name for the profile: the fraction of the soup's UMIs
            # belonging to each gene. It is renormalised again downstream, after the
            # intersection with the object's genes drops a few thousand of them.
            "est": counts / total_umis,
        })
        del raw

        with gzip.open(out_path, "wt") as fh:
            prof.to_csv(fh, index=False)

        top = prof.nlargest(3, "counts")["gene"].tolist()
        print(f"[ok]   ({i}/{len(samples)}) {sample}: {n_soup:,} empty droplets, "
              f"{total_umis:,.0f} soup UMIs, top: {', '.join(top)}", flush=True)
        rows.append(dict(sample=sample, status="ok", n_droplets=int(len(totals)),
                         n_soup_droplets=n_soup, soup_umis=int(total_umis),
                         n_genes=int(len(prof)), top_gene=top[0]))

    census = pd.DataFrame(rows)
    census_path = out_dir.parent / "soup_census.csv"
    if census_path.exists() and (args.samples or args.samples_file):
        # A partial run must not truncate the census of the full one: merge on `sample`,
        # the new rows winning.
        old = pd.read_csv(census_path)
        census = pd.concat([old[~old["sample"].isin(census["sample"])], census],
                           ignore_index=True).sort_values("sample")
    census.to_csv(census_path, index=False)

    print()
    print(census["status"].value_counts().to_string())
    print(f"\ncensus -> {census_path}")
    bad = census[census["status"].isin({"missing_raw_matrix", "no_empty_droplets"})]
    if len(bad):
        print(f"\n{len(bad)} sample(s) without a profile: {', '.join(bad['sample'])}",
              file=sys.stderr)
        print("06_2 will refuse to run until every sample of shiao.h5ad has one.",
              file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
