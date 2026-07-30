#!/usr/bin/env python
"""
02_4 metrics - kBET alone, in its own job

kBET is by far the slowest metric and the one that times out; scoring it in a
separate job means a timeout costs one metric instead of thirteen and the main
metrics.py CSV is already complete without it. This job computes kBET with the
same reference/alignment/reduce_data path as metrics.py (via _metrics_shared) and
then writes its value into the ``kBET`` row of the CSV metrics.py produced, so the
final table still has one CSV per (run, type).

Usage (run AFTER metrics.py has written --output):
    python metrics_kbet.py -u <reference.h5ad> -i <integrated.h5ad> -o <same out.csv> \
        -m harmony --type embed [-b cohort -l cell_type --max-cells 0 -v]
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils"))
import metrics_shared as shared  # noqa: E402
import scib  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description="Compute kBET only and patch it into the metrics CSV")
    p.add_argument("-u", "--uncorrected", required=True, help="reference .h5ad (matching scaling)")
    p.add_argument("-i", "--integrated", required=True, help="integrated .h5ad")
    p.add_argument("-o", "--output", required=True,
                   help="the CSV metrics.py wrote; the kBET row is patched into it")
    p.add_argument("-m", "--method", required=True)
    p.add_argument("--type", required=True, choices=shared.RESULT_TYPES)
    p.add_argument("-b", "--batch_key", default="cohort")
    p.add_argument("-l", "--label_key", default="cell_type")
    p.add_argument("--max-cells", type=int, default=0,
                   help="skip kBET above this many cells (0 = no cap, the default)")
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args()


def patch_csv(path, setup, value):
    """Write ``value`` into the kBET row of the metrics CSV.

    If the CSV metrics.py should have produced is missing, a standalone one is
    written with the full metric index and only kBET filled, so the job is not
    lost - but that is a warning: the two CSVs are meant to be one file.
    """
    if os.path.exists(path):
        df = pd.read_csv(path, index_col=0)
        col = df.columns[0]
        if "kBET" not in df.index:
            df.loc["kBET"] = np.nan
        df.loc["kBET", col] = value
        df.to_csv(path)
        return "patched"

    print(f"[warn] {path} does not exist yet; run metrics.py first. Writing a "
          "kBET-only CSV as a fallback.", flush=True)
    df = pd.DataFrame({setup: [np.nan] * len(shared.ALL_METRICS)}, index=list(shared.ALL_METRICS))
    df.loc["kBET", setup] = value
    df.to_csv(path)
    return "created"


def main():
    args = parse_args()
    type_ = args.type
    setup = f"{args.method}_{type_}"

    adata, adata_int = shared.load_and_align(
        args.uncorrected, args.integrated, args.batch_key, args.label_key, verbose=args.verbose
    )

    if args.max_cells and adata.n_obs > args.max_cells:
        print(f"[skip] {setup}: {adata.n_obs} cells > --max-cells {args.max_cells}; "
              "kBET not computed, leaving the row as NaN.", flush=True)
        # Touch nothing if the CSV already exists (kBET stays NaN there).
        if not os.path.exists(args.output):
            patch_csv(args.output, setup, np.nan)
        return 0

    shared.ensure_reference_reduced(adata, verbose=args.verbose)
    shared.reduce_integrated(adata_int, type_, n_hvgs=None, verbose=args.verbose)

    embed = shared.embed_key(type_)
    print(f"computing kBET for {setup}", flush=True)

    # Only kBET is enabled; everything else stays off so scib returns the fixed
    # metric index with kBET filled and the rest NaN.
    results = scib.me.metrics(
        adata,
        adata_int,
        batch_key=args.batch_key,
        label_key=args.label_key,
        embed=embed,
        type_=type_,
        kBET_=True,
        verbose=args.verbose,
    )
    value = float(results.loc["kBET"].iloc[0])
    print(f"  kBET = {value}", flush=True)

    action = patch_csv(args.output, setup, value)
    print(f"done -> {args.output} ({action})", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
