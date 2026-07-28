#!/usr/bin/env python
"""
02_4 metrics - merge the per-run CSVs into the single table the plotting reads

Each metrics job writes one CSV with a single column of metric values, at
``<root>/<task>/metrics/<scaling>/<hvg>/<method>_<type>.csv``. This collects them
into one table whose row index is the path of each CSV relative to ``--root`` and
whose columns are the metric names -- the exact shape ``plotSingleTaskRNA.R``
wants, since it derives task / scaling / feature space / method by splitting
that first column on ``/``.

Every CSV shares the same metric index (scib.me.metrics always returns the full
set, NaN for the metrics that do not apply to a type), so the inner merge is a
full join and no metric is dropped.

Usage:
    python merge_metrics.py -o <merged.csv> -r <root/> -i <csv> [<csv> ...]
    # or let it glob the tree itself:
    python merge_metrics.py -o <merged.csv> -r <root/> --glob '<root>/**/hvg/*.csv'
"""

import argparse
import glob as globmod
from functools import reduce

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Merge per-run metric CSVs into one table")
    p.add_argument("-o", "--output", required=True, help="merged CSV")
    p.add_argument("-r", "--root", required=True,
                   help="root the column names are made relative to (e.g. .../02_metrics/)")
    p.add_argument("-i", "--input", nargs="*", default=[], help="the per-run CSVs")
    p.add_argument("--glob", default=None, help="glob pattern to collect the CSVs instead of listing them")
    return p.parse_args()


def main():
    args = parse_args()

    files = list(args.input)
    if args.glob:
        files += sorted(globmod.glob(args.glob, recursive=True))
    files = sorted(set(files))
    if not files:
        raise SystemExit("no input CSVs given (use -i or --glob)")

    root = args.root if args.root.endswith("/") else args.root + "/"

    res_list = []
    for file in files:
        clean = file.replace(root, "").replace(".csv", "").lstrip("/")
        res = pd.read_csv(file, index_col=0)
        res.rename(columns={res.columns[0]: clean}, inplace=True)
        res_list.append(res)

    results = reduce(
        lambda left, right: pd.merge(left, right, left_index=True, right_index=True),
        res_list,
    )
    results = results.T
    results.to_csv(args.output)
    print(f"merged {len(files)} CSV(s) -> {args.output}  ({results.shape[0]} runs x {results.shape[1]} metrics)",
          flush=True)


if __name__ == "__main__":
    main()
