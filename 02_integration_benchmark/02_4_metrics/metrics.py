#!/usr/bin/env python
"""
02_4 metrics - the 13 scib metrics minus kBET, one job per (run, output type)

Scores one integrated object against its matching reference and writes a single
CSV whose one column is named ``<method>_<type>``. Twelve metrics are computed
here; kBET is scored separately by metrics_kbet.py, which patches its value into
this same CSV. Everything the two jobs share (cell-order enforcement, reference
PCA, the reduce_data options per type) lives in _metrics_shared.py.

Notes:
  * ``scib_compat`` is imported before scib (via _metrics_shared), so the job does
    not die partway through the computation on a numpy/pandas API scib 1.1.7 no
    longer finds.
  * kBET is off here (kBET_=False); it is the slowest metric and gets its own job.
  * ``reduce_data(umap=False)``: no metric reads X_umap.
  * ``trajectory_=False`` always: 13 metrics, not 14.
  * ``--hvgs 0`` -> no HVG re-selection: the reference is already HVG-restricted,
    and re-selecting on scaled (negative) values is fatal.
  * cell order is a hard assertion, not a silent obs_names rename.

Usage:
    python metrics.py -u <reference.h5ad> -i <integrated.h5ad> -o <out.csv> \
        -m harmony --type embed [-b cohort -l cell_type --organism human --hvgs 0 -v]

The reference passed with -u MUST match the scaling variant of the integrated
object; the grid's ``reference`` column already tracks that.
"""

import argparse
import os
import sys

# metrics_shared lives in utils/ and imports scib_compat before scib; importing
# it first is enough for the caller too.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils"))
import metrics_shared as shared  # noqa: E402
import scib  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description="Compute the 13 scib metrics minus kBET")
    p.add_argument("-u", "--uncorrected", required=True, help="reference .h5ad (matching scaling)")
    p.add_argument("-i", "--integrated", required=True, help="integrated .h5ad")
    p.add_argument("-o", "--output", required=True, help="output CSV")
    p.add_argument("-m", "--method", required=True, help="method name, for the column <method>_<type>")
    p.add_argument("--type", required=True, choices=shared.RESULT_TYPES, help="scib output type")
    p.add_argument("-b", "--batch_key", default="cohort")
    p.add_argument("-l", "--label_key", default="cell_type")
    p.add_argument("--organism", default="human")
    p.add_argument("--hvgs", type=int, default=0,
                   help="HVGs for scib; 0 means no re-selection (the reference is already HVG-restricted)")
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args()


def main():
    args = parse_args()
    type_ = args.type
    setup = f"{args.method}_{type_}"
    n_hvgs = args.hvgs if args.hvgs > 0 else None

    if args.verbose:
        print("=" * 70, flush=True)
        print(f"METRICS  {setup}", flush=True)
        print(f"  reference : {args.uncorrected}", flush=True)
        print(f"  integrated: {args.integrated}", flush=True)
        print(f"  batch/label: {args.batch_key} / {args.label_key}", flush=True)
        print(f"  organism  : {args.organism}   hvgs: {n_hvgs}", flush=True)
        print("=" * 70, flush=True)

    adata, adata_int = shared.load_and_align(
        args.uncorrected, args.integrated, args.batch_key, args.label_key, verbose=args.verbose
    )
    shared.ensure_reference_reduced(adata, verbose=args.verbose)
    shared.reduce_integrated(adata_int, type_, n_hvgs=n_hvgs, verbose=args.verbose)

    flags = shared.metric_flags(type_)
    embed = shared.embed_key(type_)

    # The optimal-resolution clustering NMI/ARI need is written next to the CSV.
    file_stump = os.path.splitext(args.output)[0]
    cluster_nmi = f"{file_stump}_nmi.txt"

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

    print(f"computing metrics for {setup} (kBET in a separate job)", flush=True)
    results = scib.me.metrics(
        adata,
        adata_int,
        batch_key=args.batch_key,
        label_key=args.label_key,
        embed=embed,
        type_=type_,
        organism=args.organism,
        cluster_nmi=cluster_nmi,
        nmi_method="arithmetic",
        nmi_dir=None,
        n_isolated=None,
        verbose=args.verbose,
        **flags,
    )
    results.rename(columns={results.columns[0]: setup}, inplace=True)
    results.to_csv(args.output)

    if args.verbose:
        print(results, flush=True)
    print(f"done -> {args.output}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
