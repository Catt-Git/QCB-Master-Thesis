#!/usr/bin/env python3
"""04_4: CytoTRACE2 potency per cell, run one patient at a time.

CytoTRACE2 is a per-cell predictor of differentiation potency, not a gene set. It belongs
to Route A only and has no Route B counterpart. Its value here is that it is INDEPENDENT
of the lab's lists, so it is the only non-circular evidence on the stemness axis - which
matters more than usual in this stage, because the collaborator's stemness consensus was
deliberately not built (it captured proliferation only) and no single stemness signature
is primary.

Why per patient. The score is computed within the object it is given, so running all 29
cohorts at once lets the ranking absorb batch: a patient sequenced deeper would come out
more potent than one sequenced shallow, and "stemness" would be a batch axis. Each cohort
is scored on its own and the results are concatenated afterwards, so the output cannot be
confounded with `cohort` by construction. The cost is that scores are comparable WITHIN a
patient and only ordinally across patients - which is exactly how Route A uses them, since
every readout there is standardised within (cohort, cell_type) anyway.

API, read from the package documentation rather than assumed
(https://github.com/digitalcytometry/cytotrace2, cytotrace2_python):

    cytotrace2(input_path, annotation_path=None, species="mouse", batch_size=10000,
               smooth_batch_size=1000, disable_parallelization=False,
               disable_plotting=False, max_cores=None, seed=14,
               output_dir="cytotrace2_results")

  * input is a TAB-DELIMITED .txt, GENES (rows) x CELLS (columns), first row cell ids and
    first column gene names - not an .h5ad, which the package does not read;
  * counts must be raw or CPM/TPM and NOT log-transformed, so this reads
    `layers['counts']` and never `.X`, which is scran log1p from 04_1 onward;
  * `species` defaults to "mouse" and MUST be set to "human" here. Human data is handled
    by orthology mapping onto the mouse feature set;
  * output columns: CytoTRACE2_Score, CytoTRACE2_Potency, CytoTRACE2_Relative,
    preKNN_CytoTRACE2_Score, preKNN_CytoTRACE2_Potency. Potency categories are
    Differentiated / Unipotent / Oligopotent / Multipotent / Pluripotent / Totipotent.

The package needs an environment of its own, and this is not a preference. cytotrace2-py
1.1.0.4 declares `numpy<2.0.0` as a hard pin, so `pip install cytotrace2-py` into
`benchmark-py-r` does not add a package - it silently rolls the whole stack back
(numpy 2.4 -> 1.26, pandas 3.0 -> 2.3, scanpy 1.12 -> 1.11, anndata 0.13 -> 0.12, plus
scipy and zarr) and leaves fast-array-utils and tifffile with unsatisfiable requirements.
That was done once and undone; the two stacks cannot coexist. Install it into a dedicated
env instead (`environments/cytotrace2-py.yml`, or the two commands in the README), which is
cheap because this step touches nothing else: it reads one .h5ad, writes .txt matrices, and
writes one .csv. This script checks for the package and stops with that instruction rather
than failing halfway. Everything else in 04_3 runs without it; cell_first_epi.py simply drops
the CytoTRACE2 quadrant definition if the output is missing.

Usage (in the cytotrace2-py env, NOT in benchmark-py-r):
    conda activate cytotrace2-py
    export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
    python cytotrace2_epi.py --dry-run     # export the per-patient matrices, do not score
    python cytotrace2_epi.py               # export + score + concatenate
    python cytotrace2_epi.py --cores 4     # cap the cores (see the memory note in the README)
    python cytotrace2_epi.py --write-obs   # also add the columns to shiao_epi_<run_id>.h5ad
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

# signature_common lives in the phase's utils/, as in 02_2_integration.
UTILS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils")
sys.path.insert(0, UTILS_DIR)
import signature_common as C  # noqa: E402
import sig_collections as SC  # noqa: E402

CT2_SEED = 14          # the package default, pinned here so it is on the record
CT2_SPECIES = "human"  # NOT the package default, which is "mouse"

WORK_DIR = C.EPI_DIR / f"cytotrace2_work_{C.RUN_ID}"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dry-run", action="store_true",
                   help="write the per-patient input matrices and stop (no package needed)")
    p.add_argument("--cores", type=int, default=None,
                   help="max_cores passed to cytotrace2; use 1-2 under 16 GB of memory")
    p.add_argument("--batch-size", type=int, default=10000)
    p.add_argument("--smooth-batch-size", type=int, default=1000)
    p.add_argument("--force", action="store_true", help="re-run patients already done")
    p.add_argument("--write-obs", action="store_true",
                   help="also write the columns into shiao_epi_<run_id>.h5ad (off by default: "
                        "the .csv is the canonical output and nothing already produced is rewritten)")
    p.add_argument("--keep-work", action="store_true", help="do not delete the per-patient .txt matrices")
    return p.parse_args()


def export_patient_matrix(adata, cells, path: Path) -> None:
    """Genes x cells, tab-delimited, raw counts. Written in gene chunks to stay in memory."""
    sub = adata[cells].to_memory()
    X = sub.layers["counts"]
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    # to genes x cells, integers: the package reads counts, and a '3.0' costs 2 bytes a value
    df = pd.DataFrame(np.rint(X).astype(np.int32).T, index=sub.var_names, columns=sub.obs_names)
    df.index.name = "gene"
    df.to_csv(path, sep="\t")


def main():
    args = parse_args()
    C.banner("04_4 - CytoTRACE2 potency, per patient")

    cytotrace2 = None
    if not args.dry_run:
        try:
            from cytotrace2_py.cytotrace2_py import cytotrace2  # noqa: F401
        except ImportError:
            try:
                from cytotrace2_py import cytotrace2  # noqa: F401
            except ImportError:
                sys.exit(
                    "[STOP] cytotrace2 is not installed in this environment.\n"
                    "       It pins numpy<2 and must NOT be installed into benchmark-py-r,\n"
                    "       which it would roll back to numpy 1.26 / pandas 2 / scanpy 1.11.\n"
                    "       Use the dedicated env:\n"
                    "           conda env create -f environments/cytotrace2-py.yml\n"
                    "           conda activate cytotrace2-py\n"
                    "       (cytotrace2-py 1.1.0.4 on PyPI; see also "
                    "https://github.com/digitalcytometry/cytotrace2)\n"
                    "       Re-run with --dry-run to export the per-patient matrices without it."
                )

    adata = ad.read_h5ad(C.FULL_H5AD, backed="r")
    print(f"{adata.n_obs:,} cells x {adata.n_vars:,} genes")
    assert "counts" in adata.layers, "layers['counts'] missing: CytoTRACE2 must not see log data"

    cohorts = adata.obs["cohort"].cat.categories.tolist()
    print(f"{len(cohorts)} cohorts, scored one at a time\n")

    WORK_DIR.mkdir(parents=True, exist_ok=True)
    results = []

    for i, coh in enumerate(cohorts, start=1):
        cells = adata.obs_names[adata.obs["cohort"] == coh]
        pdir = WORK_DIR / coh
        pdir.mkdir(exist_ok=True)
        mat = pdir / f"counts_{coh}.txt"
        out_csv = pdir / "cytotrace2_results.csv"

        if out_csv.exists() and not args.force:
            print(f"  [{i:>2}/{len(cohorts)}] {coh:<12} [have] {out_csv.name}")
            results.append(pd.read_csv(out_csv, index_col=0))
            continue

        if not mat.exists() or args.force:
            t0 = time.time()
            export_patient_matrix(adata, cells, mat)
            print(f"  [{i:>2}/{len(cohorts)}] {coh:<12} {len(cells):>6,} cells "
                  f"-> {mat.name} ({mat.stat().st_size / 1e6:.0f} MB, {time.time() - t0:.0f}s)")
        else:
            print(f"  [{i:>2}/{len(cohorts)}] {coh:<12} {len(cells):>6,} cells [have] {mat.name}")

        if args.dry_run:
            continue

        t0 = time.time()
        df = cytotrace2(
            str(mat),
            species=CT2_SPECIES,
            batch_size=args.batch_size,
            smooth_batch_size=args.smooth_batch_size,
            max_cores=args.cores,
            seed=CT2_SEED,
            disable_plotting=True,
            output_dir=str(pdir),
        )
        if df is None or not isinstance(df, pd.DataFrame):
            # older builds only write the file
            found = list(pdir.glob("*cytotrace2*result*.csv")) + list(pdir.glob("*.csv"))
            if not found:
                sys.exit(f"[STOP] {coh}: cytotrace2 returned nothing and wrote no csv into {pdir}")
            df = pd.read_csv(sorted(found)[0], index_col=0)

        df.to_csv(out_csv)
        results.append(df)
        print(f"       {coh:<12} scored in {time.time() - t0:.0f}s -> {out_csv.name}")

    if args.dry_run:
        print(f"\ndry run: matrices in {WORK_DIR}, nothing scored.")
        return

    # ------------------------------------------------------------- concatenate
    ct2 = pd.concat(results, axis=0)
    ct2.index.name = "cell"

    missing = set(adata.obs_names) - set(ct2.index)
    extra = set(ct2.index) - set(adata.obs_names)
    print(f"\n{len(ct2):,} cells scored; {len(missing)} unscored, {len(extra)} unknown ids")
    if extra:
        ct2 = ct2.drop(index=list(extra))

    # Reindex onto the object's order - never assume the package preserved it.
    ct2 = ct2.reindex(adata.obs_names)
    ct2.to_csv(C.CYTOTRACE_CSV)
    print(f"[write] {C.CYTOTRACE_CSV}")

    if "CytoTRACE2_Potency" in ct2:
        print("\npotency composition")
        print(ct2["CytoTRACE2_Potency"].value_counts(dropna=False).to_string())

    if "CytoTRACE2_Score" in ct2:
        by_ct = (pd.DataFrame({"score": ct2["CytoTRACE2_Score"].values,
                               "cell_type": adata.obs["cell_type"].values})
                 .groupby("cell_type", observed=True)["score"].agg(["mean", "median", "size"]))
        print("\nmean potency per cell type (the sanity read: Lumsec-prol should NOT be the "
              "only thing at the top - if it is, potency is tracking the cycle)")
        print(by_ct.sort_values("mean", ascending=False).to_string(float_format="%.3f"))
        # CytoTRACE2 is the stemness readout of the scie collection - the only one in the
        # stage that does not come from a gene list - so its summary table belongs there.
        # The per-cell csv above is NOT collection-scoped: it is a measurement, and 04_5
        # joins it only when the collection asks for it.
        C.write_table(by_ct, "cytotrace2_by_cell_type", SC.SCIE)

    if args.write_obs:
        target = C.EPI_DIR / f"shiao_epi_{C.RUN_ID}.h5ad"
        print(f"\n[write-obs] adding the CytoTRACE2 columns to {target}")
        full = ad.read_h5ad(target)
        assert (full.obs_names == ct2.index).all(), "cell order differs; refusing to write"
        for col in ct2.columns:
            full.obs[f"ct2_{col}"] = ct2[col].values
        full.write_h5ad(target)
        print("[write-obs] done")

    if not args.keep_work:
        for coh in cohorts:
            for f in (WORK_DIR / coh).glob("counts_*.txt"):
                f.unlink()
        print(f"\nper-patient count matrices deleted (results kept in {WORK_DIR}); "
              "--keep-work keeps them")

    print("\ndone.")


if __name__ == "__main__":
    main()
