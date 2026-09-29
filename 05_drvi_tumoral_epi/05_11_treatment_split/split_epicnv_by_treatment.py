#!/usr/bin/env python3
"""05_11: one .h5ad per treatment arm, cut from the post-CNV epithelial object.

Input is `shiao_epicnv_raw.h5ad` (05_2, CELL_SET=epi): every epithelial cell under the post-CNV
labels, malignant, non_malignant and borderline alike, raw counts in `.X` and `layers['counts']`.
The split is a pure row subset on `treatment`: no gene filter, no normalization, no QC is
recomputed. The only change is dropping the categories an arm does not contain.

Output: $DATA_DIR/05_tum/shiao_epicnv_raw_{BASE,PD1,RTPD1}.h5ad

Usage:
    export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
    python split_epicnv_by_treatment.py
"""
from __future__ import annotations

import os

import anndata as ad

DATA_DIR = os.environ.get("DATA_DIR", os.path.expanduser("~/Desktop/QCB-Master-Thesis/datasets"))
IN_PATH = os.path.join(DATA_DIR, "05_tum", "shiao_epicnv_raw.h5ad")
OUT_TEMPLATE = os.path.join(DATA_DIR, "05_tum", "shiao_epicnv_raw_{}.h5ad")
TREATMENTS = ["BASE", "PD1", "RTPD1"]


def main() -> None:
    adata = ad.read_h5ad(IN_PATH)
    print(f"input: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    n_written = 0
    for treatment in TREATMENTS:
        sub = adata[adata.obs["treatment"] == treatment].copy()
        for col in sub.obs.select_dtypes("category"):
            sub.obs[col] = sub.obs[col].cat.remove_unused_categories()
        out = OUT_TEMPLATE.format(treatment)
        sub.write_h5ad(out, compression="gzip")
        n_written += sub.n_obs
        status = sub.obs["cnv_status"].value_counts().to_dict()
        print(f"{treatment}: {sub.n_obs:,} cells, {sub.obs['cohort'].nunique()} cohorts, {status} -> {out}")

    assert n_written == adata.n_obs, "the three arms do not cover the input"


if __name__ == "__main__":
    main()
