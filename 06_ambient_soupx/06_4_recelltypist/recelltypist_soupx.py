#!/usr/bin/env python3
"""06_4: the ambient decision applied, and CellTypist re-run on the corrected counts.

Last step of the phase. It takes the object 06_2 assembled and the keep/drop table 06_3's
notebook wrote, produces the cell set phase 05 will actually see, and re-annotates it.

  06_amb/shiao_soupx_all_cells.h5ad  +  06_amb/ambient_keep.csv
                              |
                              v
  06_amb/shiao.h5ad     <- the drop-in replacement for datasets/shiao.h5ad

**Why the labels are recomputed and not inherited.** `cell_type` in `shiao.h5ad` was written
by `01_4/celltypist_annotation.py` on the *uncorrected* counts, and ambient RNA is exactly
the kind of signal that moves a CellTypist call: a luminal cell floating in the soup of a
tumour full of T cells picks up CD3/PTPRC reads it never transcribed, and the model has to
put it somewhere. Phase 05 then spends those labels twice - `05_1/prepare_infercnv_input.py`
uses them both to choose which cells are tested for aneuploidy (the epithelium) and to build
inferCNV's diploid baseline (T/NK and myeloid) - so a label distorted by the soup does not
just mislabel a cell, it moves the CNV call of every cell in the cohort. Correcting the
counts and keeping the old labels would have fixed the smaller of the two problems.

This is `01_4`'s procedure verbatim, not a variant of it: the same model file, the same
`majority_voting=True`, and the same temporary CP10K+log1p matrix built from raw counts,
because that is CellTypist's training normalisation and neither phase's `.X` matches it.
The one thing that differs is what goes in - corrected counts, on the cells 06_3 kept - and
that is the whole point.

The over-clustering CellTypist builds for the majority vote is computed on THIS population,
which is why the phase-01 PCA, neighbour graph and UMAP are dropped from the temporary copy
before annotating. Same reasoning as `05_1/recelltypist_nonmalignant.py`.

**What the output object is.** Byte-for-byte compatible with what phase 05 expects of
`shiao.h5ad`: `layers['counts']` are integer counts, `.X` is the same matrix (see
`assemble_soupx.py` on why the phase-01 scran normalisation cannot come along), `cell_type`
is a CellTypist label from the 58 of `Cells_Adult_Breast.pkl`, `cohort`/`sample`/`treatment`/
`response` are untouched. `fraction` is recomputed from the new labels by 01_4's rule, since
it is a FUNCTION of `cell_type` and would otherwise be stale the moment the labels change.

The phase-01 label is kept beside it as `cell_type_precorrection`, and the phase-01 fraction
as `fraction_precorrection`. **Not** `cell_type_01_4`, although that is what phase 05 calls a
CellTypist label that is no longer current: `05_2/subset_and_qc.ipynb` builds its own
`cell_type_01_4` by copying whatever `cell_type` its input object carries, which under the
shadow DATA_DIR is the corrected label written here. Two different things under one name is a
bug waiting to be read the wrong way, so this phase uses the `_precorrection` suffix it
already uses for the QC columns, and leaves `cell_type_01_4` to mean exactly what 05_2 makes
it mean.

Input : $DATA_DIR/06_amb/shiao_soupx_all_cells.h5ad   from 06_2
        $DATA_DIR/06_amb/ambient_keep.csv             from 06_3's notebook
        $DATA_DIR/Cells_Adult_Breast.pkl              the 01_4 model, read-only
Output: $DATA_DIR/06_amb/shiao.h5ad                   the drop-in
        $DATA_DIR/06_amb/cell_annotation_soupx.csv    old label, new label, kept/dropped

Local usage (benchmark-py-r):
    export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
    python3 recelltypist_soupx.py
    python3 recelltypist_soupx.py --no-majority-voting   # fast, but not what 01_4 did
    python3 recelltypist_soupx.py --keep-all             # ignore ambient_keep.csv, drop nothing
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

import celltypist
from celltypist import models

sc.settings.verbosity = 1

LABEL_KEY = "cell_type"
PRIOR_LABEL_KEY = "cell_type_precorrection"
KEEP_KEY = "ambient_keep"

# The lineage sets of 01_4/fraction_reassignment.py, verbatim. Duplicated rather than
# imported, as every phase in this repo duplicates: `fraction` is a function of `cell_type`
# and has to be recomputed here, and the rule is a lookup with no state.
# 'Lymph-*' are LYMPHATIC ENDOTHELIAL subtypes in this atlas's nomenclature, not lymphocytes.
IMMUNE_CELL_TYPES = {
    "CD4-Tem", "CD4-Th", "CD4-Th-like", "CD4-Treg", "CD4-activated", "CD4-naive",
    "CD8-Tem", "CD8-Trm", "CD8-activated", "GD", "NK", "NK-ILCs", "NKT", "T_prol",
    "Macro-IFN", "Macro-lipo", "Macro-m1", "Macro-m1-CCL", "Macro-m2", "Macro-m2-CXCL",
    "Mast", "Mono-classical", "Mono-non-classical", "Neutrophil",
    "cDC1", "cDC2", "mDC", "pDC", "mye-prol",
    "b_naive", "bmem_switched", "bmem_unswitched", "plasma_IgA", "plasma_IgG",
}
NON_IMMUNE_CELL_TYPES = {
    "LummHR-SCGB", "LummHR-active", "LummHR-major",
    "Lumsec-HLA", "Lumsec-KIT", "Lumsec-basal", "Lumsec-lac", "Lumsec-major",
    "Lumsec-myo", "Lumsec-prol", "basal",
    "Fibro-SFRP4", "Fibro-major", "Fibro-matrix", "Fibro-prematrix",
    "Lymph-immune", "Lymph-major", "Lymph-valve1", "Lymph-valve2",
    "Vas-arterial", "Vas-capillary", "Vas-venous", "pericytes", "vsmc",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--no-majority-voting", action="store_true",
                   help="skip CellTypist's over-clustering vote (fast, but not what 01_4 did)")
    p.add_argument("--keep-all", action="store_true",
                   help="do not read ambient_keep.csv: correct, re-annotate, drop no cell")
    p.add_argument("--reuse-annotation", action="store_true",
                   help="do not run CellTypist: read cell_type/celltypist_predicted back from "
                        "an existing cell_annotation_soupx.csv. For RESUMING a run whose "
                        "annotation succeeded and whose write_h5ad did not - the labels are "
                        "the ones that run computed, not a cheaper substitute for them")
    p.add_argument("--no-reannotate", action="store_true",
                   help="apply the ambient call but keep phase 01's cell_type. See the module "
                        "docstring: what this gives up is the labels that SET UP inferCNV in "
                        "05_1, not the ones 05_1 recomputes afterwards")
    p.add_argument("--in-name", default="shiao_soupx_all_cells.h5ad",
                   help="input file under 06_amb/ [%(default)s]")
    p.add_argument("--out-name", default="shiao.h5ad",
                   help="output file under 06_amb/; the default is what makes 06_amb/ usable "
                        "as a DATA_DIR for phase 05 [%(default)s]")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    data_dir = Path(os.environ["DATA_DIR"]).expanduser().resolve()
    amb_dir = data_dir / "06_amb"
    in_path = amb_dir / args.in_name
    keep_path = amb_dir / "ambient_keep.csv"
    model_path = data_dir / "Cells_Adult_Breast.pkl"
    out_path = amb_dir / args.out_name
    annot_path = amb_dir / "cell_annotation_soupx.csv"

    assert in_path.exists(), f"missing {in_path} (produced by 06_2/assemble_soupx.py)"
    if args.reuse_annotation and args.no_reannotate:
        print("--reuse-annotation and --no-reannotate are mutually exclusive", file=sys.stderr)
        return 1
    if not args.no_reannotate and not args.reuse_annotation:
        assert model_path.exists(), f"missing {model_path} (the model 01_4 used)"
    if args.reuse_annotation:
        assert annot_path.exists(), (
            f"missing {annot_path}: --reuse-annotation reads the labels back from the table a "
            "previous run of this script wrote; with no such table there is nothing to resume"
        )
    if not args.keep_all:
        assert keep_path.exists(), (
            f"missing {keep_path} (written by 06_3_ambient_qc/ambient_qc.ipynb); "
            "run the notebook, or pass --keep-all to drop no cell"
        )
    if out_path.resolve() == (data_dir / "shiao.h5ad").resolve():
        # The one accident this script could cause that nothing else could undo.
        print(f"refusing to write over the phase-01 object at {out_path}", file=sys.stderr)
        return 1

    print(f"DATA_DIR: {data_dir}")
    print(f"Input   : {in_path}")
    if args.no_reannotate:
        labels_note = "kept from phase 01 (--no-reannotate)"
    elif args.reuse_annotation:
        labels_note = f"read back from {annot_path.name} (--reuse-annotation)"
    else:
        labels_note = "recomputed with CellTypist on the corrected counts"
    print(f"Labels  : {labels_note}")
    print(f"Keep    : {keep_path if not args.keep_all else '(--keep-all, nothing dropped)'}")
    print(f"Output  : {out_path}")
    print()

    print("Loading data...", flush=True)
    adata = sc.read_h5ad(in_path)
    assert "counts" in adata.layers, "expected corrected counts in layers['counts']"
    assert "soupx_frac_removed" in adata.obs, \
        "this object did not come out of assemble_soupx.py"
    print(adata, flush=True)
    print()

    # -----------------------------------------------------------------------------------
    # The ambient call, as decided in the notebook. One row per cell of the input, so a
    # table that does not cover it is out of sync rather than partially applicable.
    # -----------------------------------------------------------------------------------
    if args.keep_all:
        keep = np.ones(adata.n_obs, dtype=bool)
    else:
        keep_df = pd.read_csv(keep_path, index_col=0)
        assert KEEP_KEY in keep_df.columns, f"{keep_path} has no '{KEEP_KEY}' column"
        missing = adata.obs_names.difference(keep_df.index)
        assert len(missing) == 0, (
            f"{len(missing):,} cells of {in_path.name} have no row in {keep_path.name}; "
            "re-run the notebook against the same object"
        )
        keep = keep_df.loc[adata.obs_names, KEEP_KEY].to_numpy().astype(bool)

    print(f"ambient call: keeping {keep.sum():,} of {adata.n_obs:,} cells "
          f"({keep.sum() / adata.n_obs:.1%}), dropping {(~keep).sum():,}")
    if (~keep).any():
        dropped_by_type = (adata.obs.loc[~keep, LABEL_KEY].value_counts().head(10))
        print("\nmost-dropped phase-01 labels:")
        print(dropped_by_type.to_string())
    print()

    prior_labels = adata.obs[LABEL_KEY].astype(str).copy()
    prior_fraction = adata.obs["fraction"].astype(str).copy()

    adata = adata[keep].copy()

    # -----------------------------------------------------------------------------------
    # 01_4's annotation, verbatim, on the corrected counts of the kept cells.
    # -----------------------------------------------------------------------------------
    if args.no_reannotate:
        # The ambient call still applies; only the labels are left alone. What this gives up
        # is stated in the module docstring and is not what 05_1 recomputes afterwards: it is
        # the labels that DECIDE 05_1's run - which cells are tested for aneuploidy, and which
        # ones form the diploid baseline. Those are read off the input object and never
        # revisited. Useful as the control arm of a comparison, not as the default.
        print("--no-reannotate: keeping phase 01's cell_type; CellTypist is not run\n")
        new_labels = prior_labels.loc[adata.obs_names]
        raw_labels = (adata.obs["celltypist_predicted"].astype(str)
                      if "celltypist_predicted" in adata.obs
                      else new_labels)
        majority_voting = None
    elif args.reuse_annotation:
        # Resuming, not shortcutting. CellTypist on 619k cells with the majority vote is the
        # expensive half of this step, and it writes its result to `annot_path` BEFORE the
        # h5ad write that is the fragile half - so a run that was interrupted between the two
        # left the labels on disk, complete and consistent with `ambient_keep.csv`. Reading
        # them back reproduces the object that run would have written; recomputing them would
        # only spend the hours again for the same answer.
        print(f"--reuse-annotation: reading labels back from {annot_path}\n", flush=True)
        prev = pd.read_csv(annot_path, index_col=0)
        for col in (LABEL_KEY, "celltypist_predicted", KEEP_KEY):
            assert col in prev.columns, f"{annot_path.name} has no '{col}' column"
        missing = adata.obs_names.difference(prev.index)
        assert len(missing) == 0, (
            f"{len(missing):,} kept cells have no row in {annot_path.name}; that table was "
            "written against a different object - re-run without --reuse-annotation"
        )
        # The table has one row per cell of the INPUT, dropped cells included, and those rows
        # carry no label. A kept cell with no label means the two calls disagree, which is the
        # one way this resume could silently produce a different cell set than the run it
        # resumes.
        assert (prev.loc[adata.obs_names, KEEP_KEY].astype(bool)).all(), (
            f"{annot_path.name} marks some of the cells kept here as dropped; it was written "
            "against a different ambient_keep.csv - re-run without --reuse-annotation"
        )
        new_labels = prev.loc[adata.obs_names, LABEL_KEY].astype(str)
        raw_labels = prev.loc[adata.obs_names, "celltypist_predicted"].astype(str)
        assert not new_labels.isna().any(), "reused cell_type has missing values"
        majority_voting = not args.no_majority_voting
    else:
        print(f"Loading CellTypist model from: {model_path}", flush=True)
        model = models.Model.load(model=str(model_path))

        print("Building temporary CP10K+log1p matrix for CellTypist input...", flush=True)
        adata_ct = adata.copy()
        adata_ct.X = adata_ct.layers["counts"].copy()
        # layers.clear(keep_x=True) rather than `del adata_ct.layers`: the latter warns that a
        # future anndata may drop `.X` with it, and `.X` is the matrix just assigned.
        adata_ct.layers.clear(keep_x=True)
        # The phase-01 PCA, neighbour graph and UMAP describe the uncorrected counts of a larger
        # population. Dropping them is what forces CellTypist to build its over-clustering on
        # this population and this matrix, which is the entire reason for the re-run.
        del adata_ct.obsm, adata_ct.obsp, adata_ct.varm
        sc.pp.normalize_total(adata_ct, target_sum=1e4)
        sc.pp.log1p(adata_ct)

        majority_voting = not args.no_majority_voting
        print(f"Annotating with CellTypist (majority_voting={majority_voting}, CPU)...", flush=True)
        predictions = celltypist.annotate(adata_ct, model=model,
                                          majority_voting=majority_voting, use_GPU=False)
        voted_col = "majority_voting" if majority_voting else "predicted_labels"
        new_labels = predictions.predicted_labels[voted_col].astype(str)
        raw_labels = predictions.predicted_labels["predicted_labels"].astype(str)
        del adata_ct

    # -----------------------------------------------------------------------------------
    # Attach. The phase-01 label is kept under its own name rather than overwritten: it is
    # the only way to say how much the correction moved. The name is deliberately NOT
    # `cell_type_01_4` - see the module docstring: 05_2 writes that column itself, from
    # whatever `cell_type` its input carries.
    # -----------------------------------------------------------------------------------
    adata.obs[PRIOR_LABEL_KEY] = pd.Categorical(prior_labels.loc[adata.obs_names].to_numpy())
    adata.obs["fraction_precorrection"] = pd.Categorical(
        prior_fraction.loc[adata.obs_names].to_numpy(), categories=["imm", "non_imm"])
    adata.obs[LABEL_KEY] = pd.Categorical(new_labels.loc[adata.obs_names].to_numpy())
    adata.obs["celltypist_predicted"] = pd.Categorical(
        raw_labels.loc[adata.obs_names].to_numpy())

    observed = set(adata.obs[LABEL_KEY].astype(str).unique())
    unclassified = observed - IMMUNE_CELL_TYPES - NON_IMMUNE_CELL_TYPES
    assert not unclassified, (
        f"{len(unclassified)} cell_type value(s) in neither lineage set: "
        f"{sorted(unclassified)}; update the sets at the top of this script"
    )
    # `fraction` is a function of `cell_type` (01_4's rule); recomputed, never inherited.
    adata.obs["fraction"] = pd.Categorical(
        np.where(adata.obs[LABEL_KEY].isin(IMMUNE_CELL_TYPES), "imm", "non_imm"),
        categories=["imm", "non_imm"],
    )

    if args.no_reannotate:
        print("labels unchanged (--no-reannotate); cell_type == cell_type_precorrection\n")
    else:
        changed = (adata.obs[LABEL_KEY].astype(str) != adata.obs[PRIOR_LABEL_KEY].astype(str))
        print(f"\nlabel changed for {changed.sum():,} of {adata.n_obs:,} cells "
              f"({changed.mean():.1%})")
        print("\nfraction, phase 01 vs after the correction:")
        print(pd.crosstab(adata.obs["fraction_precorrection"],
                          adata.obs["fraction"]).to_string())
        print("\nthe ten labels that gained or lost the most cells:")
        delta = (adata.obs[LABEL_KEY].value_counts()
                 - adata.obs[PRIOR_LABEL_KEY].value_counts().reindex(
                     adata.obs[LABEL_KEY].cat.categories, fill_value=0))
        print(delta.reindex(delta.abs().sort_values(ascending=False).index).head(10).to_string())
        print()

    # One row per cell of the INPUT object, dropped cells included: the table is the record
    # of what this phase did, and a dropped cell is a decision worth keeping.
    annot = pd.DataFrame({
        PRIOR_LABEL_KEY: prior_labels,
        "fraction_precorrection": prior_fraction,
        KEEP_KEY: keep,
    }, index=prior_labels.index)
    annot[LABEL_KEY] = adata.obs[LABEL_KEY].astype(str).reindex(annot.index)
    annot["celltypist_predicted"] = adata.obs["celltypist_predicted"].astype(str).reindex(annot.index)
    annot["fraction"] = adata.obs["fraction"].astype(str).reindex(annot.index)
    annot.to_csv(annot_path)
    print(f"annotation table -> {annot_path}")

    adata.uns.setdefault("ambient_soupx", {})
    adata.uns["ambient_soupx"].update({
        "cells_before_ambient_filter": int(len(keep)),
        "cells_kept": int(keep.sum()),
        "reannotated": not args.no_reannotate,
        "labels_reused_from": annot_path.name if args.reuse_annotation else "",
        "celltypist_model": model_path.name if not args.no_reannotate else "",
        "majority_voting": bool(majority_voting) if majority_voting is not None else False,
        "cell_type": ("re-annotated on the SoupX-corrected counts of the kept cells; the "
                      "phase-01 label is cell_type_precorrection")
                     if not args.no_reannotate else
                     ("phase 01's label, kept unchanged (--no-reannotate); "
                      "cell_type_precorrection is a copy of it"),
    })

    x_values = adata.X.data
    assert np.allclose(x_values, np.round(x_values)), ".X must hold integer counts"
    assert adata.obs[LABEL_KEY].nunique() > 1, "the annotation collapsed to one label"

    print(f"\nWriting {out_path} ...", flush=True)
    adata.write_h5ad(out_path, compression="gzip")
    print(f"done: {out_path} ({out_path.stat().st_size / 1024**3:.2f} GB)")
    print(f"\n{adata.n_obs:,} cells x {adata.n_vars:,} genes")
    print("\nnext: ../utils/link_shadow_data_dir.sh, then phase 05 with "
          "DATA_DIR=$DATA_DIR/06_amb (see the phase README)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
