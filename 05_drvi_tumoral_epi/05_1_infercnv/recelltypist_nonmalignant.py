#!/usr/bin/env python3
"""05_1 step 4: re-run CellTypist on the non-malignant cells only.

Why this exists. `cell_type` was annotated in 01_4 with Cells_Adult_Breast.pkl, a CellTypist
model trained on the Kumar et al. 2023 adult breast atlas. That atlas is NORMAL breast, so
the model has no malignant class: a tumour cell is not left unassigned, it is assigned to
whichever normal state it resembles most, and it then votes in the majority-voting step of
every neighbourhood it sits in. The labels of the epithelial compartment are therefore
contaminated twice over - once directly, and once through the voting.

Once inferCNV has separated the aneuploid cells (call_malignant.ipynb), both problems are
fixable in one move: take the malignant cells out of the vocabulary entirely, give them the
single label `malignant`, and re-run exactly the 01_4 annotation on what is left. The
non-malignant cells are then annotated by a normal-breast model against a normal-breast
population, which is the regime the model was trained for.

Two things about this that are worth being explicit on:

  * majority voting is why the re-run is not a no-op. CellTypist's over-clustering is
    recomputed on the reduced population, so a normal luminal cell that previously sat in a
    neighbourhood dominated by tumour cells is now voted on by its actual neighbours. Cells
    can and do change label; that is the point, not a bug.
  * `not_tested` cells go through the re-annotation. Those are the cells no inferCNV run
    covered: B/plasma cells (deliberately excluded from every reference, see
    prepare_infercnv_input.py) and every cell of a cohort with too little epithelium to run.
    For an immune cell that is harmless. For an EPITHELIAL cell it is not - it has simply
    not been tested for aneuploidy - so `cnv_status` travels in the output next to
    `cell_type_cnv`, and 05_2 decides what to do with those cells rather than this script
    silently promoting them to normal epithelium.

Nothing is written back into shiao.h5ad. The output is one per-cell table, which is what a
CNV call is: metadata. 05_2 joins it onto the object it builds.

Input : $DATA_DIR/shiao.h5ad                      (`layers['counts']` = raw counts)
        $DATA_DIR/05_tum/cnv_status.csv       written by call_malignant.ipynb
        $DATA_DIR/Cells_Adult_Breast.pkl          the same model as 01_4
Output: $DATA_DIR/05_tum/cell_annotation_cnv.csv
          cell, cohort, cell_type (01_4, unchanged), cnv_status,
          cell_type_cnv (the new label; 'malignant' for the malignant cells),
          celltypist_predicted_cnv (raw per-cell prediction, before voting, for QC)

## Why the re-run covers EVERY non-malignant cell, not just the epithelial ones

It is the expensive choice - ~582,000 cells and a couple of hours, most of it spent
recomputing immune labels that will barely move once 37,014 cells out of 619,693 are gone -
and it is deliberate, for two reasons.

The first is that "just the epithelial ones" is not a set this step can name. The only
definition of epithelial available is `cell_type` from 01_4, which is the contaminated
annotation this step exists to repair; gating the input on it would lock in exactly the
errors being repaired, and a cell that today carries a non-epithelial label but should carry
an epithelial one would never get the chance to change.

The second is that restricting the input would change the neighbourhood structure the vote
is computed over on top of everything else, so a cell that moved could have moved because
the tumour cells stopped voting or because the immune cells did.

## What this run does NOT let you claim  (measured, not assumed)

An earlier version of this docstring said that running the same method on the same
population minus the malignant cells makes every label change attributable to that removal
and to nothing else. **That is wrong, and the output of the first run disproves it.**

29.5% of the 582,679 non-malignant cells came back with a different label, but:

  * 95.3% of those changes stay INSIDE the same lineage (immune -> immune, epithelial ->
    epithelial); only 4.7% cross one;
  * 26,327 of them - 15.3% of all changes - are CD4 <-> CD8 swaps, which no amount of
    removing epithelial tumour cells can cause. A T cell does not change co-receptor.

The cause is the majority vote itself. CellTypist has no neighbourhood graph to reuse here
(this script drops the phase-01 one on purpose), so it builds one and runs leiden at
resolution 30 over it. That partition is stochastic and its boundaries move between runs;
each cluster then takes its majority label, and a cluster whose boundary shifted can flip
wholesale. The removal of 37,014 cells perturbs the PCA and the graph enough to reshuffle
boundaries everywhere, including in compartments the tumour never touched.

So the per-cell labels here are usable, but the CHURN is not a measurement of contamination.
Quoting the 29.5% as "what the tumour cells were doing to the annotation" would be wrong.

What would make it a measurement: compute the over-clustering ONCE on the full object, then
pass the same partition to both runs via `celltypist.annotate(over_clustering=...)`, so the
only thing that differs is which cells vote inside each fixed cluster. That is a change to
this script, not a re-run of it, and it is not done yet.

`--scope non-immune` exists for a fast re-run (the `fraction == 'non_imm'` compartment,
~140k cells) and is NOT the default: it buys speed by giving up the comparison above.

Local usage (benchmark-py-r):
    export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
    python3 recelltypist_nonmalignant.py
    python3 recelltypist_nonmalignant.py --no-majority-voting   # per-cell only, much faster
    python3 recelltypist_nonmalignant.py --scope non-immune     # faster, but see above
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

import celltypist
from celltypist import models

sc.settings.verbosity = 1

MALIGNANT_LABEL = "malignant"
STATUS_KEY = "cnv_status"
LABEL_KEY = "cell_type"
BATCH_KEY = "cohort"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--no-majority-voting", action="store_true",
                   help="skip CellTypist's over-clustering vote (fast, but not what 01_4 did)")
    p.add_argument("--scope", choices=("all", "non-immune"), default="all",
                   help="which non-malignant cells to re-annotate; 'all' (default) reproduces "
                        "01_4's population minus the malignant cells, see the module docstring")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    data_dir = Path(os.environ["DATA_DIR"]).expanduser().resolve()
    in_path = data_dir / "shiao.h5ad"
    cnv_dir = data_dir / "05_tum"
    status_path = cnv_dir / "cnv_status.csv"
    model_path = data_dir / "Cells_Adult_Breast.pkl"
    out_path = cnv_dir / "cell_annotation_cnv.csv"

    for path, hint in ((in_path, "phase 01"),
                       (status_path, "call_malignant.ipynb"),
                       (model_path, "01_4")):
        assert path.exists(), f"missing {path} (produced by {hint})"

    print(f"DATA_DIR: {data_dir}")
    print(f"Input   : {in_path}")
    print(f"Status  : {status_path}")
    print(f"Output  : {out_path}")
    print()

    status = pd.read_csv(status_path, index_col=0)
    assert STATUS_KEY in status.columns, f"{status_path} has no '{STATUS_KEY}' column"

    print("Loading data...", flush=True)
    adata = sc.read_h5ad(in_path)
    assert "counts" in adata.layers, "expected raw counts in layers['counts']"

    # Every cell must have a status. call_malignant.ipynb writes one row per cell of
    # shiao.h5ad, 'not_tested' included, so a missing row means the two are out of sync.
    missing = adata.obs_names.difference(status.index)
    assert len(missing) == 0, (
        f"{len(missing):,} cells of shiao.h5ad have no row in {status_path.name}; "
        "re-run call_malignant.ipynb against the same object"
    )
    adata.obs[STATUS_KEY] = status.loc[adata.obs_names, STATUS_KEY].to_numpy()
    print(adata.obs[STATUS_KEY].value_counts().to_string())
    print()

    is_malignant = (adata.obs[STATUS_KEY] == "malignant").to_numpy()

    # Which non-malignant cells get re-annotated. Anything left out keeps its 01_4 label, and
    # that is recorded in the output rather than silently implied.
    if args.scope == "all":
        reannotate = ~is_malignant
    else:
        assert "fraction" in adata.obs, "obs['fraction'] is needed by --scope non-immune"
        reannotate = (~is_malignant) & (adata.obs["fraction"] == "non_imm").to_numpy()

    print(f"malignant: {is_malignant.sum():,} cells -> single label '{MALIGNANT_LABEL}'")
    print(f"re-annotated (--scope {args.scope}): {reannotate.sum():,} cells")
    kept_as_is = (~is_malignant) & (~reannotate)
    if kept_as_is.any():
        print(f"keeping the 01_4 label unchanged: {kept_as_is.sum():,} cells")
    print()

    # ---------------------------------------------------------------------------------
    # The 01_4 annotation, verbatim, on the non-malignant subset. CellTypist expects
    # CP10K + log1p (its training normalization), NOT the scran .X carried by this object,
    # so the temporary matrix is rebuilt from raw counts exactly as celltypist_annotation.py
    # does, used, and thrown away.
    # ---------------------------------------------------------------------------------
    print(f"Loading CellTypist model from: {model_path}", flush=True)
    model = models.Model.load(model=str(model_path))

    print("Building temporary CP10K+log1p matrix for CellTypist input...", flush=True)
    adata_ct = adata[reannotate].copy()
    adata_ct.X = adata_ct.layers["counts"].copy()
    # The phase-01 PCA, neighbour graph and UMAP were computed on the population that still
    # contained the tumour cells. Dropping them is what forces CellTypist to build its
    # over-clustering on the reduced population, which is the entire reason for this re-run.
    # layers.clear(keep_x=True) rather than `del adata_ct.layers`: the latter warns that a
    # future anndata may drop `.X` with it, and `.X` is the matrix that was just assigned.
    adata_ct.layers.clear(keep_x=True)
    del adata_ct.obsm, adata_ct.obsp, adata_ct.varm
    sc.pp.normalize_total(adata_ct, target_sum=1e4)
    sc.pp.log1p(adata_ct)

    majority_voting = not args.no_majority_voting
    print(f"Annotating with CellTypist (majority_voting={majority_voting}, CPU)...", flush=True)
    predictions = celltypist.annotate(
        adata_ct,
        model=model,
        majority_voting=majority_voting,
        use_GPU=False,
    )
    voted_col = "majority_voting" if majority_voting else "predicted_labels"
    new_labels = predictions.predicted_labels[voted_col].astype(str)
    raw_labels = predictions.predicted_labels["predicted_labels"].astype(str)
    del adata_ct

    # ---------------------------------------------------------------------------------
    # Assemble the per-cell table.
    # ---------------------------------------------------------------------------------
    out = pd.DataFrame(index=adata.obs_names)
    out[BATCH_KEY] = adata.obs[BATCH_KEY].astype(str).to_numpy()
    out[LABEL_KEY] = adata.obs[LABEL_KEY].astype(str).to_numpy()
    out[STATUS_KEY] = adata.obs[STATUS_KEY].to_numpy()

    # Start from the 01_4 label, overwrite the malignant cells with the single label, then
    # overwrite the re-annotated ones. Under --scope all the middle group is empty; under
    # --scope non-immune it is the immune cells, which keep what 01_4 gave them.
    out["cell_type_cnv"] = out[LABEL_KEY]
    out["celltypist_predicted_cnv"] = out[LABEL_KEY]
    out.loc[is_malignant, "cell_type_cnv"] = MALIGNANT_LABEL
    out.loc[is_malignant, "celltypist_predicted_cnv"] = MALIGNANT_LABEL
    out.loc[new_labels.index, "cell_type_cnv"] = new_labels.to_numpy()
    out.loc[raw_labels.index, "celltypist_predicted_cnv"] = raw_labels.to_numpy()
    out["reannotated"] = reannotate

    assert (out.loc[is_malignant, "cell_type_cnv"] == MALIGNANT_LABEL).all()
    assert (out.loc[~is_malignant, "cell_type_cnv"] != MALIGNANT_LABEL).all(), \
        "the model produced a label literally called 'malignant'; rename MALIGNANT_LABEL"
    assert out["cell_type_cnv"].notna().all()

    print()
    print("cell_type_cnv distribution:")
    print(out["cell_type_cnv"].value_counts().to_string())
    print()
    print("How the re-annotation moved the non-malignant cells (01_4 -> now), top 20 changes:")
    nm = out.loc[reannotate]
    changed = nm[nm[LABEL_KEY] != nm["cell_type_cnv"]]
    print(f"{len(changed):,} of {len(nm):,} re-annotated cells changed label "
          f"({len(changed) / max(len(nm), 1):.1%})")
    print(changed.groupby([LABEL_KEY, "cell_type_cnv"]).size()
          .sort_values(ascending=False).head(20).to_string())

    out.index.name = "cell"
    out.to_csv(out_path)
    print()
    print(f"Wrote {out_path} ({len(out):,} cells)")
    print("next: 05_2_subsetting reads this table instead of `cell_type` from shiao.h5ad")
    return 0


if __name__ == "__main__":
    sys.exit(main())
