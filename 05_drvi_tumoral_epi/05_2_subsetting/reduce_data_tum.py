"""
05_2 feature selection + dimensionality reduction of the malignant subset
(mirrors 04_1/reduce_data_epi.py), plus the DRVI input for 05_3.

The HVGs of 04_1 cannot be reused, and this is the step that makes phase 05 worth running.
They were selected on 74,441 epithelial cells of which roughly half are diploid, so a large
part of that budget is spent on genes that separate malignant from normal epithelium - the
very contrast that is CONSTANT once the subset is malignant-only. Selecting again inside the
tumour is what lets DRVI spend its latent dimensions on states of the tumour rather than on
tumour-versus-normal. Same argument 04_1 makes against 01_5 and 03_1, one compartment down.

Everything else matches 01_5 (same flavor, same 2,000 genes, same PCA), so the objects stay
comparable across phases.

No (re)normalization happens here: .X is consumed as-is from the previous step.

Input : $DATA_DIR/05_tum/<prefix>_norm_cc.h5ad
        .X       = scran log-normalized (from scran_norm_tum.py)
        .layers['counts'] = raw integer counts
        .obs     = 'cohort' (integration batch_key), 'cell_type' (post-CNV label),
                   'cell_type_01_4' (pre-CNV CellTypist label), 'cnv_status', ...
Output: $DATA_DIR/05_tum/<prefix>_reduced.h5ad
        adds .var['highly_variable'] (binary, batch-aware), .obsm['X_pca'],
        neighbors graph, .obsm['X_umap']
        $DATA_DIR/05_tum/<prefix>_hvg_2k_list.csv
        $DATA_DIR/05_tum/<prefix>_hvg_2k.h5ad   <- what DRVI trains on in 05_3

Usage:
This is step 3 of 4 (`reduce`) of subsetting_all.sh.

export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
./subsetting_all.sh            # the whole chain, resuming
./subsetting_all.sh reduce     # only this step, through the wrapper

Standalone (DATA_DIR exported by hand):
python3 reduce_data_tum.py
"""

from __future__ import annotations
import os
import matplotlib
matplotlib.use("Agg")  # headless: figures go to disk via save=, never shown
import numpy as np
import pandas as pd

# scib's hvg_batch uses np.in1d, removed in numpy>=2.4 in favour of np.isin.
if not hasattr(np, "in1d"):
    np.in1d = np.isin

import anndata2ri
from rpy2.robjects import conversion, default_converter

def _activate():
    conversion.set_conversion(default_converter + anndata2ri.converter)

def _deactivate():
    conversion.set_conversion(default_converter)

anndata2ri.activate = _activate
anndata2ri.deactivate = _deactivate

import scanpy as sc
import scib

import cell_set as C

C.banner("05_2 feature selection + reduction")
IN_PATH = C.path("_norm_cc.h5ad")
OUT_PATH = C.path("_reduced.h5ad")
HVG_CSV_PATH = C.path("_hvg_2k_list.csv")
HVG_H5AD_PATH = C.path("_hvg_2k.h5ad")

try:
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
except NameError:
    SCRIPT_DIR = os.getcwd()
PHASE_DIR = os.path.dirname(SCRIPT_DIR)  # 05_drvi_tumoral_epi/
FIG_DIR = os.environ.get("FIG_DIR", os.path.join(PHASE_DIR, "figures", "05_2_reduce_data"))
os.makedirs(FIG_DIR, exist_ok=True)
sc.settings.figdir = FIG_DIR
sc.settings.set_figure_params(dpi=300, facecolor="white")

adata = sc.read_h5ad(IN_PATH)
print(f"Loaded {adata.n_obs:,} cells x {adata.n_vars:,} genes", flush=True)

assert C.BATCH_KEY in adata.obs, f"missing batch key {C.BATCH_KEY!r}"
assert C.LABEL_KEY in adata.obs, f"missing label key {C.LABEL_KEY!r}"
assert "counts" in adata.layers, "expected raw counts in .layers['counts']"

cells_per_batch = adata.obs[C.BATCH_KEY].value_counts()
cells_per_batch = cells_per_batch[cells_per_batch > 0]
print(f"{C.BATCH_KEY}: {cells_per_batch.size} batches, "
      f"min {cells_per_batch.min():,} cells ({cells_per_batch.idxmin()})", flush=True)
assert cells_per_batch.min() >= C.MIN_CELLS_PER_COHORT, (
    f"a cohort under {C.MIN_CELLS_PER_COHORT} cells reached this step "
    f"({cells_per_batch.min()}); subset_and_qc.ipynb should have dropped it"
)

# Batch-aware HVG + PCA + neighbours + UMAP, parameters identical to 01_5 and 04_1.
scib.preprocessing.reduce_data(
    adata,
    batch_key=C.BATCH_KEY,
    flavor="cell_ranger",
    n_top_genes=C.N_HVGS,
    n_bins=20,
    pca=True,
    pca_comps=50,
    svd_solver="arpack",
    overwrite_hvg=True,
    neighbors=True,
    use_rep="X_pca",
    umap=True,
)
n_hvg = int(adata.var["highly_variable"].sum())
print(f"HVGs selected: {n_hvg}", flush=True)
assert n_hvg == C.N_HVGS, f"expected {C.N_HVGS} HVGs, got {n_hvg}"

hvg_genes = adata.var_names[adata.var["highly_variable"]]
pd.Series(hvg_genes).to_csv(HVG_CSV_PATH, index=False, header=False)
print(f"Wrote {len(hvg_genes)} HVGs to {HVG_CSV_PATH}", flush=True)

# Overlap with the earlier selections. The one that matters for this phase is 04_1: it is the
# direct measure of how much of the epithelial HVG budget was being spent on the
# malignant-versus-normal contrast rather than on structure inside the tumour. A LOW overlap
# is the expected and desired outcome - it is the quantitative version of the argument in the
# docstring. Purely informative, skipped when a list is not on disk.
for label, rel in (("01_5 full-object", "shiao_hvg_2k_unintegrated_list.csv"),
                   ("03_1 non-immune ", os.path.join("03_nonimm", "shiao_nonimm_hvg_2k_list.csv")),
                   ("04_1 epithelial ", os.path.join("04_epi", "shiao_epi_hvg_2k_list.csv"))):
    p = os.path.join(C.data_dir(), rel)
    if os.path.exists(p):
        other = set(pd.read_csv(p, header=None)[0].astype(str))
        shared = len(set(hvg_genes) & other)
        print(f"Overlap with the {label} HVGs: {shared}/{C.N_HVGS} "
              f"({100 * shared / C.N_HVGS:.1f}%)", flush=True)

# Diagnostic figures. `cell_type` is constant on the malignant subset, so the label panels
# are drawn on the pre-CNV label instead - it is the only non-constant annotation here, and
# what it shows is which normal state each tumour cell was being mistaken for.
suffix = f"_{C.compartment()}"
label_for_plots = C.LABEL_KEY if adata.obs[C.LABEL_KEY].nunique() > 1 else C.PRIOR_LABEL_KEY
print(f"colouring the label panels by {label_for_plots!r} "
      f"({adata.obs[label_for_plots].nunique()} levels)", flush=True)

sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, save=f"{suffix}_elbow.png")
sc.pl.pca(adata, color=label_for_plots, save=f"{suffix}_{label_for_plots}.png")
sc.pl.umap(adata, color=label_for_plots, save=f"{suffix}_{label_for_plots}.png")
sc.pl.umap(adata, color=C.BATCH_KEY, save=f"{suffix}_{C.BATCH_KEY}.png")
if "phase" in adata.obs:
    sc.pl.umap(adata, color="phase", save=f"{suffix}_phase.png")
if C.STATUS_KEY in adata.obs and adata.obs[C.STATUS_KEY].nunique() > 1:
    sc.pl.umap(adata, color=C.STATUS_KEY, save=f"{suffix}_cnv_status.png")
for col in ("cnv_score", "cnv_corr"):
    if col in adata.obs:
        sc.pl.umap(adata, color=col, save=f"{suffix}_{col}.png")

adata.write_h5ad(OUT_PATH, compression="gzip")
print(f"Wrote {OUT_PATH}", flush=True)

# The DRVI input: the same object restricted to the 2,000 HVGs. Done here rather than in a
# separate script because this is the step that made the selection, so the .h5ad and the flag
# in .var cannot drift apart.
print("\nBuilding the DRVI input (2,000 HVGs)...", flush=True)
mask = adata.var["highly_variable"].to_numpy()
hvg = adata[:, mask].copy()

assert hvg.n_vars == C.N_HVGS, f"got {hvg.n_vars} genes after subsetting"
assert hvg.n_obs == adata.n_obs, "cells were lost during the gene subset"

# DRVI reads .layers['counts'] and ignores .X entirely; if anything upstream replaced that
# layer the model would train for hours on the wrong input. Fail here instead.
counts = hvg.layers["counts"]
counts_values = counts.data if hasattr(counts, "data") else np.asarray(counts).ravel()
assert np.allclose(counts_values, np.floor(counts_values)), (
    "layers['counts'] contains non-integer values: it is no longer raw counts. "
    "DRVI reads this layer directly, so stop here rather than train on normalized data."
)
assert counts_values.min() >= 0, "layers['counts'] contains negative values"

x_values = hvg.X.data if hasattr(hvg.X, "data") else np.asarray(hvg.X).ravel()
assert not np.allclose(x_values, np.floor(x_values)), (
    ".X looks like integer counts, but it should be scran log-normalized expression. "
    "Check whether .X and layers['counts'] were swapped upstream."
)

hvg.write_h5ad(HVG_H5AD_PATH, compression="gzip")
print(f"Wrote {HVG_H5AD_PATH} "
      f"({hvg.n_obs:,} x {hvg.n_vars:,}, "
      f"{os.path.getsize(HVG_H5AD_PATH) / 1024 ** 3:.2f} GB on disk)", flush=True)
