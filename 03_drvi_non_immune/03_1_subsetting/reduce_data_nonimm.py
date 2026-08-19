"""
03_1 feature selection + dimensionality reduction of the non-immune subset
(mirrors 01_5/scib_reduce_data.py), plus the DRVI input for 03_2.

The HVGs of phase 01 cannot be reused: they were selected on an object where
443,083 of 619,693 cells were immune, so the 2,000 genes in
`shiao_hvg_2k_unintegrated_list.csv` are dominated by immune variation. Selecting
again on the 176,610 non-immune cells is the whole point of this subsetting phase
- it is what lets DRVI spend its latent dimensions on epithelial / stromal /
vascular biology instead of on lymphocyte programs.

Everything else matches 01_5 (same flavor, same 2,000 genes, same PCA), so the
non-immune object stays comparable with the full one.

No (re)normalization happens here: .X is consumed as-is from the previous step.

Input : $DATA_DIR/03_nonimm/shiao_nonimm_norm_cc.h5ad
        .X       = scran log-normalized (from scran_norm_nonimm.py)
        .layers['counts'] = raw integer counts
        .var     = gene symbols as var_names, Ensembl in .var['gene_ids']
        .obs     = 'cohort' (34 patients, integration batch_key), 'cell_type'
                   (CellTypist majority_voting label_key), etc.
Output: $DATA_DIR/03_nonimm/shiao_nonimm_reduced.h5ad
        adds .var['highly_variable'] (binary, batch-aware), .obsm['X_pca'],
        neighbors graph, .obsm['X_umap']
        $DATA_DIR/03_nonimm/shiao_nonimm_hvg_2k_list.csv
        gene symbols of the selected HVGs, one per line, no header
        $DATA_DIR/03_nonimm/shiao_nonimm_hvg_2k.h5ad
        the same object restricted to those 2,000 genes, raw counts preserved in
        .layers['counts'] -- this is what DRVI trains on in 03_2

Usage:
This is step 3 of 4 (`reduce`) of subsetting_all.sh, which is the intended way
to run it: the wrapper chains the four scripts in order, skips the steps whose
output already exists, and tees everything to 03_1_subsetting/logs/.

export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
./subsetting_all.sh            # the whole chain, resuming
./subsetting_all.sh reduce     # only this step, through the wrapper

Standalone, outside the chain (no logging, no resume, DATA_DIR must be
exported by hand):
python3 reduce_data_nonimm.py
"""

from __future__ import annotations
import os
import matplotlib
matplotlib.use("Agg")  # headless backend: figures are only written to disk via save=, never shown/blocking
import numpy as np
import pandas as pd

# scib's hvg_batch (called by reduce_data with overwrite_hvg=True) uses
# np.in1d, removed in numpy>=2.4 in favour of np.isin. Restore the alias so
# scib keeps working on modern numpy.
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

N_HVGS = 2000
BATCH_KEY = "cohort"
LABEL_KEY = "cell_type"

# Resolve DATA_DIR from environment, never hardcoded
DATA_DIR = os.environ["DATA_DIR"]
NONIMM_DIR = os.path.join(DATA_DIR, "03_nonimm")
IN_PATH = os.path.join(NONIMM_DIR, "shiao_nonimm_norm_cc.h5ad")
OUT_PATH = os.path.join(NONIMM_DIR, "shiao_nonimm_reduced.h5ad")
HVG_CSV_PATH = os.path.join(NONIMM_DIR, "shiao_nonimm_hvg_2k_list.csv")
HVG_H5AD_PATH = os.path.join(NONIMM_DIR, "shiao_nonimm_hvg_2k.h5ad")

# Figure dir anchored to the script location, fallback to cwd for interactive use
try:
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
except NameError:
    SCRIPT_DIR = os.getcwd()
PHASE_DIR = os.path.dirname(SCRIPT_DIR)  # 03_drvi_non_immune/
FIG_DIR = os.environ.get("FIG_DIR", os.path.join(PHASE_DIR, "figures", "03_1_reduce_data"))
os.makedirs(FIG_DIR, exist_ok=True)
sc.settings.figdir = FIG_DIR
sc.settings.set_figure_params(dpi=300, facecolor="white")

adata = sc.read_h5ad(IN_PATH)
print(f"Loaded {adata.n_obs:,} cells x {adata.n_vars:,} genes", flush=True)

assert BATCH_KEY in adata.obs, f"missing batch key {BATCH_KEY!r}"
assert LABEL_KEY in adata.obs, f"missing label key {LABEL_KEY!r}"
assert "counts" in adata.layers, "expected raw counts in .layers['counts']"

# scib's hvg_batch selects per batch and then merges, so a batch with very few
# cells contributes a noisy ranking. Report the smallest ones rather than fail:
# the decision of whether to drop them was already taken in subset_and_qc.ipynb.
cells_per_batch = adata.obs[BATCH_KEY].value_counts()
print(f"{BATCH_KEY}: {cells_per_batch.size} batches, "
      f"min {cells_per_batch.min():,} cells ({cells_per_batch.idxmin()})", flush=True)

# Batch-aware feature selection + dimensionality reduction.
# batch_key='cohort' -> HVG selection done per patient then merged (scib hvg_batch);
# 'cohort' is the technical integration batch, not biological.
# flavor='cell_ranger', n_top_genes=2000, n_bins=20 -> scib/Luecken default HVG config,
# identical to 01_5 so the gene sets are comparable.
# pca_comps=50, svd_solver='arpack' -> deterministic PCA.
# neighbors + umap computed on X_pca
scib.preprocessing.reduce_data(
    adata,
    batch_key=BATCH_KEY,
    flavor="cell_ranger",
    n_top_genes=N_HVGS,
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
assert n_hvg == N_HVGS, f"expected {N_HVGS} HVGs, got {n_hvg}"

# Save the selected HVGs for reuse outside this AnnData: plain list of gene symbols,
# one per line, no header (same format as 01_5, so the 02 helpers can read it too).
hvg_genes = adata.var_names[adata.var["highly_variable"]]
pd.Series(hvg_genes).to_csv(HVG_CSV_PATH, index=False, header=False)
print(f"Wrote {len(hvg_genes)} HVGs to {HVG_CSV_PATH}", flush=True)

# How much of the selection is shared with the full-object HVGs of 01_5. Purely
# informative (a low overlap is the expected, desired outcome here), and skipped
# if the phase 01 list isn't around.
full_hvg_path = os.path.join(DATA_DIR, "shiao_hvg_2k_unintegrated_list.csv")
if os.path.exists(full_hvg_path):
    full_hvg = set(pd.read_csv(full_hvg_path, header=None)[0].astype(str))
    shared = len(set(hvg_genes) & full_hvg)
    print(f"Overlap with the 01_5 full-object HVGs: {shared}/{N_HVGS} "
          f"({100 * shared / N_HVGS:.1f}%)", flush=True)

# Diagnostic figures
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, save="_nonimm_elbow_unintegrated.png")
sc.pl.pca(adata, color=LABEL_KEY, save="_nonimm_cell_type_unintegrated.png")
sc.pl.umap(adata, color=LABEL_KEY, save="_nonimm_cell_type_unintegrated.png")
sc.pl.umap(adata, color=BATCH_KEY, save="_nonimm_cohort_unintegrated.png")
if "phase" in adata.obs:
    sc.pl.umap(adata, color="phase", save="_nonimm_phase_unintegrated.png")

adata.write_h5ad(OUT_PATH, compression="gzip")
print(f"Wrote {OUT_PATH}", flush=True)


# The DRVI input: the same object restricted to the 2,000 HVGs.
#
# Done here rather than in a separate script because this is the step that made
# the selection, so the .h5ad and the flag in .var cannot drift apart. The gene
# subset is taken with a boolean mask on the flag, not by indexing with the CSV
# order, so the genes keep the order of the source object and var['highly_variable'],
# varm['PCs'] and the CSV all stay consistent.
#
# obsm['X_pca'] and uns['pca'] are carried over deliberately: the PCA above was
# computed on exactly these 2,000 genes, so it stays valid after the subset.
print("\nBuilding the DRVI input (2,000 HVGs)...", flush=True)
mask = adata.var["highly_variable"].to_numpy()
hvg = adata[:, mask].copy()

assert hvg.n_vars == N_HVGS, f"got {hvg.n_vars} genes after subsetting"
assert hvg.n_obs == adata.n_obs, "cells were lost during the gene subset"

# The check this block exists for: DRVI reads .layers['counts'] and ignores .X
# entirely. If anything upstream replaced that layer with normalized values, the
# model would train for hours on the wrong input. Fail here instead.
counts = hvg.layers["counts"]
counts_values = counts.data if hasattr(counts, "data") else np.asarray(counts).ravel()
assert np.allclose(counts_values, np.floor(counts_values)), (
    "layers['counts'] contains non-integer values: it is no longer raw counts. "
    "DRVI reads this layer directly, so stop here rather than train on normalized data."
)
assert counts_values.min() >= 0, "layers['counts'] contains negative values"

# The mirror image: .X must NOT be raw counts (it is scran log-normalized).
x_values = hvg.X.data if hasattr(hvg.X, "data") else np.asarray(hvg.X).ravel()
assert not np.allclose(x_values, np.floor(x_values)), (
    ".X looks like integer counts, but it should be scran log-normalized expression. "
    "Check whether .X and layers['counts'] were swapped upstream."
)

hvg.write_h5ad(HVG_H5AD_PATH, compression="gzip")
print(f"Wrote {HVG_H5AD_PATH} "
      f"({hvg.n_obs:,} x {hvg.n_vars:,}, "
      f"{os.path.getsize(HVG_H5AD_PATH) / 1024 ** 3:.2f} GB on disk)", flush=True)
