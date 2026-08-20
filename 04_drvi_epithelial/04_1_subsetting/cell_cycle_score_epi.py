"""
04_1 cell cycle scoring of the epithelial subset (mirrors 01_4/cell_cycle_score.py).

Re-scored rather than inherited from phase 01: score_genes_cell_cycle compares the
S / G2M gene sets against a reference set of genes sampled from expression bins
computed on the *current* object. Those bins are built from the mean expression of
the ~74k epithelial cells here, of 620k mixed cells there, so the same gene sets
give different scores. Inheriting `phase` would mean carrying a score calibrated
on a population that no longer exists.

Nothing else from 01_4 is repeated: CellTypist is not re-run (`cell_type` is
inherited and sub-annotation belongs to 04_2), and `fraction_reassignment` is
meaningless here because `fraction` is constant ('non_imm') by construction.

Input : $DATA_DIR/04_epi/shiao_epi_norm.h5ad
        Output of scran_norm_epi.py. .X = scran log-normalized (scoring
        requires log-norm, not raw); gene symbols as var_names.
Output: $DATA_DIR/04_epi/shiao_epi_norm_cc.h5ad
        Adds .obs['S_score'], .obs['G2M_score'], .obs['phase'].

Requires the Tirosh/Regev cell cycle gene list (97 genes: 43 S + 54 G2M),
one gene symbol per line, path via $CC_GENES env var.

Usage:
This is step 2 of 4 (`cc`) of subsetting_all.sh, which is the intended way to
run it: the wrapper chains the four scripts in order, skips the steps whose
output already exists, and tees everything to 04_1_subsetting/logs/. It also
defaults $CC_GENES to the file below, so only DATA_DIR has to be exported.

export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
./subsetting_all.sh            # the whole chain, resuming
./subsetting_all.sh cc         # only this step, through the wrapper

Standalone, outside the chain (no logging, no resume, both variables must be
exported by hand):
export CC_GENES=$DATA_DIR/regev_lab_cell_cycle_genes.txt
python3 cell_cycle_score_epi.py
"""

from __future__ import annotations
import os
import numpy as np
import scanpy as sc
import anndata as ad

sc.settings.verbosity = 1

DATA_DIR = os.environ["DATA_DIR"]
EPI_DIR = os.path.join(DATA_DIR, "04_epi")
IN_PATH = os.path.join(EPI_DIR, "shiao_epi_norm.h5ad")
OUT_PATH = os.path.join(EPI_DIR, "shiao_epi_norm_cc.h5ad")

# Cell cycle gene list from env. Tirosh/Regev standard file:
# 97 lines, first 43 = S phase, remaining 54 = G2/M phase.
CC_GENES = os.environ["CC_GENES"]

print("Loading data...", flush=True)
adata = sc.read_h5ad(IN_PATH)
print(adata, flush=True)

# Scoring must run on log-normalized expression (scran output in .X), not raw counts.
# Heuristic: raw counts are integers; log-norm values are not.
x_sample = adata.X[:100].data if hasattr(adata.X[:100], "data") else np.asarray(adata.X[:100]).ravel()
x_sample = x_sample[np.isfinite(x_sample)]
assert not np.allclose(x_sample, np.round(x_sample)), \
    ".X looks like raw integer counts; cell cycle scoring expects log-normalized .X"

# Load and split the gene list.
with open(CC_GENES) as f:
    cell_cycle_genes = [x.strip() for x in f if x.strip()]
assert len(cell_cycle_genes) == 97, \
    f"expected 97 cell cycle genes (43 S + 54 G2M), got {len(cell_cycle_genes)}"
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]

# Keep only genes present in the dataset (var_names are gene symbols).
# The overlap can be lower than in 01_4: the subset went through a second
# filter_genes pass, so cell cycle genes expressed only outside the epithelium are gone.
# Report the survivors; a much smaller set would weaken the score and must be documented.
s_present = [g for g in s_genes if g in adata.var_names]
g2m_present = [g for g in g2m_genes if g in adata.var_names]
print(f"S genes:   {len(s_present)}/{len(s_genes)} present", flush=True)
print(f"G2M genes: {len(g2m_present)}/{len(g2m_genes)} present", flush=True)
missing = sorted(set(s_genes + g2m_genes) - set(s_present) - set(g2m_present))
if missing:
    print(f"Missing cell cycle genes (not in var_names): {missing}", flush=True)

print("Scoring cell cycle...", flush=True)
# score_genes_cell_cycle: per-cell mean expression of S and G2M gene sets relative to
# a randomly sampled reference set of comparably-expressed genes (Tirosh et al. 2016).
# Phase assigned by dominant score; G1 when both scores are non-positive.
# random_state=0 makes the reference sampling reproducible (documented in Methods).
sc.tl.score_genes_cell_cycle(
    adata,
    s_genes=s_present,
    g2m_genes=g2m_present,
    random_state=0,
)

print("Cell cycle phase distribution:", flush=True)
print(adata.obs["phase"].value_counts(), flush=True)

# Post-condition checks
for col in ("S_score", "G2M_score", "phase"):
    assert col in adata.obs, f"scoring did not add {col}"

print("Saving...", flush=True)
adata.write_h5ad(OUT_PATH, compression="gzip")
print(f"Wrote {OUT_PATH}", flush=True)
