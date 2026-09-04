"""
05_2 cell cycle scoring of the malignant subset (mirrors 04_1/cell_cycle_score_epi.py).

Re-scored rather than inherited, for the reason that applies at every subsetting step:
`score_genes_cell_cycle` compares the S / G2M gene sets against reference genes sampled
from expression bins computed on the CURRENT object. Those bins come from ~36k malignant
cells here and from 74,441 mixed epithelial cells in 04_1, so the same gene lists give
different scores. Inheriting `phase` would mean carrying a score calibrated on a population
that no longer exists.

This matters more here than it did in 04. The malignant subset is enriched for proliferating
cells - 18,943 of its cells carried the `Lumsec-prol` label from 01_4 - so the bins, and
therefore the baseline the score is measured against, shift substantially. The cell cycle is
also a NAMED confounder of the `scie` collection, which 05_9 tests explicitly; that test is
only meaningful if the score it uses was calibrated on the cells it is testing.

Input : $DATA_DIR/05_tum/<prefix>_norm.h5ad
        .X = scran log-normalized (scoring requires log-norm, not raw).
Output: $DATA_DIR/05_tum/<prefix>_norm_cc.h5ad
        Adds .obs['S_score'], .obs['G2M_score'], .obs['phase'].

Requires the Tirosh/Regev cell cycle gene list (97 genes: 43 S + 54 G2M), one symbol per
line, path via $CC_GENES.

Usage:
This is step 2 of 4 (`cc`) of subsetting_all.sh, which defaults $CC_GENES for you.

export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
./subsetting_all.sh            # the whole chain, resuming
./subsetting_all.sh cc         # only this step, through the wrapper

Standalone (both variables exported by hand):
export CC_GENES=$DATA_DIR/regev_lab_cell_cycle_genes.txt
python3 cell_cycle_score_tum.py
"""

from __future__ import annotations
import os
import numpy as np
import scanpy as sc

import cell_set as C

sc.settings.verbosity = 1

C.banner("05_2 cell cycle scoring")
IN_PATH = C.path("_norm.h5ad")
OUT_PATH = C.path("_norm_cc.h5ad")
CC_GENES = os.environ["CC_GENES"]

print("Loading data...", flush=True)
adata = sc.read_h5ad(IN_PATH)
print(adata, flush=True)

# Scoring must run on log-normalized expression, not raw counts. Heuristic: raw counts are
# integers, log-norm values are not.
x_sample = adata.X[:100].data if hasattr(adata.X[:100], "data") else np.asarray(adata.X[:100]).ravel()
x_sample = x_sample[np.isfinite(x_sample)]
assert not np.allclose(x_sample, np.round(x_sample)), \
    ".X looks like raw integer counts; cell cycle scoring expects log-normalized .X"

with open(CC_GENES) as f:
    cell_cycle_genes = [x.strip() for x in f if x.strip()]
assert len(cell_cycle_genes) == 97, \
    f"expected 97 cell cycle genes (43 S + 54 G2M), got {len(cell_cycle_genes)}"
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]

# The overlap can be lower than in 04_1: this subset went through another filter_genes pass,
# so a cell cycle gene expressed only in the non-malignant epithelium is gone. Report the
# survivors; a much smaller set would weaken the score and has to be documented.
s_present = [g for g in s_genes if g in adata.var_names]
g2m_present = [g for g in g2m_genes if g in adata.var_names]
print(f"S genes:   {len(s_present)}/{len(s_genes)} present", flush=True)
print(f"G2M genes: {len(g2m_present)}/{len(g2m_genes)} present", flush=True)
missing = sorted(set(s_genes + g2m_genes) - set(s_present) - set(g2m_present))
if missing:
    print(f"Missing cell cycle genes (not in var_names): {missing}", flush=True)

print("Scoring cell cycle...", flush=True)
sc.tl.score_genes_cell_cycle(
    adata,
    s_genes=s_present,
    g2m_genes=g2m_present,
    random_state=C.SEED,
)

print("Cell cycle phase distribution:", flush=True)
print(adata.obs["phase"].value_counts(), flush=True)

# How the cycling fraction compares with the epithelial object of 04, when it is around.
# Not a check with a pass/fail - it is the number to quote when 05_9 asks how much of a
# 'stemness' readout is proliferation.
if C.PRIOR_LABEL_KEY in adata.obs:
    print("\nphase by pre-CNV CellTypist label:", flush=True)
    print(adata.obs.groupby(C.PRIOR_LABEL_KEY, observed=True)["phase"]
          .value_counts(normalize=True).unstack().round(3).to_string(), flush=True)

for col in ("S_score", "G2M_score", "phase"):
    assert col in adata.obs, f"scoring did not add {col}"

print("\nSaving...", flush=True)
adata.write_h5ad(OUT_PATH, compression="gzip")
print(f"Wrote {OUT_PATH}", flush=True)
