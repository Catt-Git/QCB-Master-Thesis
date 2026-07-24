"""
01_4_cell_cycle: cell cycle phase scoring (Tirosh et al. 2016 signature).

Input : $DATA_DIR/all_samples_combined_scrublet_norm.h5ad
        Output of 01_3. .X = scran log-normalized (scoring requires log-norm, not raw);
        .layers['counts'] = raw counts; gene symbols as var_names.
Output: $DATA_DIR/all_samples_combined_scrublet_norm_cc.h5ad
        Adds .obs['S_score'], .obs['G2M_score'], .obs['phase'].

Requires the Tirosh/Regev cell cycle gene list (97 genes: 43 S + 54 G2M),
one gene symbol per line, path via $CC_GENES env var.

Local usage:
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
export CC_GENES=$DATA_DIR/regev_lab_cell_cycle_genes.txt
python3 cell_cycle_score.py
"""

from __future__ import annotations
import os
import numpy as np
import scanpy as sc
import anndata as ad

sc.settings.verbosity = 1

DATA_DIR = os.environ["DATA_DIR"]
IN_PATH = os.path.join(DATA_DIR, "all_samples_combined_scrublet_norm.h5ad")
OUT_PATH = os.path.join(DATA_DIR, "all_samples_combined_scrublet_norm_cc.h5ad")

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
# Report survivors: the standard file uses dated symbols (MLF1IP, FAM64A, HN1 ->
# CENPU, PIMREG, JPT1 in current HGNC), which may not match a recent reference and
# will drop silently here. Low overlap weakens the score and must be documented.
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