"""
03_1 leiden clustering at the optimal resolution (mirrors 01_5/scib_clustering.py).
Produces the definitive non-immune object.

Kept separate from reduce_data_nonimm.py because this is the stochastic,
compute-heavy step and may be rerun on its own.

One deliberate difference from 01_5: the resolution grid is the full scib
0.1-2.0, not capped at 1.0. The cap in 01_5 was tuned on the full object and its
48 CellTypist labels; here only 18 labels are observed but they carry finer
sub-structure (fibroblast, luminal and vascular subtypes), so the NMI optimum can
sit above 1.0 and cutting the grid there would hide it.

Input : $DATA_DIR/03_nonimm/shiao_nonimm_reduced.h5ad
        expects the neighbors graph already computed (from reduce_data_nonimm.py)
Output: $DATA_DIR/03_nonimm/shiao_nonimm.h5ad
        adds per-resolution leiden columns 'optscib_nonimm_leiden_<res>', the
        selected 'optscib_nonimm_leiden' and the resolution/NMI profile in
        .uns['optscib_nonimm_leiden_profile']
        $DATA_DIR/03_nonimm/shiao_nonimm_leiden_resolution_profile.csv
        the same profile as a table (resolution, score)

Usage:
This is step 4 of 4 (`cluster`) of subsetting_all.sh, which is the intended way
to run it: the wrapper chains the four scripts in order, skips the steps whose
output already exists, and tees everything to 03_1_subsetting/logs/. Since this
step is the one most likely to be rerun on its own, note that re-running it
needs --force, otherwise the existing shiao_nonimm.h5ad marks it as done.

export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
./subsetting_all.sh                    # the whole chain, resuming
./subsetting_all.sh cluster            # only this step, through the wrapper
./subsetting_all.sh --force cluster    # rerun it over an existing output

Standalone, outside the chain (no logging, no resume, DATA_DIR must be
exported by hand):
python3 clustering_nonimm.py
"""

from __future__ import annotations
import os
import numpy as np
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

np.random.seed(0)
sc.settings.verbosity = 1

LABEL_KEY = "cell_type"
CLUSTER_KEY = "optscib_nonimm_leiden"

DATA_DIR = os.environ["DATA_DIR"]
NONIMM_DIR = os.path.join(DATA_DIR, "03_nonimm")
IN_PATH = os.path.join(NONIMM_DIR, "shiao_nonimm_reduced.h5ad")
OUT_PATH = os.path.join(NONIMM_DIR, "shiao_nonimm.h5ad")  # definitive non-immune object
PROFILE_CSV_PATH = os.path.join(NONIMM_DIR, "shiao_nonimm_leiden_resolution_profile.csv")

adata = sc.read_h5ad(IN_PATH)
print(f"Loaded {adata.n_obs:,} cells x {adata.n_vars:,} genes", flush=True)

assert "neighbors" in adata.uns, "neighbors graph missing; run reduce_data_nonimm.py first"
assert LABEL_KEY in adata.obs, f"missing label key {LABEL_KEY!r}"
print(f"{LABEL_KEY}: {adata.obs[LABEL_KEY].nunique()} observed labels", flush=True)

# Optimal-resolution leiden. Sweeps the resolution grid and selects the value
# maximizing NMI vs label_key='cell_type'. return_all=True keeps every
# per-resolution column ('optscib_nonimm_leiden_<res>') for the visualization grid.
# Full scib grid 0.1..2.0 (see the module docstring for why it is not capped).
# flavor='igraph', n_iterations=2 -> scanpy-recommended fast leiden backend.
resolutions = scib.clustering.get_resolutions(n=20)  # [0.1, 0.2, ..., 2.0]
print(f"Resolution grid: {resolutions[0]}..{resolutions[-1]} ({len(resolutions)} values)", flush=True)
res_max, score_max, score_all = scib.clustering.cluster_optimal_resolution(
    adata,
    label_key=LABEL_KEY,
    cluster_key=CLUSTER_KEY,
    resolutions=resolutions,
    verbose=True,
    return_all=True,
    flavor="igraph",
    n_iterations=2,
)

assert CLUSTER_KEY in adata.obs, "cluster_optimal_resolution did not write the cluster column"
print(f"\nOptimal resolution: {res_max} (NMI {score_max:.4f}), "
      f"{adata.obs[CLUSTER_KEY].nunique()} clusters", flush=True)

# scib returns the resolution/NMI profile but stores nothing, so keep it: it is the
# evidence behind the chosen resolution and the visualization notebook plots it.
# Written both into the object (dict of arrays, which h5ad handles) and as a CSV.
adata.uns[f"{CLUSTER_KEY}_profile"] = {
    "resolution": score_all["resolution"].to_numpy(),
    "nmi": score_all["score"].to_numpy(),
}
score_all.to_csv(PROFILE_CSV_PATH, index=False)
print(f"Wrote the resolution profile to {PROFILE_CSV_PATH}", flush=True)

# If the optimum lands on the last grid point, the grid is the binding constraint
# rather than the data: the true optimum may sit above it.
if res_max == resolutions[-1]:
    print(f"WARNING: the optimum is the largest resolution tested ({res_max}); "
          "extend the grid (get_resolutions(n=..., max=...)) before trusting it", flush=True)

adata.write_h5ad(OUT_PATH, compression="gzip")
print(f"Wrote {OUT_PATH}", flush=True)
