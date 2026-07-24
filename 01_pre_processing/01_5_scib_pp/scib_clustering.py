"""
01_5b scib clustering (optimal resolution)
Leiden clustering of the reduced (unintegrated) object at the resolution that
maximizes NMI against the CellTypist label_key. Kept separate from reduce_data
because this is the stochastic, compute-heavy step and may be rerun on its own.

Input : $DATA_DIR/all_samples_combined_scrublet_norm_cc_annotated_reduced.h5ad
        expects neighbors graph already computed (from 01_5)
Output: $DATA_DIR/shiao.h5ad
        adds per-resolution leiden columns 'optscib_unintegrated_leiden_<res>'
        and the selected 'optscib_unintegrated_leiden'

Local usage:
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
python3 scib_clustering.py
"""

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

DATA_DIR = os.environ["DATA_DIR"]
IN_PATH = os.path.join(DATA_DIR, "all_samples_combined_scrublet_norm_cc_annotated_reduced.h5ad")
OUT_PATH = os.path.join(DATA_DIR, "shiao.h5ad")  # definitive unintegrated main object

adata = sc.read_h5ad(IN_PATH)
print(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes")

# Optimal-resolution leiden. Sweeps a resolution grid and selects the value
# maximizing NMI vs label_key='cell_type'. return_all=True keeps every
# per-resolution column ('optscib_unintegrated_leiden_<res>') for the 01_6 grid.
# Sweep capped at 1.0 (default scib grid is 0.1..2.0); the low end is where this
# dataset's optimum has previously fallen.
# flavor='igraph', n_iterations=2 -> scanpy-recommended fast leiden backend.
resolutions = [r for r in scib.clustering.get_resolutions(n=20) if r <= 1.0]  # [0.1, 0.2, ..., 1.0]
scib.clustering.cluster_optimal_resolution(
    adata,
    label_key="cell_type",
    cluster_key="optscib_unintegrated_leiden",
    resolutions=resolutions,
    verbose=True,
    return_all=True,
    flavor="igraph",
    n_iterations=2,
)

adata.write_h5ad(OUT_PATH, compression="gzip")
print(f"Wrote {OUT_PATH}")