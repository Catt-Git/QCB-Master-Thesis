"""
03_1 scran normalization of the non-immune subset (mirrors 01_3).

The size factors written in phase 01 were estimated on the full 619,693-cell
object and are meaningless here: scran pools cells to deconvolve per-cell
factors, so removing 443,083 immune cells changes every estimate. The subset is
therefore re-normalized from the raw counts restored by `subset_and_qc.ipynb`,
with the same parameters as 01_3 so the two objects stay comparable.

Unlike 01_3 this runs locally: 176k cells instead of 620k, no SLURM wrapper.

Input : $DATA_DIR/03_nonimm/shiao_nonimm_raw.h5ad
        Output of subset_and_qc.ipynb. .X = raw counts (int), identical to
        .layers['counts']; zero-count genes already removed.
Output: $DATA_DIR/03_nonimm/shiao_nonimm_norm.h5ad
        .X = scran log-normalized, float32; .layers['counts'] unchanged;
        .obs['size_factors'] added by scib; .raw explicitly cleared.

Usage:
This is step 1 of 4 (`norm`) of subsetting_all.sh, which is the intended way to
run it: the wrapper chains the four scripts in order, skips the steps whose
output already exists, and tees everything to 03_1_subsetting/logs/.

export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
./subsetting_all.sh            # the whole chain, resuming
./subsetting_all.sh norm       # only this step, through the wrapper

Standalone, outside the chain (no logging, no resume, DATA_DIR must be
exported by hand):
python3 scran_norm_nonimm.py
"""

from __future__ import annotations
import os
import numpy as np
import scanpy as sc
import anndata as ad
import scipy.sparse as _sp

import anndata2ri
from rpy2.robjects import conversion, default_converter

if not hasattr(anndata2ri, "activate"):
    _a2r_converter = default_converter + anndata2ri.converter

    def _activate():
        conversion.set_conversion(_a2r_converter)

    def _deactivate():
        conversion.set_conversion(default_converter)

    anndata2ri.activate = _activate
    anndata2ri.deactivate = _deactivate

# scib 1.1.7's preprocessing.normalize divides adata.X by the scran size factors via
# `sparse_matrix.multiply(dense_vector)`. scipy.sparse.multiply always returns a
# coo_matrix regardless of input format, and current anndata (>=0.11) rejects any
# .X that isn't strictly CSR/CSC ("Only CSR and CSC matrices are supported"), which
# crashes scib on modern anndata/scipy. Coerce multiply's output back to CSR on the
# concrete classes (not the private base class, which has moved across scipy versions).
for _cls in (_sp.csr_matrix, _sp.csc_matrix):
    if not getattr(_cls.multiply, "_csr_coerced", False):
        _orig_multiply = _cls.multiply

        def _multiply_tocsr(self, other, _orig=_orig_multiply):
            result = _orig(self, other)
            return result.tocsr() if _sp.issparse(result) else result

        _multiply_tocsr._csr_coerced = True
        _cls.multiply = _multiply_tocsr

import scib
import random
import matplotlib
matplotlib.use("Agg")  # headless backend: figures are written to disk, never shown
import matplotlib.pyplot as plt
import seaborn as sns

# Seed all RNGs read by scib's leiden preclustering. scib.preprocessing.normalize
# exposes no random_state argument, so we seed the global state instead.
# Same SEED as 01_3 so the two normalizations differ only by the input cells.
SEED = 0
random.seed(SEED)
np.random.seed(SEED)
sc.settings.verbosity = 1

DATA_DIR = os.environ["DATA_DIR"]
NONIMM_DIR = os.path.join(DATA_DIR, "03_nonimm")
IN_PATH = os.path.join(NONIMM_DIR, "shiao_nonimm_raw.h5ad")
OUT_PATH = os.path.join(NONIMM_DIR, "shiao_nonimm_norm.h5ad")

print("Loading data...", flush=True)
adata = sc.read_h5ad(IN_PATH)
print(adata, flush=True)

# Sanity checks on the input object (do not silently proceed if violated)
assert "counts" in adata.layers, "Expected raw counts in .layers['counts']"
# .X must be raw integer counts at this point (scran normalizes from raw). This is
# the check that catches the classic mistake of re-normalizing an already
# scran-normalized subset: shiao.h5ad carries log-normalized values in .X, and
# subset_and_qc.ipynb is the step that restores the counts.
x_sample = adata.X[:100].data if hasattr(adata.X[:100], "data") else np.asarray(adata.X[:100]).ravel()
assert np.allclose(x_sample, np.round(x_sample)), (
    ".X does not look like raw integer counts; subset_and_qc.ipynb must restore "
    ".X = .layers['counts'] before this step"
)

# The subset was re-filtered with filter_genes upstream; scran's computeSumFactors
# fails on all-zero genes, so assert rather than refilter (a silent refilter here
# would make the gene set of the saved object disagree with the notebook's report).
n_zero_genes = int(np.asarray((adata.X.sum(axis=0) == 0)).sum())
assert n_zero_genes == 0, f"{n_zero_genes} zero-count genes present; run sc.pp.filter_genes(min_cells=1) first"

# Old size factors must not survive the subset: they were estimated on the full object.
assert "size_factors" not in adata.obs, (
    "obs['size_factors'] from phase 01 is still present; it was estimated on the "
    "full 620k-cell object and must be dropped by subset_and_qc.ipynb"
)

print("Normalizing data (scran via scib)...", flush=True)
# scib scran wrapper: quickCluster -> computeSumFactors -> logNormCounts.
#   min_mean=0.1     : computeSumFactors min.mean; standard scib default for droplet data
#   log=True         : return log1p-transformed normalized expression in .X
#   precluster=True  : run quickCluster (leiden) before size-factor estimation
#   cluster_method='leiden' : leiden preclustering (louvain deprecated)
#   sparsify=False   : keep .X storage as-is, no forced conversion
# Identical to 01_3. Raw counts remain in .layers['counts'].
scib.preprocessing.normalize(
    adata,
    min_mean=0.1,
    log=True,
    precluster=True,
    cluster_method="leiden",
    sparsify=False,
)

# scib.preprocessing.normalize() leaves two footprints we don't want in the output:
#  1. adata.raw = adata, a full second copy of the just-normalized matrix that
#     nothing downstream reads.
#  2. adata.X ends up float64 (the R size_factors vector is double precision).
# See the same comment in 01_3/scran_norm.py.
adata.raw = None
if adata.X.dtype != np.float32:
    print(f"Downcasting .X from {adata.X.dtype} to float32", flush=True)
    adata.X = adata.X.astype(np.float32)

# Post-condition checks
assert "size_factors" in adata.obs, "scib did not write size_factors"
assert np.all(np.isfinite(adata.obs["size_factors"].values)), "non-finite size factors present"
assert "counts" in adata.layers, "raw counts layer lost during normalization"
assert adata.raw is None, ".raw should not be set (see comment above)"
assert adata.X.dtype == np.float32, ".X should be float32"

print("Saving normalized data...", flush=True)
os.makedirs(NONIMM_DIR, exist_ok=True)
adata.write_h5ad(OUT_PATH, compression="gzip")
print(f"Wrote {OUT_PATH}", flush=True)


print("Generating normalization diagnostics...", flush=True)
# Figures live inside the repo as lightweight QC evidence, anchored to this
# script: 03_drvi_non_immune/figures/03_1_normalization/.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PHASE_DIR = os.path.dirname(SCRIPT_DIR)  # 03_drvi_non_immune/
FIG_DIR = os.environ.get("FIG_DIR", os.path.join(PHASE_DIR, "figures", "03_1_normalization"))
os.makedirs(FIG_DIR, exist_ok=True)

# Guard: scran size factors must be strictly positive (<= 0 breaks the log and
# flags pathological cells that survived QC). Report before plotting.
sf = adata.obs["size_factors"].to_numpy()
n_nonpos = int((sf <= 0).sum())
print(f"size_factors: min={sf.min():.4f}, max={sf.max():.4f}, n(<=0)={n_nonpos}", flush=True)
if n_nonpos > 0:
    print(f"WARNING: {n_nonpos} cells with non-positive size factors", flush=True)

# Subsample for scatter readability (plot only; data untouched)
rng = np.random.default_rng(0)
n_plot = min(30000, adata.n_obs)
idx = rng.choice(adata.n_obs, size=n_plot, replace=False)
sub = adata.obs.iloc[idx]

# 1. size_factors vs total_counts (log-log).
#    Expectation: monotone positive relation; deviations flag odd cells/clusters.
#    Colored by treatment, the only low-cardinality biological covariate available
#    at this stage (cell_type has 18 levels, cohort 34).
fig, ax = plt.subplots(figsize=(6, 6))
hue = "treatment" if "treatment" in sub.columns else None
sns.scatterplot(data=sub, x="total_counts", y="size_factors",
                hue=hue, s=4, alpha=0.3, linewidth=0, ax=ax)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Total raw counts (library size)")
ax.set_ylabel("scran size factors")
ax.set_title("Size factors vs library size (non-immune)")
ax.grid(True, which="both", ls="-", alpha=0.2)
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "scran_size_factors_vs_library_size_unintegrated.png"),
            dpi=300, bbox_inches="tight")
plt.close(fig)

# 2. Distribution of size factors: expect unimodal, positive, no extreme tails.
fig, ax = plt.subplots(figsize=(6, 4))
sns.histplot(sf, bins=50, kde=False, ax=ax)
ax.set_xlabel("scran size factors")
ax.set_ylabel("Number of cells")
ax.set_title("Distribution of scran size factors (non-immune)")
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "scran_size_factors_distribution_unintegrated.png"),
            dpi=300, bbox_inches="tight")
plt.close(fig)

print(f"Saved diagnostics to {FIG_DIR}", flush=True)
