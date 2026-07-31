"""
02_1_prepare: per-batch z-scoring of the HVG object.

Half the benchmark grid runs on scaled input, so this produces the scaled twin of
`shiao_hvg_2k.h5ad`. Each cohort is centred and scaled independently, which is
what scIB does and what the methods expect.

Why this code does not call scib.preprocessing.scale_batch?
`scale_batch` splits the object into 34 AnnData copies, scales each, and stitches
them back with `scib.utils.merge_adata`, which is a wrapper around
`anndata.AnnData.concatenate` - removed in anndata 0.13. The function raises
AttributeError before scaling anything, and pinning anndata back below 0.11 to
recover it would conflict with scanpy 1.12.
This code basically reimplements the same logic.

Everything inherited from 01_5/01_6 - the PCA, the neighbour graph, the UMAP -
describes the *unscaled* matrix, so all of it is dropped here. The PCA and the
UMAP are then rebuilt on the scaled matrix; the graph is not kept.

The PCA, because `scib.metrics.pcr()` reuses `obsm['X_pca']` together with
`uns['pca']['variance']` whenever both exist, and is called with
`recompute_pca=False`. Carrying the 01_5 PCA over would silently hand PCR batch
and cell-cycle conservation a baseline computed on a different matrix. It is
recomputed with the parameters scib would use internally (50 components, arpack),
so the stored result is what scib would have computed anyway - once, and
identically for every metrics job that reads it.

The UMAP, because 02_3_plot_method_umap draws the unintegrated half of its QC
panels straight from the reference's `obsm['X_umap']`, and a scaled run has to be
shown against the scaled reference: per-batch z-scoring is itself a weak batch
correction, so the unscaled layout would credit the method with mixing the
scaling had already done. No scib metric reads it (`metrics.py` runs
`reduce_data(umap=False)`), so for 02_4 it is inert either way. Parameters are
the ones 01_6 used, taken from `uns['neighbors']['params']` of the unscaled
object rather than assumed, so the two unintegrated baselines differ by the
scaling alone and not by how they were laid out.

The neighbour graph behind that UMAP is computed and then discarded: nothing
downstream reads it off the reference, and on 619,693 cells it would add hundreds
of MB to a file that is already 1.1 GB.

Note for an existing `shiao_hvg_2k_scaled.h5ad` that predates the UMAP step: do
not re-run this script to acquire it. Re-running rewrites .X, and that file is
the exact input the scaled integrations consumed - deterministic or not, it must
not be re-derived under finished runs. The UMAP was appended to it in place
instead, reading only its stored PCA. Rebuilding from scratch is for a clean
run, where no integration depends on the old file.

Input : $DATA_DIR/shiao_hvg_2k.h5ad         (619,693 x 2,000, sparse)
Output: $DATA_DIR/shiao_hvg_2k_scaled.h5ad  (dense by construction, 1.1 GB)

Local usage:
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
python3 scale_batch.py
"""

from __future__ import annotations

import os
import sys

import numpy as np
import scanpy as sc

sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils")
)
from h5ad_compat import write_h5ad_compat  # noqa: E402

DATA_DIR = os.environ["DATA_DIR"]
IN_PATH = os.path.join(DATA_DIR, "shiao_hvg_2k.h5ad")
OUT_PATH = os.path.join(DATA_DIR, "shiao_hvg_2k_scaled.h5ad")

BATCH_KEY = "cohort"
N_PCS = 50
SEED = 0

# UMAP settings, matching uns['neighbors']['params'] of the unscaled object so the
# two unintegrated baselines are laid out identically.
N_NEIGHBORS = 15
UMAP_METRIC = "euclidean"

# scanpy scales to a *sample* standard deviation of 1 (its mean/variance helper
# applies the ddof=1 correction), so a check written with the numpy default
# ddof=0 finds sqrt((n-1)/n) instead of 1 and fails on correct data.
STD_DDOF = 1
MEAN_TOL = 1e-3
STD_TOL = 1e-2


print(f"Loading {IN_PATH} ...", flush=True)
adata = sc.read_h5ad(IN_PATH)
print(f"Loaded {adata.n_obs:,} cells x {adata.n_vars:,} genes", flush=True)

assert BATCH_KEY in adata.obs, f"missing batch key {BATCH_KEY!r}"
assert "counts" in adata.layers, "expected raw counts in .layers['counts']"

obs_names_before = adata.obs_names.to_numpy().copy()
var_names_before = adata.var_names.to_numpy().copy()

# Per-batch scaling

batches = adata.obs[BATCH_KEY].cat.categories
print(f"\nScaling {len(batches)} batches of {BATCH_KEY!r} ...", flush=True)

# Preallocated output. Row positions are the only thing used to place results,
# so the cell order of the input is preserved by construction.
scaled = np.empty((adata.n_obs, adata.n_vars), dtype=np.float32)
written = np.zeros(adata.n_obs, dtype=bool)

batch_values = adata.obs[BATCH_KEY].to_numpy()
for i, batch in enumerate(batches, start=1):
    rows = np.flatnonzero(batch_values == batch)
    if rows.size == 0:
        raise ValueError(f"batch {batch!r} is an empty category; drop it upstream")

    block = adata.X[rows]
    block = block.toarray() if hasattr(block, "toarray") else np.asarray(block)
    block = block.astype(np.float32, copy=False)

    # scanpy defaults, the same ones scIB's scale_batch uses: zero_center=True
    # and no max_value clipping. Genes with zero variance in a batch come out as
    # a column of zeros, which is what sc.pp.scale does everywhere else.
    scaled[rows] = sc.pp.scale(block, zero_center=True, max_value=None, copy=True)
    written[rows] = True

    print(
        f"  [{i:2d}/{len(batches)}] {str(batch):<20} {rows.size:>7,} cells",
        flush=True,
    )

assert written.all(), (
    f"{int((~written).sum())} cells were never assigned to a batch; "
    f"check for NaN in obs[{BATCH_KEY!r}]"
)

adata.X = scaled
del scaled

# Post scaling checks

assert np.array_equal(adata.obs_names.to_numpy(), obs_names_before), (
    "cell order changed during scaling"
)
assert np.array_equal(adata.var_names.to_numpy(), var_names_before), (
    "gene order changed during scaling"
)
assert np.isfinite(adata.X).all(), "scaling produced non-finite values"

# Verify the scaling actually happened, per batch rather than globally: a global
# mean of zero is also what you would get from scaling the whole matrix at once,
# which is a different (and wrong) operation.
print("\nVerifying per-batch moments ...", flush=True)
for batch in batches:
    rows = np.flatnonzero(batch_values == batch)
    block = adata.X[rows]
    mean = np.abs(block.mean(axis=0)).max()
    # Genes constant within a batch are legitimately all-zero after scaling, so
    # their standard deviation is 0 by construction and is excluded here.
    std = block.std(axis=0, ddof=STD_DDOF)
    varying = std > 0
    std_error = np.abs(std[varying] - 1).max() if varying.any() else 0.0
    assert mean < MEAN_TOL, f"batch {batch!r}: max |mean| = {mean:.2e}"
    assert std_error < STD_TOL, f"batch {batch!r}: max |std - 1| = {std_error:.2e}"
print(f"  all {len(batches)} batches: |mean| < {MEAN_TOL}, |std - 1| < {STD_TOL}", flush=True)

# Everything derived from the unscaled matrix is now wrong for this object, so it
# all goes. The PCA and the UMAP are rebuilt on the scaled matrix further down;
# the graph is not kept. Dropping first means a rebuild that fails leaves nothing
# stale behind to be mistaken for a valid layout.
for key in ["X_umap"]:
    if key in adata.obsm:
        del adata.obsm[key]
        print(f"\nDropped obsm[{key!r}] (computed on the unscaled matrix)", flush=True)
for key in ["neighbors", "umap"]:
    if key in adata.uns:
        del adata.uns[key]
        print(f"Dropped uns[{key!r}]", flush=True)
for key in list(adata.obsp):
    del adata.obsp[key]
    print(f"Dropped obsp[{key!r}]", flush=True)

print(f"\nRecomputing PCA ({N_PCS} components, arpack) on the scaled matrix ...", flush=True)
# mask_var=None forces all 2,000 genes to be used. Without it scanpy would fall
# back to var['highly_variable'], which is all-True here, so the result is the
# same - but only by accident, and that is not a thing to leave implicit.
sc.pp.pca(adata, n_comps=N_PCS, svd_solver="arpack", random_state=SEED, mask_var=None)

assert "pca" in adata.uns and "variance" in adata.uns["pca"], "PCA metadata not written"
assert adata.obsm["X_pca"].shape == (adata.n_obs, N_PCS), "unexpected X_pca shape"

# The unintegrated UMAP of *this* object, for 02_3's QC panels.

print(f"\nComputing the unintegrated UMAP (n_neighbors={N_NEIGHBORS}, "
      f"metric={UMAP_METRIC}, on X_pca) ...", flush=True)
sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, metric=UMAP_METRIC,
                use_rep="X_pca", random_state=SEED)
sc.tl.umap(adata, random_state=SEED)

assert adata.obsm["X_umap"].shape == (adata.n_obs, 2), "unexpected X_umap shape"
assert np.isfinite(adata.obsm["X_umap"]).all(), "UMAP produced non-finite coordinates"

# Keep the layout, drop the graph that produced it: no downstream step reads a
# neighbour graph off the reference, and it would cost hundreds of MB here.
adata.uns["unintegrated_umap"] = {
    "computed_on": "X_pca",
    "n_neighbors": N_NEIGHBORS,
    "metric": UMAP_METRIC,
    "random_state": SEED,
    "note": ("unintegrated layout of the per-batch scaled matrix; the neighbour "
             "graph behind it is not persisted"),
}
for key in ["neighbors", "umap"]:
    adata.uns.pop(key, None)
for key in list(adata.obsp):
    del adata.obsp[key]
print("  layout kept in obsm['X_umap']; neighbour graph discarded", flush=True)

# Save the scaled object. The output is dense by construction (much larger).

print(f"\nWriting {OUT_PATH} ...", flush=True)
print(f"  .X is dense: {adata.n_obs:,} x {adata.n_vars:,} float32 "
      f"= {adata.X.nbytes / 1024 ** 3:.2f} GB in memory", flush=True)
# compression_opts=1 rather than the default 4: these are scaled floats, which
# barely compress, and the higher levels cost several minutes for a few percent.
write_h5ad_compat(adata, OUT_PATH, compression="gzip", compression_opts=1)
print(
    f"Done: {adata.n_obs:,} x {adata.n_vars:,}, "
    f"{os.path.getsize(OUT_PATH) / 1024 ** 3:.2f} GB on disk",
    flush=True,
)
