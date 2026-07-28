"""
02_2 integration: method implementations for the Python side.

One function per method, all with the same shape:

    method(adata, batch_key, label_key=None, seed=0) -> AnnData

so the dispatcher in run_integration.py can call them uniformly. Each returns
the object carrying whatever the metrics need for its output type:

    knn   (bbknn)                 -> corrected neighbour graph in .obsp / .uns
    embed (harmony, scvi, scanvi) -> obsm['X_emb']
    full + embed (scanorama)      -> corrected .X and obsm['X_emb']

BBKNN, Harmony, scVI and scANVI are thin wrappers over scib.integration: they
modify the object in place and preserve cell order, so nothing extra is needed.

Scanorama is reimplemented here and NOT taken from scib.integration, because
`anndata.AnnData.concatenate` was removed in anndata 0.13. This version recreates
what the scib wrapper did, but preserves the input cell order by preallocating the
output arrays and writing each batch's corrected block back to its original row.

scib_compat must be imported by the caller before scib; this module assumes it.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import scib


def bbknn(adata, batch_key, label_key=None, seed=0):
    """BBKNN: batch-balanced kNN graph. Output type: knn.

    scib.integration.bbknn recomputes PCA and returns a copy with the corrected
    neighbour graph in .obsp / .uns['neighbors']. It switches to
    neighbors_within_batch=25 above 1e5 cells, which is our regime.
    """
    return scib.integration.bbknn(adata, batch_key)


def harmony(adata, batch_key, label_key=None, seed=0):
    """Harmony (harmony-pytorch). Output type: embed.

    scib.integration.harmony imports `from harmony import harmonize`.
    It computes its own PCA on .X and writes the corrected
    embedding to obsm['X_emb']; .X is left untouched. It ignores any precomputed
    X_pca, so on the scaled variant the PCA it uses is the scaled one, as intended.
    """
    return scib.integration.harmony(adata, batch_key)


def scvi(adata, batch_key, label_key=None, seed=0):
    """scVI. Output type: embed.

    Reads raw counts from layers['counts'] (scib calls setup_anndata(layer=
    'counts')), so the input must be the unscaled object with that layer intact.
    Writes the latent space to obsm['X_emb'].
    """
    import scvi as scvi_tools

    scvi_tools.settings.seed = seed
    return scib.integration.scvi(adata, batch_key)


def scanvi(adata, batch_key, label_key=None, seed=0):
    """scANVI. Output type: embed.

    Semi-supervised: needs the cell-type labels as well as the batch. Like scVI
    it reads layers['counts'] and writes obsm['X_emb'].
    """
    if label_key is None:
        raise ValueError("scanvi requires label_key (the cell-type annotation)")

    import scvi as scvi_tools

    scvi_tools.settings.seed = seed
    return scib.integration.scanvi(adata, batch_key, label_key)


def scanorama(adata, batch_key, label_key=None, seed=0):
    """Scanorama, reimplemented to preserve cell order. Output types: full + embed.

    Corrects each batch with scanorama.correct_scanpy and writes the results back
    to preallocated arrays by original row position, so the cell order of the
    input is preserved by construction. The corrected expression becomes .X
    (full), the joint embedding becomes obsm['X_emb'] (embed).
    """
    try:
        import scanorama
    except ModuleNotFoundError as error:
        raise ModuleNotFoundError(
            "scanorama is not installed in this environment"
        ) from error

    if not isinstance(adata.obs[batch_key].dtype, pd.CategoricalDtype):
        adata.obs[batch_key] = adata.obs[batch_key].astype("category")
    batches = adata.obs[batch_key].cat.categories
    batch_values = adata.obs[batch_key].to_numpy()

    # Split by batch, remembering each batch's original row positions.
    splits = []
    row_index = []
    for batch in batches:
        rows = np.flatnonzero(batch_values == batch)
        if rows.size == 0:
            raise ValueError(f"batch {batch!r} is an empty category; drop it upstream")
        row_index.append(rows)
        splits.append(adata[rows].copy())

    # return_dimred=True writes obsm['X_scanorama'] (the joint embedding) on each
    # returned object alongside the corrected .X.
    corrected = scanorama.correct_scanpy(splits, return_dimred=True)

    emb_dim = corrected[0].obsm["X_scanorama"].shape[1]
    corrected_x = np.empty((adata.n_obs, adata.n_vars), dtype=np.float32)
    corrected_emb = np.empty((adata.n_obs, emb_dim), dtype=np.float32)
    written = np.zeros(adata.n_obs, dtype=bool)

    for rows, part in zip(row_index, corrected):
        block = part.X
        block = block.toarray() if sp.issparse(block) else np.asarray(block)
        corrected_x[rows] = block.astype(np.float32, copy=False)
        corrected_emb[rows] = part.obsm["X_scanorama"].astype(np.float32, copy=False)
        written[rows] = True

    assert written.all(), f"{int((~written).sum())} cells never corrected by scanorama"

    adata.X = corrected_x
    adata.obsm["X_emb"] = corrected_emb
    return adata


# Registry consumed by run_integration.py. label_key is only used by scanvi.
METHODS = {
    "bbknn": bbknn,
    "harmony": harmony,
    "scvi": scvi,
    "scanvi": scanvi,
    "scanorama": scanorama,
}
