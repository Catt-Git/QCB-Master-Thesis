"""
02 utils - shared preparation used by 02_4_metrics' metrics.py

Everything the metric job does *before* calling ``scib.me.metrics`` lives here, so
every (run, type) reduces the data in exactly the same way and no two jobs can
disagree on the graph or embedding they score.

Three things happen, in order:

1. ``load_and_align`` reads the reference and the integrated object and enforces
   the property the metrics silently rely on: the two objects are the same cells,
   in the same order. scib compares them cell by cell while checking the names
   only as a *set*, so a reordering passes its own validation and corrupts every
   score. The scib-pipeline prototype papered over this by renaming ``obs_names``
   from the batch suffix; here a mismatch is a hard error instead, because on this
   benchmark the integration dispatchers already assert order and a mismatch means
   something upstream is wrong. Because the order is guaranteed identical, the
   batch and label columns are copied from the reference onto the integrated
   object, which is both correct and removes the prototype's fragile
   value-counts re-matching.

2. ``ensure_reference_reduced`` gives the reference a PCA and a neighbour graph if
   it does not already carry them. On the real unscaled reference the PCA from
   01_5 is present and is kept untouched (all 21 jobs then share one baseline);
   the scaled reference gets its fresh PCA at preparation time. The smoke fixture
   carries ``obsm['X_pca']`` but no ``uns['pca']`` (no stored variance), so it is
   recomputed here with the parameters scib itself uses.

3. ``reduce_integrated`` runs ``scib.preprocessing.reduce_data`` on the integrated
   object with the options dictated by the output type, ``umap=False`` throughout:
   no scib metric reads ``X_umap`` and the visual QC lives in 02_3.

``import scib_compat`` must precede ``import scib`` anywhere in this phase; this
module does it, so importing it first is enough for the caller too.
"""

from __future__ import annotations

import functools
import inspect
import os
import sys

import numpy as np
import scanpy as sc

# scib_compat is a sibling module in utils/; it must be imported before scib to
# restore the numpy / pandas APIs scib 1.1.7 still expects. Do it here so the
# caller only has to import this module.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scib_compat  # noqa: F401,E402  (must precede scib)
import scib  # noqa: E402

# PCA parameters scib uses internally, applied when the reference has no stored
# PCA of its own (the smoke fixture, and any future reference without one).
PCA_N_COMPS = 50
SEED = 0

# scib output types, and the metric flags each one enables.
RESULT_TYPES = ("full", "embed", "knn")

# The twelve metrics of this benchmark, in the order scib returns them, so the
# per-type "expected" lists and the smoke report speak the same language as the
# CSVs. merge_metrics.SCORED_METRICS repeats this list and must stay in step.
ALL_METRICS = (
    "NMI_cluster/label", "ARI_cluster/label", "ASW_label", "ASW_label/batch",
    "PCR_batch", "cell_cycle_conservation", "isolated_label_F1",
    "isolated_label_silhouette", "graph_conn", "iLISI", "cLISI",
    "hvg_overlap",
)

# Metrics expected to produce a finite value for each output type. Everything not
# listed is NaN *by construction* for that type and must not be read as a failure.
EXPECTED = {
    "full": [
        "NMI_cluster/label", "ARI_cluster/label", "ASW_label", "ASW_label/batch",
        "PCR_batch", "cell_cycle_conservation", "isolated_label_F1",
        "isolated_label_silhouette", "graph_conn", "iLISI", "cLISI",
        "hvg_overlap",
    ],
    "embed": [
        "NMI_cluster/label", "ARI_cluster/label", "ASW_label", "ASW_label/batch",
        "PCR_batch", "cell_cycle_conservation", "isolated_label_F1",
        "isolated_label_silhouette", "graph_conn", "iLISI", "cLISI",
    ],
    # knn (BBKNN): only the graph-based metrics; silhouette, PCR, cell cycle and
    # HVG conservation are unavailable, so isolated_label_silhouette and both ASW
    # variants fall away with them.
    "knn": [
        "NMI_cluster/label", "ARI_cluster/label", "isolated_label_F1",
        "graph_conn", "iLISI", "cLISI",
    ],
}


def use_igraph_leiden(n_iterations: int = 2, verbose: bool = False) -> bool:
    """Make scib's optimised clustering use igraph's leiden instead of leidenalg.

    On this benchmark leidenalg spends ~20 h on the ten resolutions of a
    619,693-cell graph, before a single metric is computed; igraph's
    implementation is the one scanpy itself recommends (the ``FutureWarning``
    printed in every log) and is orders of magnitude faster.

    It has to be swapped *here* rather than by pre-computing the clusterings:
    ``scib.me.metrics`` calls ``cluster_optimal_resolution(..., force=True)``,
    which reclusters every resolution unconditionally, so existing
    ``cluster_<res>`` columns are overwritten instead of reused.
    ``cluster_optimal_resolution`` defaults to ``sc.tl.leiden`` and looks that
    attribute up when it runs, so rebinding it is enough and needs no change to
    the metric call itself.

    What this changes: both flavors optimise the same objective (leidenalg's
    ``RBConfigurationVertexPartition`` *is* modularity with a resolution
    parameter), but the implementations and the stopping criterion differ
    (``n_iterations=2`` instead of iterating to convergence), so the partitions
    are similar but not identical and NMI/ARI move slightly. Nothing else does:
    the other eleven metrics read the graph or the embedding, not the clustering.
    That makes this a benchmark-wide switch - every run must be scored with the
    same flavor, or the NMI/ARI rows cannot be compared across methods.

    :returns: True if the swap happened, False on a scanpy too old to take
        ``flavor`` (the caller then runs with leidenalg, slowly but correctly).
    """
    if "flavor" not in inspect.signature(sc.tl.leiden).parameters:
        print("[cluster] WARNING: this scanpy's leiden has no `flavor` argument; "
              "falling back to leidenalg. Expect the clustering to take hours.",
              flush=True)
        return False

    original = sc.tl.leiden
    if getattr(original, "_igraph_flavor", False):  # already patched
        return True

    @functools.wraps(original)  # keeps __name__ == "leiden", which scib prints
    def leiden_igraph(adata, **kwargs):
        kwargs.setdefault("flavor", "igraph")
        kwargs.setdefault("n_iterations", n_iterations)
        kwargs.setdefault("directed", False)  # required by the igraph flavor
        return original(adata, **kwargs)

    leiden_igraph._igraph_flavor = True
    sc.tl.leiden = leiden_igraph
    if verbose:
        print(f"[cluster] leiden flavor=igraph, n_iterations={n_iterations}", flush=True)
    return True


def metric_flags(type_: str) -> dict:
    """The scib.me.metrics flags for one output type.

    Every flag matches the scib-pipeline defaults, minus HVG conservation for
    ``embed`` and minus silhouette/PCR/cell-cycle/HVG for ``knn``. Trajectory
    conservation is off for the whole benchmark.
    """
    if type_ not in RESULT_TYPES:
        raise ValueError(f"unknown output type {type_!r}, expected one of {RESULT_TYPES}")

    flags = dict(
        silhouette_=True, nmi_=True, ari_=True, pcr_=True, cell_cycle_=True,
        isolated_labels_=True, hvg_score_=True, graph_conn_=True, lisi_graph_=True,
        trajectory_=False,
    )
    if type_ == "embed":
        flags["hvg_score_"] = False
    elif type_ == "knn":
        flags["silhouette_"] = False
        flags["pcr_"] = False
        flags["cell_cycle_"] = False
        flags["hvg_score_"] = False
    return flags


def embed_key(type_: str) -> str:
    """The ``embed`` argument scib.me.metrics expects for one output type."""
    return "X_emb" if type_ == "embed" else "X_pca"


def load_and_align(uncorrected, integrated, batch_key, label_key, verbose=False):
    """Read reference and integrated object and enforce identical cell order.

    :returns: ``(adata, adata_int)`` with ``batch_key``/``label_key`` on
        ``adata_int`` guaranteed to match the reference cell for cell.
    :raises ValueError: if the two objects are not the same cells in the same
        order.
    """
    if verbose:
        print(f"[align] reference : {uncorrected}", flush=True)
    adata = sc.read_h5ad(uncorrected)
    if verbose:
        print(f"[align] integrated: {integrated}", flush=True)
    adata_int = sc.read_h5ad(integrated)

    if adata.n_obs != adata_int.n_obs:
        raise ValueError(
            f"cell count differs: reference has {adata.n_obs}, integrated has "
            f"{adata_int.n_obs}. The two objects must be the same cells."
        )

    ref_names = adata.obs_names.to_numpy()
    int_names = adata_int.obs_names.to_numpy()
    if not np.array_equal(ref_names, int_names):
        # Distinguish "different cells" from "same cells, reordered": both are
        # fatal, but they point at different upstream bugs.
        if set(ref_names) == set(int_names):
            first = int(np.argmax(ref_names != int_names))
            raise ValueError(
                "cell order differs between reference and integrated object: same "
                f"cells, permuted. First mismatch at position {first}: reference "
                f"{ref_names[first]!r} vs integrated {int_names[first]!r}. The "
                "metrics compare cells positionally; fix the integration output "
                "before scoring it (see check_integrations.py)."
            )
        missing = len(set(ref_names) - set(int_names))
        extra = len(set(int_names) - set(ref_names))
        raise ValueError(
            "reference and integrated object are not the same cells: "
            f"{missing} in reference only, {extra} in integrated only."
        )

    # Order is identical, so the reference's batch/label columns apply position
    # for position. Copy them on rather than trusting whatever the integration
    # round-trip preserved (some outputs drop or recategorise them).
    for key in (batch_key, label_key):
        if key not in adata.obs:
            raise ValueError(f"reference is missing obs[{key!r}]")
        adata_int.obs[key] = adata.obs[key].to_numpy()
        adata_int.obs[key] = adata_int.obs[key].astype("category")

    return adata, adata_int


def ensure_reference_reduced(adata, verbose=False):
    """Give the reference a PCA (with stored variance) and a neighbour graph.

    A no-op when both are already present, which is the case for the real
    reference and keeps the baseline identical across all jobs. Recomputes only
    what is missing, so the smoke fixture (X_pca but no ``uns['pca']``) gets a
    proper PCA here.
    """
    has_pca = "X_pca" in adata.obsm and "pca" in adata.uns
    if not has_pca:
        if verbose:
            print(f"[ref] recomputing PCA ({PCA_N_COMPS} comps, arpack)", flush=True)
        sc.pp.pca(adata, n_comps=PCA_N_COMPS, svd_solver="arpack", random_state=SEED)

    has_graph = "neighbors" in adata.uns and "connectivities" in adata.obsp
    if not has_graph:
        if verbose:
            print("[ref] recomputing neighbour graph", flush=True)
        sc.pp.neighbors(adata, random_state=SEED)

    return adata


def reduce_integrated(adata_int, type_, n_hvgs=None, verbose=False):
    """Run scib.preprocessing.reduce_data on the integrated object by type.

    ``umap=False`` throughout. ``full``/``embed`` recompute a PCA and a graph on
    the relevant representation; ``knn`` leaves the corrected graph (in ``obsp``)
    as the method produced it.
    """
    embed = embed_key(type_)
    if type_ == "knn":
        neighbors, pca = False, False
    else:
        neighbors, pca = True, True

    if verbose:
        print(f"[reduce] type={type_} use_rep={embed} neighbors={neighbors} "
              f"pca={pca} n_top_genes={n_hvgs}", flush=True)

    scib.preprocessing.reduce_data(
        adata_int,
        n_top_genes=n_hvgs,
        neighbors=neighbors,
        use_rep=embed,
        pca=pca,
        umap=False,
    )
    return adata_int
