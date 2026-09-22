#!/usr/bin/env python
"""
05_3b Harmony run: the same step as 05_3, with a linear integration instead of DRVI.

This is a CONTROL on 05_3, not a second main line, and it is deliberately a dead end: it
writes an embedding and figures and nothing else. 05_4 - 05_8 are not run on it and no
`shiao_tum_<run_id>.h5ad` is written by default, because what those steps do - read a latent
dimension and ask which gene programme loads on it - is a question about a DRVI decoder.
Harmony has no decoder. It rotates and shifts a PCA so the cohorts overlap, so its
dimensions are corrected principal components: they have loadings, but the loadings belong
to the PCA *before* the correction, and the correction is a different affine map per cohort.
Reading a gene programme off one of them would be reading the loadings of a basis the
coordinates are no longer expressed in. That is the whole reason this script stops at the
embedding.

What it IS for is the question 05 keeps running into: **is the structure DRVI finds in these
cells a property of the cells, or a property of DRVI?** A linear method with a different
objective, on the same cells and the same 2,000 genes, either recovers the same axes or does
not, and either answer is worth having. Concretely, three things come out of here:

  1. **The UMAPs**, drawn with the same keys, the same palettes and the same seeded
     permutation as `05_3_drvi_run/run_drvi_tum.py`, so the only thing that differs between
     figures/05_3_drvi_tum_64_nomt/ and figures/05_3b_harmony_tum_64_nomt/ is the space.
  2. **The side-by-side panels** (`umap_vs_drvi_*.png`), which put the two spaces next to
     each other on the keys that decide whether an integration worked: `cohort` (did the
     batches mix), `optscib_tum_leiden` (did the states survive the mixing), and `cnv_score`
     (is either space just drawing aneuploidy).
  3. **The two diagnostic tables**, which are the point of running this at all:

     `harmony_dim_stats.tsv`      per dimension: variance, the eta-squared of `cohort` (how
                                  much of that dimension is still batch after correction) and
                                  the Spearman correlation with sequencing depth, mito/ribo
                                  fraction and CNV burden. This is the 05 diagnosis - the one
                                  that found DR 1 of the DRVI space to be depth and DR 10 to
                                  be a real axis - run on Harmony's dimensions instead.
     `harmony_vs_drvi_corr.tsv`   the 64 x 64 Spearman matrix between these dimensions and
                                  the DRVI ones, plus the best match for each DRVI dimension.
                                  A DRVI axis that a linear method also finds is an axis of
                                  the data; one with no Harmony counterpart above the noise
                                  is either something only a non-linear decoder can see, or
                                  something DRVI made up, and the table does not settle which
                                  - it says which dimensions the question is worth asking of.

**64 components, not scanpy's 50.** Harmony has no latent size of its own: it corrects a PCA
and hands back a matrix of the same width. 64 is chosen to equal the DRVI run it is compared
against (`drvi_tum_64_nomt`), so the two spaces have the same number of coordinates and the
64 x 64 matrix above is square. It is NOT the number `02_2_integration` used - that one runs
`scib.integration.harmony`, which takes scanpy's default 50 - so this run is comparable with
05_3 and not with the phase-02 benchmark. `--n-comps 50` gets the other one.

**The input is scaled before the PCA, and that is a choice.** `.X` of the 05_2 object is
scran log-normalised counts on 2,000 HVGs. Harmony's own recipe - Seurat's, and every
harmony-pytorch example - z-scores the genes before the PCA, so that a few high-variance
genes do not own the first components; `scib.integration.harmony` does not, because the scib
benchmark controls scaling as a separate variant of the input. This script scales by default
because it is running the method as the method is meant to be run, and `--no-scale` turns it
off for a run comparable with 02_2's. The flag is in the log and nowhere in the run id, so if
both are wanted they need `--fig-dir` to be given, or the second overwrites the first.

None of this touches 05_3. Different directory, different run id (`harmony_tum_64_nomt`
against `drvi_tum_64_nomt`, built by the same `cell_set.run_id()` with `method='harmony'`),
different figure and table folders. The DRVI model, embedding and downstream object are read
where they exist and never written.

`CELL_SET` and `HVG_SET` are honoured exactly as in 05_3, through `05_2/cell_set.py`: unset
is the malignant subset on the panel without the mitochondrial genes, which is what this run
is for; `CELL_SET=epi` is the epithelial control set.

Outputs, named from the run id `harmony_<compartment>_<n_comps><hvg_tag>`:

    05_tum/embed_<run_id>.h5ad                        corrected embedding + per-dim stats
    figures/05_3b_<run_id>/                           the UMAPs and the comparison panels
    tables/05_3b_<run_id>/harmony_dim_stats.tsv       the confound diagnosis
    tables/05_3b_<run_id>/harmony_vs_drvi_corr.tsv    the 64 x 64 matrix and the best matches

`--downstream` additionally writes `shiao_tum_<run_id>.h5ad` - the 05_2 object with the
corrected space in `obsm['X_harmony']` - for anything later that wants genes beside these
coordinates. Off by default: it is ~400 MB and a few minutes of gzip, and nothing in this
phase reads it.

Resuming is the default, as everywhere in this phase: an output already on disk is reported
as [have] and reused. --overwrite recomputes everything.

Usage:
  export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
  conda run -n benchmark-py-r python3 run_harmony_tum.py                 # 64 comps, scaled
  conda run -n benchmark-py-r python3 run_harmony_tum.py --n-comps 50    # as 02_2 sizes it
  conda run -n benchmark-py-r python3 run_harmony_tum.py --no-scale      # as 02_2 inputs it
  conda run -n benchmark-py-r python3 run_harmony_tum.py --overwrite
  CELL_SET=epi conda run -n benchmark-py-r python3 run_harmony_tum.py    # the control set

The environment is `benchmark-py-r`: it is the only one holding both `harmony` (harmony-
pytorch, the implementation scib wraps) and `drvi`, and this script needs the second only to
read the DRVI embedding's var table for the comparison.
"""

from __future__ import annotations

import argparse
import os
import sys
import traceback
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless: figures are written to disk, never shown
import matplotlib.pyplot as plt  # noqa: E402

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
from scipy.stats import rankdata  # noqa: E402

# The CELL_SET -> prefix mapping and every threshold of this phase live in one module, next
# to the scripts that wrote the input. Imported rather than duplicated for the reason its own
# docstring gives: a prefix recomputed independently in each script is a silent-overwrite bug
# waiting to happen.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "05_2_subsetting"))
import cell_set as C  # noqa: E402

BATCH_KEY = C.BATCH_KEY              # 'cohort', the batch to correct, as in 05_2 and 05_3
LABEL_KEY = C.LABEL_KEY              # 'cell_type': the POST-CNV label, constant under `tum`
GROUP_KEY = C.PRIOR_LABEL_KEY        # 'cell_type_01_4': the pre-CNV label, what groups
LEIDEN_KEY = "optscib_tum_leiden"    # 05_2's clustering, as named by clustering_tum.py

# These four lists are COPIES of 05_3's, and they are copied rather than imported for one
# reason: importing run_drvi_tum pulls in scvi-tools and drvi, which is half a minute and a
# CUDA probe, to read four literals. They must stay in step with that file - the whole point
# of this run is that the only difference between the two figure folders is the space.
UMAP_QC_KEYS = {
    "cell_type": LABEL_KEY,          # constant under `tum`, skipped; the point of `epi`
    "leiden": LEIDEN_KEY,            # the primary grouping; only once attach_leiden ran
    "cell_type_01_4": GROUP_KEY,     # the secondary one, a landmark
    "cohort": BATCH_KEY,
    "treatment": "treatment",
    "response": "response",
    "phase": "phase",
    "cnv_score": "cnv_score",
    "cnv_corr": "cnv_corr",
    "n_genes_by_counts": "n_genes_by_counts",
    "total_counts": "total_counts",
    "mito": "pct_counts_mt",
    "ribo": "pct_counts_ribo",
    "size_factors": "size_factors",
}
UMAP_COMBINED_KEYS = [LABEL_KEY, LEIDEN_KEY, GROUP_KEY, BATCH_KEY, "treatment", "response",
                      "phase"]
HEATMAP_KEYS = [LEIDEN_KEY, GROUP_KEY, BATCH_KEY, "treatment", "response", "phase", LABEL_KEY]
PANEL_KEYS = [LEIDEN_KEY, GROUP_KEY]

# The keys the two spaces are put side by side on. Three questions, in this order: did the
# batches mix, did the states survive the mixing, and is either space just drawing how
# aneuploid a cell is.
COMPARE_KEYS = [BATCH_KEY, LEIDEN_KEY, GROUP_KEY, "cnv_score"]

# The continuous covariates every dimension is scored against in harmony_dim_stats.tsv.
# Depth first: it is what the 05 diagnosis found sitting on DRVI's first dimension, and a
# corrected PCA has no more defence against it than an uncorrected one - Harmony corrects
# BATCH, and sequencing depth varies within a cohort as much as between.
CONFOUND_KEYS = ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo",
                 "size_factors", "cnv_score", "cnv_corr"]

UMAP_SEED = 0


def parse_args():
    p = argparse.ArgumentParser(
        description="Run Harmony on the malignant epithelial subset, as a control on 05_3")
    p.add_argument("-n", "--n-comps", type=int, default=None,
                   help="principal components to correct; the run id follows it. Defaults "
                        "to $N_LATENT (32 when unset), so one export lines this run up with "
                        "the DRVI run of the same size")
    p.add_argument("--seed", type=int, default=0,
                   help="seed for the PCA and for harmonize [default: 0]")
    p.add_argument("--no-scale", dest="scale", action="store_false",
                   help="do not z-score the genes before the PCA. Off by default: scaling is "
                        "Harmony's own recipe, and NOT scaling is what 02_2's scib wrapper "
                        "does. See the module docstring")
    p.add_argument("--theta", type=float, default=2.0,
                   help="harmonize's diversity penalty; higher mixes the cohorts harder "
                        "[default: 2.0, harmony-pytorch's own]")
    p.add_argument("--downstream", action="store_true",
                   help="also write shiao_<set>_<run_id>.h5ad: the 05_2 object (all genes) "
                        "with the corrected space in obsm['X_harmony']. ~400 MB, a few "
                        "minutes, and nothing in this phase reads it")
    p.add_argument("--drvi-run-id", default=None,
                   help="the DRVI run to compare against [default: the same cell set, HVG "
                        "panel and size, i.e. drvi_<compartment>_<n_comps><hvg_tag>]")
    p.add_argument("--data-dir", default=os.environ.get("DATA_DIR"),
                   help="directory holding the datasets [default: $DATA_DIR]")
    p.add_argument("--fig-dir", default=None,
                   help="where the figures go [default: the repo's "
                        "05_drvi_tumoral_epi/figures/05_3b_<run_id>]")
    p.add_argument("--overwrite", action="store_true",
                   help="recompute and rewrite everything instead of reusing what is on disk")
    return p.parse_args()


def savefig(fig_dir, run_id, name, fig=None, dpi=300):
    """Save a matplotlib figure (the current one by default) into fig_dir.

    05_3's helper, verbatim: the run id is appended to the name so figures from different
    runs stay distinguishable once they are pulled out of their folder, and the figure is
    closed rather than shown because hundreds are drawn in one headless process.
    """
    fig = plt.gcf() if fig is None else fig
    path = Path(fig_dir) / f"{name}_{run_id}.png"
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {path}", flush=True)
    return path


def attach_leiden(embed, full_h5ad, key=LEIDEN_KEY):
    """Add 05_2's leiden label to `embed.obs`, reading only `obs` out of the h5ad.

    05_3's helper, verbatim, and needed here for the same reason: the 2,000-gene input was
    written by reduce_data_tum.py BEFORE clustering_tum.py ran, so the clustering is not in
    it. It is worth going and getting - `cell_type` is constant on this subset and
    `cell_type_01_4` comes from a normal-breast CellTypist model, so leiden is the only
    grouping of these cells that was computed on these cells. Any failure is reported and
    swallowed: it costs figures, not the run.
    """
    if key in embed.obs:
        return True
    full_h5ad = Path(full_h5ad)
    if not full_h5ad.exists():
        print(f"[skip] {key}: {full_h5ad.name} is not on disk", flush=True)
        return False
    try:
        import h5py
        try:
            from anndata.io import read_elem
        except ImportError:                      # anndata < 0.11
            from anndata.experimental import read_elem
        with h5py.File(full_h5ad, "r") as f:
            obs = read_elem(f["obs"])
        if key not in obs:
            print(f"[skip] {key}: not a column of {full_h5ad.name}", flush=True)
            return False
        # By name, never by position: a reordering upstream must not pair a cell with
        # another cell's cluster.
        missing = embed.obs_names.difference(obs.index)
        if len(missing):
            print(f"[skip] {key}: {len(missing):,} cells of the embedding are not in "
                  f"{full_h5ad.name}", flush=True)
            return False
        embed.obs[key] = obs.loc[embed.obs_names, key].values
        print(f"[have] {key}: {embed.obs[key].nunique()} clusters attached from "
              f"{full_h5ad.name}", flush=True)
        return True
    except Exception as exc:
        print(f"[skip] {key}: {type(exc).__name__}: {exc}", flush=True)
        return False


def groupable(embed, keys):
    """The keys of `keys` that are in `embed.obs` and have more than one level.

    05_3's helper, verbatim. A constant column is not an error here, it is the shape of this
    subset: `cell_type` is `malignant` for every cell under CELL_SET=tum.
    """
    out = []
    for key in keys:
        if key not in embed.obs:
            print(f"[skip] grouping {key!r}: not an obs column", flush=True)
        elif embed.obs[key].astype(str).nunique() < 2:
            print(f"[skip] grouping {key!r}: constant "
                  f"({embed.obs[key].astype(str).iloc[0]!r})", flush=True)
        else:
            out.append(key)
    return out


def run_harmony(adata, n_comps, batch_key, seed, scale, theta):
    """PCA, then Harmony on it. Returns the corrected matrix, (n_obs, n_comps).

    This is `scib.integration.harmony`'s body with the two things that wrapper hardcodes
    made explicit - the number of components and whether the genes are z-scored first - and
    nothing else changed: the same `harmonize` from harmony-pytorch, on the PCA of `.X`, with
    `adata.obs` as the batch frame.

    The work happens on a copy. Scaling rewrites `.X` in place and densifies it (42,096 x
    2,000 float32 is ~340 MB, which is affordable and a surprise to nobody reading a
    traceback later), and the caller's object is still needed afterwards with its log-
    normalised values intact.
    """
    from harmony import harmonize

    work = adata.copy()
    if scale:
        # zero_center=True and max_value=10 are the scanpy/Seurat defaults. The clip matters
        # on this subset: a handful of cells carry very large values for a few genes and
        # without it they own a principal component on their own.
        print(f"scaling {work.n_vars:,} genes (zero-centred, clipped at 10) ...", flush=True)
        sc.pp.scale(work, max_value=10)
    else:
        print("[note] --no-scale: the PCA runs on the scran log-normalised values as they "
              "are, which is what scib.integration.harmony does", flush=True)

    print(f"PCA -> {n_comps} components (seed {seed}) ...", flush=True)
    sc.tl.pca(work, n_comps=n_comps, svd_solver="arpack", random_state=seed)
    var_ratio = work.uns["pca"]["variance_ratio"]
    print(f"    {100 * var_ratio.sum():.1f}% of the variance in {n_comps} components "
          f"(PC1 {100 * var_ratio[0]:.1f}%)", flush=True)

    # harmonize wants the batch column as a frame; handing it the whole `obs` is what scib
    # does, but an unused categorical level anywhere in it has been known to trip the
    # grouping, so it gets exactly the one column it reads, with the empty levels dropped.
    batch_df = pd.DataFrame(
        {batch_key: pd.Categorical(work.obs[batch_key].astype(str))},
        index=work.obs_names,
    )
    n_batches = batch_df[batch_key].nunique()

    use_gpu = False
    try:
        import torch
        use_gpu = bool(torch.cuda.is_available())
    except Exception:
        pass

    print(f">>> harmonize: {work.n_obs:,} cells, {n_batches} {batch_key}s, theta={theta}, "
          f"{'GPU' if use_gpu else 'CPU'}", flush=True)
    corrected = harmonize(
        work.obsm["X_pca"],
        batch_df,
        batch_key=batch_key,
        theta=theta,
        random_state=seed,
        use_gpu=use_gpu,
    )
    corrected = np.asarray(corrected, dtype=np.float32)
    assert corrected.shape == (adata.n_obs, n_comps), corrected.shape
    assert np.isfinite(corrected).all(), "harmonize returned non-finite values"

    # The uncorrected PCA is returned too: the dim-stats table scores both, so that a
    # dimension correlating with depth can be read as "Harmony did not introduce this" or
    # "Harmony did".
    return corrected, np.asarray(work.obsm["X_pca"], dtype=np.float32), var_ratio


def _spearman_columns(a, b):
    """Spearman correlation between every column of `a` and every column of `b`.

    Rank each column once, then a Pearson correlation of the ranks - which is what Spearman
    is - as a single matrix product. `scipy.stats.spearmanr` would rank the same columns
    64 x 7 times over 42,096 cells; this is the same numbers in one pass.
    """
    def z(m):
        r = np.apply_along_axis(rankdata, 0, np.asarray(m, dtype=np.float64))
        r -= r.mean(axis=0)
        sd = r.std(axis=0)
        sd[sd == 0] = 1.0            # a constant column correlates with nothing, not NaN
        return r / sd
    za, zb = z(a), z(b)
    return (za.T @ zb) / za.shape[0]


def _eta_squared(values, groups):
    """Fraction of a dimension's variance that lies BETWEEN the levels of `groups`.

    The batch readout of this table. Harmony's whole job is to make the cohorts overlap, so
    a dimension where `cohort` still explains a large share of the variance is one the
    correction did not reach. It is a description, not a benchmark metric - the scIB
    battery is 02_4's job and is deliberately not run here.
    """
    values = np.asarray(values, dtype=np.float64)
    total = values.var() * values.size
    if total <= 0:
        return 0.0
    codes = pd.Categorical(groups).codes
    between = 0.0
    grand = values.mean()
    for c in np.unique(codes):
        v = values[codes == c]
        between += v.size * (v.mean() - grand) ** 2
    return float(between / total)


def dim_stats_table(embed, raw_pca, batch_key, confound_keys, var_ratio):
    """Per-dimension variance, batch residual and confound correlations. The point of 05_3b.

    One row per corrected dimension, and for each one:

      var, var_frac      how much of the corrected space this dimension carries. Harmony
                         returns the components in the PCA's order but does not preserve its
                         orthogonality, so these do not have to decrease monotonically and
                         they do not sum the way a PCA's do. `pca_var_ratio` beside them is
                         the variance the UNCORRECTED component had, for reference.
      eta2_<batch>       the share of this dimension still explained by the batch key AFTER
                         the correction, and `eta2_<batch>_pca` the same before it. The pair
                         is the readout: a large drop is the correction working on that
                         dimension, no drop is a dimension Harmony could not reach.
      rho_<covariate>    Spearman correlation with each of CONFOUND_KEYS, and the same on the
                         uncorrected PCA. Depth is the one to read first: the 05 diagnosis
                         found DRVI's first dimension to be sequencing depth, and Harmony has
                         no more defence against it - it corrects BATCH, and depth varies
                         inside a cohort as much as between cohorts.

    The table is sorted by dimension, not by anything interesting, so it lines up row-for-row
    with the UMAP grid drawn from the same matrix.
    """
    X = np.asarray(embed.X, dtype=np.float64)
    n_dims = X.shape[1]
    rows = {"dim": np.arange(n_dims)}

    rows["var"] = X.var(axis=0)
    rows["var_frac"] = rows["var"] / rows["var"].sum()
    rows["pca_var_ratio"] = np.asarray(var_ratio[:n_dims], dtype=np.float64)

    groups = embed.obs[batch_key].astype(str).values
    rows[f"eta2_{batch_key}"] = np.array([_eta_squared(X[:, i], groups) for i in range(n_dims)])
    rows[f"eta2_{batch_key}_pca"] = np.array(
        [_eta_squared(raw_pca[:, i], groups) for i in range(n_dims)])

    present = [k for k in confound_keys
               if k in embed.obs and pd.api.types.is_numeric_dtype(embed.obs[k])]
    missing = [k for k in confound_keys if k not in present]
    if missing:
        print(f"[skip] confounds not in obs or not numeric: {missing}", flush=True)
    if present:
        cov = embed.obs[present].to_numpy(dtype=np.float64)
        rho = _spearman_columns(X, cov)
        rho_pca = _spearman_columns(raw_pca, cov)
        for j, key in enumerate(present):
            rows[f"rho_{key}"] = rho[:, j]
            rows[f"rho_{key}_pca"] = rho_pca[:, j]

    return pd.DataFrame(rows).set_index("dim")


def drvi_comparison(embed, drvi_embed):
    """Spearman between every Harmony dimension and every DRVI dimension, plus best matches.

    Returns (matrix, per-DRVI-dimension summary). Both objects are realigned on cells BY
    NAME first: they were produced by two scripts, from two readings of the same input, and
    assuming the row order matches is the one mistake here that would produce a plausible
    wrong answer rather than a crash.

    The summary is the half to read. For each DRVI dimension it gives the Harmony dimension
    it correlates with most strongly and by how much, so the question "is DR 10 an axis of
    the data or an artefact of the decoder" has a first answer: a DRVI dimension with a
    strong linear counterpart is one a completely different method also found. The converse
    does not follow - a DRVI dimension with no match may be non-linear structure Harmony
    cannot represent, or it may be nothing - so this column narrows the question, it does not
    close it. `vanished` is carried over from the DRVI embedding so the dimensions DRVI
    itself did not use can be dropped from the reading.
    """
    shared = embed.obs_names.intersection(drvi_embed.obs_names)
    if len(shared) != embed.n_obs or len(shared) != drvi_embed.n_obs:
        print(f"[warn] the two embeddings share {len(shared):,} cells of "
              f"{embed.n_obs:,} (harmony) and {drvi_embed.n_obs:,} (drvi); "
              f"the comparison uses the intersection", flush=True)
    order = embed.obs_names[embed.obs_names.isin(shared)]
    H = np.asarray(embed[order].X, dtype=np.float64)
    D = np.asarray(drvi_embed[order].X, dtype=np.float64)

    rho = _spearman_columns(H, D)                     # (harmony dims, drvi dims)
    drvi_names = [str(v) for v in drvi_embed.var_names]
    matrix = pd.DataFrame(rho, index=[f"H{i}" for i in range(H.shape[1])],
                          columns=[f"D{n}" for n in drvi_names])

    best = np.abs(rho).argmax(axis=0)
    summary = pd.DataFrame({
        "drvi_dim": drvi_names,
        "drvi_title": drvi_embed.var["title"].astype(str).values
        if "title" in drvi_embed.var else drvi_names,
        "vanished": drvi_embed.var["vanished"].values
        if "vanished" in drvi_embed.var else False,
        "best_harmony_dim": best,
        "best_rho": rho[best, np.arange(rho.shape[1])],
        "best_abs_rho": np.abs(rho[best, np.arange(rho.shape[1])]),
    })
    return matrix, summary.sort_values("best_abs_rho", ascending=False)


def _heatmap(df, title, xlabel, ylabel, cmap="vlag", center=0.0, figsize=None):
    """One seaborn heatmap with the repo's defaults, or a matplotlib one if seaborn is absent."""
    try:
        import seaborn as sns
        fig, ax = plt.subplots(figsize=figsize or (max(6, 0.32 * df.shape[1] + 3),
                                                   max(4, 0.28 * df.shape[0] + 2)))
        sns.heatmap(df, cmap=cmap, center=center, ax=ax, cbar_kws={"shrink": 0.6},
                    linewidths=0, xticklabels=True, yticklabels=True)
    except ImportError:
        fig, ax = plt.subplots(figsize=figsize or (10, 8))
        im = ax.imshow(df.values, cmap=cmap, aspect="auto")
        ax.set_xticks(range(df.shape[1]), df.columns, rotation=90, fontsize=6)
        ax.set_yticks(range(df.shape[0]), df.index, fontsize=6)
        fig.colorbar(im, ax=ax, shrink=0.6)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.tick_params(labelsize=6)
    plt.tight_layout()
    return fig


def draw_figures(adata, embed, fig_dir, run_id, heatmap_keys, stats):
    """05_3's plotting section, with the DRVI-specific plots replaced by their equivalents.

    Everything that can be identical to 05_3 is identical - the same keys in the same order,
    the same seeded permutation so overplotting hides the same cells, the same palettes
    carried over from the input object - because the comparison between the two figure
    folders is the deliverable. What cannot be identical is the last third: DRVI's dimension
    stats, its vanished-dimension plot and its OOD/IND interpretability scores are all
    properties of a decoder Harmony does not have. In their place go the dimension grid on
    the UMAP, the dimension-by-group heatmaps and the confound heatmap, which are the three
    of those plots that mean anything for a corrected PCA.
    """
    sc.settings.figdir = fig_dir

    order = np.random.default_rng(UMAP_SEED).permutation(embed.n_obs)
    embed_plot = ad.AnnData(
        obs=embed.obs.iloc[order].copy(),
        obsm={"X_umap": embed.obsm["X_umap"][order]},
        uns={k: v for k, v in adata.uns.items() if k.endswith("_colors")},
    )

    for label, col in UMAP_QC_KEYS.items():
        if col not in embed_plot.obs:
            print(f"[skip] umap {label}: no obs column {col!r}", flush=True)
            continue
        if embed_plot.obs[col].dtype.kind not in "fiu" and \
                embed_plot.obs[col].astype(str).nunique() < 2:
            print(f"[skip] umap {label}: {col!r} is constant "
                  f"({embed_plot.obs[col].astype(str).iloc[0]!r})", flush=True)
            continue
        sc.pl.umap(embed_plot, color=col, show=False, save=f"_{label}_{run_id}.png")
        plt.close("all")

    combined = [k for k in UMAP_COMBINED_KEYS
                if k in embed_plot.obs and embed_plot.obs[k].astype(str).nunique() > 1]
    with plt.rc_context({"figure.figsize": (7, 7)}):
        sc.pl.umap(embed_plot, color=combined, ncols=2, wspace=0.8, hspace=0.25,
                   show=False, save=f"_combined_{run_id}.png")
        plt.close("all")

    for key in PANEL_KEYS:
        if key not in embed_plot.obs or embed_plot.obs[key].astype(str).nunique() < 2:
            print(f"[skip] panels {key}: missing or constant", flush=True)
            continue
        groups = embed_plot.obs[key].astype("category")
        labels = [c for c in groups.cat.categories if (groups == c).any()]
        ncols = 4
        nrows = int(np.ceil(len(labels) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4, nrows * 3.4))
        axes = np.atleast_1d(axes).flatten()
        for ax, label in zip(axes, labels):
            n = int((groups == label).sum())
            sc.pl.umap(embed_plot, color=key, groups=[label], title=f"{label} ({n:,})",
                       ax=ax, show=False, legend_loc="none", na_in_legend=False, size=4)
            ax.set_xlabel("")
            ax.set_ylabel("")
        for ax in axes[len(labels):]:
            ax.set_visible(False)
        plt.tight_layout()
        savefig(fig_dir, run_id, f"umap_per_group_{key}", fig=fig)
    del embed_plot

    # Every corrected dimension painted on the UMAP: DRVI's plot_latent_dims_in_umap, for a
    # space with no `vanished` column to filter on, so all of them are drawn.
    n_dims = embed.n_vars
    dim_ad = ad.AnnData(
        X=np.zeros((embed.n_obs, 1), dtype=np.float32),
        obs=pd.DataFrame({f"dim_{i}": np.asarray(embed.X[:, i]) for i in range(n_dims)},
                         index=embed.obs_names),
        obsm={"X_umap": embed.obsm["X_umap"]},
    )
    ncols = 8
    nrows = int(np.ceil(n_dims / ncols))
    with plt.rc_context({"figure.figsize": (2.2, 2.2)}):
        sc.pl.umap(dim_ad, color=[f"dim_{i}" for i in range(n_dims)], ncols=ncols,
                   cmap="RdBu_r", vcenter=0, frameon=False, show=False,
                   save=f"_dims_{run_id}.png")
        plt.close("all")
    del dim_ad, nrows

    # Which dimension responds to which group: DRVI's plot_latent_dims_in_heatmap, hand-rolled
    # because that helper reads var columns only a DRVI embedding has. Each dimension is
    # z-scored over all cells first, so the colour is "how far this group sits from the
    # average cell on this axis" and the dimensions are comparable across the rows.
    X = np.asarray(embed.X, dtype=np.float64)
    Z = (X - X.mean(axis=0)) / np.where(X.std(axis=0) == 0, 1.0, X.std(axis=0))
    for key in heatmap_keys:
        levels = embed.obs[key].astype(str)
        means = pd.DataFrame(Z, index=embed.obs_names).groupby(levels.values).mean()
        means.columns = [f"H{i}" for i in range(n_dims)]
        fig = _heatmap(means, f"mean z-scored dimension by {key}  ({run_id})",
                       "harmony dimension", key)
        savefig(fig_dir, run_id, f"dims_in_heatmap_{key}", fig=fig)

    # The confound panel: the diagnosis table as a picture. Corrected and uncorrected side by
    # side, so a correlation can be read as something Harmony left alone or something it
    # introduced.
    rho_cols = [c for c in stats.columns if c.startswith("rho_") and not c.endswith("_pca")]
    if rho_cols:
        block = stats[rho_cols + [f"eta2_{BATCH_KEY}"]].copy()
        block.columns = [c.replace("rho_", "") for c in rho_cols] + [f"eta2 {BATCH_KEY}"]
        block.index = [f"H{i}" for i in block.index]
        fig = _heatmap(block, f"dimension vs covariate (Spearman; last column is eta2)  "
                              f"({run_id})", "covariate", "harmony dimension",
                       figsize=(max(6, 0.7 * block.shape[1] + 3), max(5, 0.22 * n_dims + 2)))
        savefig(fig_dir, run_id, "dims_vs_confounds", fig=fig)

    # Did the correction actually move the batch out? One point per dimension, eta2 before
    # against eta2 after: everything below the diagonal is a dimension Harmony cleaned.
    if f"eta2_{BATCH_KEY}_pca" in stats:
        fig, ax = plt.subplots(figsize=(5, 5))
        ax.scatter(stats[f"eta2_{BATCH_KEY}_pca"], stats[f"eta2_{BATCH_KEY}"], s=18)
        lim = max(stats[f"eta2_{BATCH_KEY}_pca"].max(), stats[f"eta2_{BATCH_KEY}"].max()) * 1.1
        ax.plot([0, lim], [0, lim], "--", lw=1, color="grey")
        ax.set_xlim(0, lim)
        ax.set_ylim(0, lim)
        ax.set_xlabel(f"eta2({BATCH_KEY}) on the uncorrected PCA")
        ax.set_ylabel(f"eta2({BATCH_KEY}) after harmony")
        ax.set_title(f"batch variance removed per dimension  ({run_id})")
        plt.tight_layout()
        savefig(fig_dir, run_id, "eta2_before_after", fig=fig)


def draw_comparison(embed, drvi_embed, matrix, summary, fig_dir, run_id, keys):
    """The side-by-side panels and the Harmony-vs-DRVI correlation heatmap.

    The UMAPs are the two spaces on the same key, same palette, same seeded permutation, in
    one figure - which is the only honest way to look at two UMAPs, because a UMAP is not
    comparable to another one drawn separately and at a different scale. The keys are
    COMPARE_KEYS: cohort (did the batches mix), leiden and cell_type_01_4 (did the states
    survive it), cnv_score (is either space just drawing aneuploidy).
    """
    rng = np.random.default_rng(UMAP_SEED)
    order = rng.permutation(embed.n_obs)
    shared = embed.obs_names[order]

    left = ad.AnnData(obs=embed.obs.loc[shared].copy(),
                      obsm={"X_umap": embed[shared].obsm["X_umap"]})
    # The DRVI embedding is realigned on the SAME cells in the SAME order, by name: the two
    # panels must be the same cells, or the comparison is between two different pictures.
    d = drvi_embed[shared]
    right = ad.AnnData(obs=left.obs.copy(), obsm={"X_umap": np.asarray(d.obsm["X_umap"])})

    for key in keys:
        if key not in left.obs:
            print(f"[skip] comparison {key}: not an obs column", flush=True)
            continue
        is_num = pd.api.types.is_numeric_dtype(left.obs[key])
        if not is_num and left.obs[key].astype(str).nunique() < 2:
            print(f"[skip] comparison {key}: constant", flush=True)
            continue
        fig, axes = plt.subplots(1, 2, figsize=(15, 6.5))
        sc.pl.umap(left, color=key, ax=axes[0], show=False, size=4,
                   title=f"harmony  -  {key}", legend_loc="right margin" if not is_num else None)
        sc.pl.umap(right, color=key, ax=axes[1], show=False, size=4,
                   title=f"drvi  -  {key}", legend_loc="right margin" if not is_num else None)
        plt.tight_layout()
        savefig(fig_dir, run_id, f"umap_vs_drvi_{key}", fig=fig)

    # |Spearman| rather than the signed value: the sign of a latent dimension is arbitrary in
    # both methods, so it carries no information and a signed map would just be harder to
    # read. The vanished DRVI dimensions are kept in the picture, labelled, because their
    # column being empty is itself the expected result and worth seeing.
    fig = _heatmap(matrix.abs(), f"|Spearman| harmony vs drvi dimensions  ({run_id})",
                   "drvi dimension", "harmony dimension", cmap="magma", center=None)
    savefig(fig_dir, run_id, "corr_vs_drvi", fig=fig)

    # The summary as a picture: how well each DRVI dimension is matched by SOME linear
    # dimension, in descending order, with the vanished ones marked.
    s = summary.copy()
    fig, ax = plt.subplots(figsize=(max(6, 0.22 * len(s) + 2), 4))
    colors = ["lightgrey" if v else "tab:blue" for v in s["vanished"]]
    ax.bar(range(len(s)), s["best_abs_rho"], color=colors)
    ax.set_xticks(range(len(s)), s["drvi_title"], rotation=90, fontsize=5)
    ax.set_ylabel("|Spearman| with its best harmony dimension")
    ax.set_xlabel("drvi dimension (grey = vanished in DRVI)")
    ax.set_title(f"how much of each DRVI axis a linear method also finds  ({run_id})")
    plt.tight_layout()
    savefig(fig_dir, run_id, "best_match_per_drvi_dim", fig=fig)


def main():
    args = parse_args()

    if not args.data_dir:
        sys.exit("set DATA_DIR (or pass --data-dir) to the directory holding the datasets")
    os.environ["DATA_DIR"] = str(args.data_dir)   # cell_set reads it from the environment

    C.banner("05_3b Harmony run")

    n_comps = args.n_comps or C.n_latent()
    # Built by cell_set.run_id() with method='harmony', which is the only place the string
    # exists. `harmony_tum_64_nomt` cannot collide with 05_3's `drvi_tum_64_nomt`, with 04's
    # `drvi_epi_*` or with 02_2's whole-dataset runs.
    run_id = C.run_id(n_comps, method="harmony")
    drvi_run_id = args.drvi_run_id or C.run_id(n_comps)

    phase_dir = Path(__file__).resolve().parent.parent   # 05_drvi_tumoral_epi/
    tum_dir = C.tum_dir()                                # $DATA_DIR/05_tum

    input_h5ad = C.hvg_path(".h5ad")            # the 2,000-gene object 05_2 wrote
    full_h5ad = C.path(".h5ad")                 # the definitive 05_2 object (leiden lives here)
    drvi_embed_h5ad = tum_dir / f"embed_{drvi_run_id}.h5ad"

    embed_h5ad = tum_dir / f"embed_{run_id}.h5ad"
    downstream_h5ad = C.path(f"_{run_id}.h5ad")

    fig_dir = Path(args.fig_dir) if args.fig_dir else phase_dir / "figures" / f"05_3b_{run_id}"
    tab_dir = phase_dir / "tables" / f"05_3b_{run_id}"
    fig_dir.mkdir(parents=True, exist_ok=True)
    tab_dir.mkdir(parents=True, exist_ok=True)

    print(f"run id      {run_id}", flush=True)
    print(f"input       {input_h5ad}")
    print(f"embedding   {embed_h5ad}")
    print(f"compare to  {drvi_embed_h5ad}")
    print(f"figures     {fig_dir}")
    print(f"tables      {tab_dir}")
    print(f"n_comps {n_comps}, scale {args.scale}, theta {args.theta}, seed {args.seed}",
          flush=True)

    if not input_h5ad.exists():
        sys.exit(f"missing {input_h5ad}: run 05_2_subsetting/subsetting_all.sh first")

    sc.settings.set_figure_params(dpi=300, facecolor="white")
    sc.settings.figdir = fig_dir

    adata = sc.read_h5ad(input_h5ad)
    print(f"{adata.n_obs:,} cells x {adata.n_vars:,} HVGs", flush=True)

    # The same sanity checks 05_3 makes on what 05_2 produced: the right cell set, no cell
    # without a CNV call. `compartment` is what tells this object from 04's and 03's.
    expected = C.compartment()
    assert adata.obs[C.COMPARTMENT_KEY].astype(str).nunique() == 1, \
        "not a single-compartment object"
    assert adata.obs[C.COMPARTMENT_KEY].astype(str).iloc[0] == expected, \
        f"compartment is not {expected!r}: this is not the CELL_SET={C.cell_set()} object"
    assert (adata.obs[C.STATUS_KEY].astype(str) != "not_tested").all(), \
        "not_tested cells in the input: 05_2 should have dropped them"
    print("cnv_status  :", sorted(adata.obs[C.STATUS_KEY].astype(str).unique().tolist()))
    print(f"{adata.obs[BATCH_KEY].nunique()} cohorts, "
          f"{adata.obs[GROUP_KEY].nunique()} pre-CNV labels", flush=True)

    raw_pca_npy = tum_dir / f"pca_{run_id}.npy"     # kept only for the before/after columns

    if args.overwrite or not embed_h5ad.exists():
        corrected, raw_pca, var_ratio = run_harmony(
            adata, n_comps, BATCH_KEY, args.seed, args.scale, args.theta)

        embed = ad.AnnData(corrected, obs=adata.obs.copy())
        embed.var_names = [str(i) for i in range(n_comps)]
        embed.var["title"] = [f"H {i + 1}" for i in range(n_comps)]
        embed.uns["harmony"] = {
            "n_comps": n_comps, "scaled": bool(args.scale), "theta": float(args.theta),
            "seed": int(args.seed), "batch_key": BATCH_KEY,
            "pca_variance_ratio": np.asarray(var_ratio, dtype=np.float32),
        }

        print("Dimension reduction ...", flush=True)
        # Same neighbourhood as 05_3 builds on the DRVI latent space (15 neighbours, the
        # whole space as the representation), so the two UMAPs differ by the space and not
        # by how the graph over it was built.
        sc.pp.neighbors(embed, n_neighbors=15, use_rep="X", n_pcs=embed.X.shape[1])
        sc.tl.umap(embed)
        sc.pp.pca(embed)

        np.save(raw_pca_npy, raw_pca)
        print(f"[write] {embed_h5ad}", flush=True)
        embed.write_h5ad(embed_h5ad)
    else:
        print(f"[have] {embed_h5ad}, reading it back", flush=True)
        embed = sc.read_h5ad(embed_h5ad)
        var_ratio = np.asarray(embed.uns.get("harmony", {}).get(
            "pca_variance_ratio", np.full(n_comps, np.nan)))
        if raw_pca_npy.exists():
            raw_pca = np.load(raw_pca_npy)
        else:
            print(f"[warn] {raw_pca_npy.name} is gone, so the before/after columns of the "
                  f"stats table cannot be recomputed without redoing the PCA; "
                  f"recomputing it", flush=True)
            _, raw_pca, var_ratio = run_harmony(adata, n_comps, BATCH_KEY, args.seed,
                                                args.scale, args.theta)
            np.save(raw_pca_npy, raw_pca)

    # 05_2's clustering, for the heatmaps and the side-by-side panels. Not in the 2,000-gene
    # input (clustering_tum.py runs after reduce_data_tum.py), so it is fetched from the
    # definitive object; failing to get it costs figures and nothing else.
    attach_leiden(embed, full_h5ad)
    heatmap_keys = groupable(embed, HEATMAP_KEYS)

    print(">>> dimension statistics", flush=True)
    stats = dim_stats_table(embed, raw_pca, BATCH_KEY, CONFOUND_KEYS, var_ratio)
    stats_path = tab_dir / "harmony_dim_stats.tsv"
    stats.to_csv(stats_path, sep="\t", float_format="%.5f")
    print(f"[write] {stats_path}", flush=True)
    # The two numbers worth having in the log rather than only in a file.
    print(f"    eta2({BATCH_KEY}): {stats[f'eta2_{BATCH_KEY}_pca'].mean():.3f} on the PCA "
          f"-> {stats[f'eta2_{BATCH_KEY}'].mean():.3f} after harmony (mean over dimensions)",
          flush=True)
    depth_col = "rho_n_genes_by_counts"
    if depth_col in stats:
        worst = stats[depth_col].abs().idxmax()
        print(f"    depth: the dimension most correlated with n_genes_by_counts is H{worst} "
              f"(rho {stats.loc[worst, depth_col]:+.3f})", flush=True)

    print(">>> figures", flush=True)
    try:
        draw_figures(adata, embed, fig_dir, run_id, heatmap_keys, stats)
    except Exception:
        print("[warn] the figures section failed; the artifacts above are complete",
              flush=True)
        traceback.print_exc()

    # The comparison against DRVI. Optional in the strict sense - if that run is not on disk
    # this one is still a complete Harmony run - but it is why the script exists, so a
    # missing DRVI embedding is a loud skip rather than a quiet one.
    if not drvi_embed_h5ad.exists():
        print(f"[skip] {drvi_embed_h5ad.name} is not on disk, so the DRVI comparison - the "
              f"point of this step - was not drawn.\n"
              f"       Run 05_3_drvi_run/run_drvi_tum.py --n-latent {n_comps} first, or pass "
              f"--drvi-run-id for a run that is there.", flush=True)
    else:
        print(f">>> comparison against {drvi_run_id}", flush=True)
        try:
            drvi_embed = sc.read_h5ad(drvi_embed_h5ad)
            matrix, summary = drvi_comparison(embed, drvi_embed)
            matrix.to_csv(tab_dir / "harmony_vs_drvi_corr.tsv", sep="\t", float_format="%.5f")
            summary.to_csv(tab_dir / "harmony_vs_drvi_best_match.tsv", sep="\t",
                           index=False, float_format="%.5f")
            print(f"[write] {tab_dir / 'harmony_vs_drvi_corr.tsv'}", flush=True)
            print(f"[write] {tab_dir / 'harmony_vs_drvi_best_match.tsv'}", flush=True)
            live = summary[~summary["vanished"].astype(bool)] \
                if "vanished" in summary else summary
            print(f"    {len(live)} non-vanished DRVI dimensions; "
                  f"{(live['best_abs_rho'] > 0.5).sum()} have a harmony dimension at "
                  f"|rho| > 0.5, {(live['best_abs_rho'] > 0.3).sum()} at > 0.3", flush=True)
            print(live.head(10).to_string(index=False), flush=True)
            draw_comparison(embed, drvi_embed, matrix, summary, fig_dir, run_id, COMPARE_KEYS)
        except Exception:
            print("[warn] the DRVI comparison failed; the harmony run above is complete",
                  flush=True)
            traceback.print_exc()

    # Off by default: nothing in this phase reads it, and it is ~400 MB of gzip.
    if args.downstream:
        if not full_h5ad.exists():
            print(f"[warn] --downstream given but {full_h5ad} is missing; not written",
                  flush=True)
        elif not args.overwrite and downstream_h5ad.exists():
            print(f"[have] {downstream_h5ad} (--overwrite to rewrite)", flush=True)
        else:
            print(f"[read] {full_h5ad}", flush=True)
            full = sc.read_h5ad(full_h5ad)
            assert set(embed.obs_names) == set(full.obs_names), \
                "the embedding and the 05_2 object disagree on cells"
            full.obsm["X_harmony"] = np.asarray(embed[full.obs_names].X, dtype=np.float32)
            print(f"Writing {downstream_h5ad.name} (gzip, a few minutes) ...", flush=True)
            full.write_h5ad(downstream_h5ad, compression="gzip")
            print(f"[write] {downstream_h5ad}", flush=True)
            del full

    print(f"[ok] {run_id} -> {embed_h5ad}", flush=True)
    print(f"next: read {tab_dir / 'harmony_vs_drvi_best_match.tsv'} and the "
          f"umap_vs_drvi_* figures in {fig_dir}", flush=True)


if __name__ == "__main__":
    main()
