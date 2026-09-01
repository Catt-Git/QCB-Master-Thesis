#!/usr/bin/env python3
"""04_9: the same coordinate system question, asked of Harmony instead of DRVI.

Route A (04_5) asks where the cells that prior knowledge calls stem-like, immune-evasive or
hybrid-EMT sit along the axes of a latent space. Which space that is, is an ARGUMENT of the
question and not part of it: A1 - A5 of `cell_first_epi.py` never touch an embedding. This
step produces the alternative spaces, and `--embedding harmony` on 04_5 reads Route A off
one of them.

WHY. The phase's justification for DRVI is written in the README: the additive decoder is
what makes Route B possible, and "without it any of the higher-scoring phase-02 methods
would have done". That claim is about Route B. It says nothing about Route A, and Route A is
the route that delivers the cell assignment. Running it on Harmony is the control for that
half: if a batch-corrected PCA puts the same cells on an axis just as sharply, then Route A
is not evidence for DRVI and the whole weight of the choice sits on Route B - which is worth
knowing and worth writing down. If instead the association is CONCENTRATED on DRVI (one
dimension carrying a readout) and SMEARED across Harmony (ten PCs each carrying a bit), that
is the disentanglement claim measured on this dataset rather than cited from the paper.
`compare_embeddings_epi.py` computes exactly that contrast.

WHAT IS NOT HERE. Route B and therefore Route C. They read DRVI's additive decoder off
`embed.varm`, which no method below has: the honest linear analogue would be the PCA
loadings, and Harmony corrects the embedding and not the loadings, so a "gene programme" for
a harmonised PC would be the programme of the PC BEFORE correction. The embeddings written
here carry no `varm` at all, so 04_6 fails on them immediately instead of quietly producing
something that reads like a Route B result.

RE-RUN, NOT REUSED FROM PHASE 02. The phase-02 Harmony ran on 619,693 cells over the global
HVGs, which are dominated by immune and stromal genes; 04_1 re-did the whole pre-processing
on the epithelial subset for that reason and 04_2 retrained DRVI rather than reusing 03's.
Same argument here: Harmony is re-run on `shiao_epi_hvg_2k.h5ad`, the 2,000 batch-aware HVGs
of this compartment - the same feature set DRVI was trained on, so the two spaces differ by
the method and not by the genes.

THE CALL is `scib.integration.harmony`'s, reproduced in three lines rather than imported, the
way `02_2_integration/integration_methods.py` reproduces scVI's:

    sc.tl.pca(adata)                                            -> 50 components
    harmonize(adata.obsm['X_pca'], adata.obs, batch_key=batch)  -> the corrected embedding

Everything that affects the result is kept identical to scib - the same 50 components,
`batch_key='cohort'` as everywhere in this phase, and harmony-pytorch's own defaults for
theta, sigma, ridge_lambda and the two iteration caps. The two additions are `svd_solver`
and `random_state`, pinned so the run is reproducible; scib leaves both to scanpy's defaults
and a thesis table should not depend on that.

THE PCA ARM comes free. Harmony corrects the PCA it just computed, so writing that PCA out
uncorrected costs one extra UMAP and gives Route A a null arm: it separates "the integration
matters for this readout" from "any 50-dimensional linear space would have done".

Outputs, per embedding written, all under $DATA_DIR/04_epi/:

    embed_harmony_epi_50.h5ad   .X = cells x dimensions, .obs = the metadata,
    embed_pca_epi_50.h5ad       .var = title / order / vanished / variance_ratio

which is DRVI's embedding format, so 04_5 reads either without a branch. `vanished` is a
DRVI concept and is written all-False: a PCA dimension is never "vanished", it is only ever
low-variance, and `variance_ratio` is the honest column for that.

Figures go to figures/04_9_embedding_control/<method>/ - one folder per method, `harmony/`
and `pca/`, so a third one added later needs no reorganisation - and mirror 04_2's set, minus
the four that only exist for a decoder (the OOD/IND interpretability panels and the latent
dimension statistics): the UMAPs per metadata key, the per-cell-type panel, the dimensions on
the UMAP and the dimensions against each categorical key. Same names as 04_2's, so a Harmony
figure and its DRVI counterpart differ only by the run id in the filename. The PCA elbow
stays at the root of that folder: it describes the space BOTH arms come from.

Usage:
    export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
    cd 04_drvi_epithelial/04_9_embedding_control
    python run_harmony_epi.py                  # Harmony + the uncorrected PCA arm
    python run_harmony_epi.py --no-pca-arm     # Harmony alone
    python run_harmony_epi.py --overwrite      # recompute instead of reusing what is on disk
    python run_harmony_epi.py --n-dims 64      # parity with DRVI's 64 instead of scib's 50

Runtime: ~74k cells x 2,000 genes, CPU. The PCA is a couple of minutes, harmonize a handful,
and each UMAP a few more; the whole script is well under an hour and needs no scheduler.
Environment `benchmark-py-r`, which already carries harmony-pytorch 0.1.8 from phase 02.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")           # headless: figures are written to disk, never shown
import matplotlib.pyplot as plt  # noqa: E402

import anndata as ad            # noqa: E402
import numpy as np              # noqa: E402
import pandas as pd             # noqa: E402
import scanpy as sc             # noqa: E402
import seaborn as sns           # noqa: E402

UTILS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils")
sys.path.insert(0, UTILS_DIR)
import signature_common as C    # noqa: E402

BATCH_KEY = "cohort"            # as in phase 02, 03, 04_1 and 04_2
LABEL_KEY = "cell_type"         # CellTypist, 10 labels observed in this compartment

# The 04_2 keys, in the same order and with the same names, so the figures are comparable
# panel by panel. `compartment` and `fraction` are constant after the 04_1 subset and
# `dataset_origin` is left out as everywhere in Part 2.
UMAP_QC_KEYS = {
    "cell_type": "cell_type",
    "cohort": "cohort",
    "treatment": "treatment",
    "response": "response",
    "phase": "phase",
    "n_genes_by_counts": "n_genes_by_counts",
    "total_counts": "total_counts",
    "mito": "pct_counts_mt",
    "ribo": "pct_counts_ribo",
    "size_factors": "size_factors",
}
UMAP_COMBINED_KEYS = ["cell_type", "cohort", "treatment", "response", "phase"]
HEATMAP_KEYS = ["cell_type", "cohort", "treatment", "response", "phase"]
UMAP_SEED = 0


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--n-dims", type=int, default=C.LINEAR_N_DIMS,
                   help=f"PCA components, and so Harmony dimensions [default: {C.LINEAR_N_DIMS}, "
                        "scib's]. The run id follows it")
    p.add_argument("--seed", type=int, default=C.SEED,
                   help=f"PCA and harmonize random_state [default: {C.SEED}]")
    p.add_argument("--n-jobs", type=int, default=-1, help="harmonize n_jobs [default: -1, all cores]")
    p.add_argument("--no-pca-arm", action="store_true",
                   help="write Harmony only, without the uncorrected PCA null arm")
    p.add_argument("--no-figures", action="store_true", help="write the embeddings and stop")
    p.add_argument("--figures-only", action="store_true",
                   help="redraw the figures from the embeddings already on disk, computing "
                        "nothing: the written .h5ad carries its own UMAP, so this recovers a "
                        "deleted or restyled figure set without re-running Harmony")
    p.add_argument("--overwrite", action="store_true",
                   help="recompute and rewrite instead of reusing what is on disk")
    return p.parse_args()


def as_embedding(X, obs, prefix: str, variance_ratio=None) -> ad.AnnData:
    """Wrap a cells x dimensions matrix in DRVI's embedding format.

    The three columns 04_5 reads are `title`, `order` and `vanished`. `vanished` is False
    everywhere and is written rather than omitted: `analysis_dimensions()` then needs no
    branch, and the column being present and empty is a clearer statement than its absence -
    this space HAS no vanished dimensions, as opposed to nobody having looked.
    """
    n = X.shape[1]
    var = pd.DataFrame(
        {
            "title": [f"{prefix} {i + 1}" for i in range(n)],
            "order": np.arange(n),
            "vanished": np.zeros(n, dtype=bool),
        },
        index=[f"{prefix}_{i + 1}" for i in range(n)],
    )
    if variance_ratio is not None:
        var["variance_ratio"] = np.asarray(variance_ratio, dtype=float)[:n]
    return ad.AnnData(np.asarray(X, dtype=np.float32), obs=obs.copy(), var=var)


def draw_variance_ratio(vr, fig_root: Path) -> None:
    """The elbow of the PCA Harmony corrects. Drawn from `var['variance_ratio']` on a redraw,
    so it does not need the PCA object to still be in memory."""
    vr = np.asarray(vr, dtype=float)
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.plot(np.arange(1, len(vr) + 1), 100 * vr, "o-", ms=3, lw=1)
    ax.set_yscale("log")
    ax.set_xlabel("component")
    ax.set_ylabel("variance explained (%)")
    ax.set_title(f"PCA on the 2,000 HVGs of 04_1, {len(vr)} components\n"
                 "the space Harmony corrects", fontsize=9)
    sns.despine(ax=ax)
    fig_root.mkdir(parents=True, exist_ok=True)
    path = fig_root / f"pca_variance_ratio_epi_{len(vr)}.png"
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {path}")


def draw_figures(embed: ad.AnnData, emb: C.Embedding, fig_dir: Path, seed: int) -> None:
    """04_2's figure set, minus everything that only exists for a decoder."""
    fig_dir.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = str(fig_dir)
    sc.set_figure_params(dpi_save=300, frameon=False)

    if "X_umap" in embed.obsm:
        # The UMAP is computed once and written into the embedding, so a redraw is a redraw
        # and not a second, slightly different layout.
        print(f"\nreusing the UMAP stored in the embedding ({embed.n_obs:,} cells)", flush=True)
    else:
        print(f"\nneighbours + UMAP on the {embed.n_vars} dimensions", flush=True)
        sc.pp.neighbors(embed, n_neighbors=15, use_rep="X", n_pcs=embed.n_vars,
                        random_state=seed)
        sc.tl.umap(embed, random_state=seed)

    # Plot on a shuffled copy: drawn in file order the last cohort sits on top of every
    # other and the batch panel lies about the mixing. Same guard as run_drvi_epi.py.
    order = np.random.default_rng(UMAP_SEED).permutation(embed.n_obs)
    plot = ad.AnnData(obs=embed.obs.iloc[order].copy(),
                      obsm={"X_umap": embed.obsm["X_umap"][order]})

    for label, col in UMAP_QC_KEYS.items():
        if col not in plot.obs:
            print(f"  [skip] {col}: not in .obs")
            continue
        sc.pl.umap(plot, color=col, show=False, save=f"_{label}_{emb.run_id}.png")

    combined = [k for k in UMAP_COMBINED_KEYS if k in plot.obs]
    if combined:
        sc.pl.umap(plot, color=combined, ncols=2, wspace=0.8, show=False,
                   save=f"_combined_{emb.run_id}.png")

    # One panel per label: three of them are 81% of the compartment and own the single panel.
    cats = [c for c in plot.obs[LABEL_KEY].astype("category").cat.categories
            if (plot.obs[LABEL_KEY] == c).any()]
    ncol = 3
    nrow = int(np.ceil(len(cats) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.8 * nrow), squeeze=False)
    for a, ct in zip(axes.ravel(), cats):
        sc.pl.umap(plot, color=LABEL_KEY, groups=[ct], title=ct, ax=a, show=False,
                   legend_loc=None, size=4)
    for a in axes.ravel()[len(cats):]:
        a.axis("off")
    fig.suptitle(f"{emb.title}: one panel per cell type ({emb.run_id})", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(fig_dir / f"umap_per_cell_type_{emb.run_id}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {fig_dir / f'umap_per_cell_type_{emb.run_id}.png'}")

    # Every dimension on the UMAP, DRVI's plot_latent_dims_in_umap by hand.
    dims = embed.var["title"].tolist()
    plot.obs[dims] = pd.DataFrame(np.asarray(embed.X)[order], columns=dims,
                                  index=plot.obs_names)
    sc.pl.umap(plot, color=dims, ncols=5, show=False, save=f"_dims_{emb.run_id}.png")

    # And how the dimensions respond to each categorical key: the mean of a dimension in a
    # group, z-scored ACROSS groups so a dimension with a large scale does not own the
    # colour bar. DRVI's plot_latent_dims_in_heatmap, reimplemented for a space that has no
    # DRVI object behind it.
    L = pd.DataFrame(np.asarray(embed.X), index=embed.obs_names, columns=dims)
    for key in HEATMAP_KEYS:
        if key not in embed.obs:
            continue
        m = L.groupby(embed.obs[key].astype(str).values, observed=True).mean()
        m = (m - m.mean()) / m.std(ddof=0).replace(0, np.nan)
        m = m.fillna(0.0)
        fig, ax = plt.subplots(figsize=(0.22 * len(dims) + 4, 0.32 * len(m) + 2.5))
        sns.heatmap(m.astype(float), cmap="vlag", center=0, ax=ax,
                    cbar_kws={"label": "mean of the dimension in the group, z across groups",
                              "shrink": 0.6})
        ax.set_title(f"{emb.title} ({emb.run_id}): dimensions x {key}", fontsize=10)
        plt.setp(ax.get_xticklabels(), rotation=90, fontsize=5)
        plt.setp(ax.get_yticklabels(), rotation=0, fontsize=7)
        path = fig_dir / f"dims_in_heatmap_{key}_{emb.run_id}.png"
        fig.savefig(path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"[fig] {path}")


def main():
    args = parse_args()
    C.banner("04_9 - the coordinate systems Route A is read against, other than DRVI")

    if args.n_dims != C.LINEAR_N_DIMS:
        # The registry names the files, so a size it does not know about would be written to
        # a path 04_5 would then not find. Fail here, with the one-line fix, rather than there.
        raise SystemExit(
            f"--n-dims {args.n_dims} does not match LINEAR_N_DIMS = {C.LINEAR_N_DIMS} in "
            "utils/signature_common.py.\nThe run ids in the EMBEDDINGS registry are built "
            "from that constant; change it there and both follow.")

    fig_root = C.PHASE_DIR / "figures" / "04_9_embedding_control"
    print(f"input     {C.HVG_H5AD}")
    print(f"batch     {BATCH_KEY}")
    print(f"n_dims    {args.n_dims}   (scib's default, i.e. the phase-02 Harmony)")
    print(f"seed      {args.seed}")
    print(f"outputs   {C.EPI_DIR}")

    wanted = ["harmony"] + ([] if args.no_pca_arm else ["pca"])
    todo = [C.get_embedding(n) for n in wanted]

    if args.figures_only:
        C.banner("figures only: nothing is recomputed")
        for e in todo:
            if not e.embed_h5ad.exists():
                print(f"[skip] {e.embed_h5ad.name} is not on disk")
                continue
            embed = ad.read_h5ad(e.embed_h5ad)
            print(f"\n{e.title}: {embed.n_obs:,} x {embed.n_vars} from {e.embed_h5ad.name}")
            if "variance_ratio" in embed.var:
                draw_variance_ratio(embed.var["variance_ratio"].values, fig_root)
            draw_figures(embed, e, fig_root / e.name, args.seed)
        print("\ndone.")
        return
    missing = [e for e in todo if args.overwrite or not e.embed_h5ad.exists()]
    for e in todo:
        if e not in missing:
            print(f"[have] {e.embed_h5ad.name} (--overwrite to recompute)")
    if not missing:
        print("\nnothing to do.")
        return

    # ------------------------------------------------------------------ input
    adata = ad.read_h5ad(C.HVG_H5AD)
    print(f"\n{adata.n_obs:,} cells x {adata.n_vars:,} genes "
          "(the 2,000 batch-aware HVGs of 04_1, the DRVI training feature set)")
    assert BATCH_KEY in adata.obs, f"missing obs key: {BATCH_KEY}"
    print(f"{adata.obs[BATCH_KEY].nunique()} batches in {BATCH_KEY}")

    # Route A asserts that the embedding and the all-genes object hold the same cells in the
    # same order. Checked here, where it can still be fixed, instead of at that assert.
    full = ad.read_h5ad(C.FULL_H5AD, backed="r")
    assert (full.obs_names == adata.obs_names).all(), \
        "the HVG object and shiao_epi.h5ad are not in the same cell order"
    full.file.close()
    print(f"cell order matches {C.FULL_H5AD.name}")

    # `.X` must be the log-normalised matrix: scib runs its PCA on whatever .X holds.
    xmax = float(adata.X.max())
    assert xmax < 50, f".X looks like counts (max {xmax:.1f}); the PCA must run on log1p data"
    print(f".X max {xmax:.3f} -> log-normalised, as scib's harmony expects")

    # -------------------------------------------------------------------- PCA
    C.banner(f"PCA, {args.n_dims} components")
    sc.tl.pca(adata, n_comps=args.n_dims, svd_solver="arpack", random_state=args.seed)
    vr = adata.uns["pca"]["variance_ratio"]
    print(f"variance explained: {100 * vr.sum():.1f}% over {args.n_dims} components "
          f"(PC1 {100 * vr[0]:.1f}%, PC{args.n_dims} {100 * vr[-1]:.2f}%)")

    draw_variance_ratio(vr, fig_root)

    # ---------------------------------------------------------------- Harmony
    built = {}
    if any(e.name == "harmony" for e in missing):
        C.banner("Harmony (harmony-pytorch, scib.integration.harmony's call)")
        from harmony import harmonize
        X = harmonize(adata.obsm["X_pca"], adata.obs, batch_key=BATCH_KEY,
                      random_state=args.seed, n_jobs=args.n_jobs, use_gpu=False)
        built["harmony"] = as_embedding(X, adata.obs, C.get_embedding("harmony").dim_prefix, vr)
        shift = np.linalg.norm(np.asarray(X) - adata.obsm["X_pca"], axis=1)
        print(f"correction size: median ||harmony - pca|| = {np.median(shift):.3f} "
              f"(mean {shift.mean():.3f})")

    if any(e.name == "pca" for e in missing):
        built["pca"] = as_embedding(adata.obsm["X_pca"], adata.obs,
                                    C.get_embedding("pca").dim_prefix, vr)

    # ----------------------------------------------------------------- output
    for e in missing:
        embed = built[e.name]
        C.banner(f"{e.title}  ->  {e.embed_h5ad.name}")
        print(e.description)
        if not args.no_figures:
            draw_figures(embed, e, fig_root / e.name, args.seed)
        embed.write_h5ad(e.embed_h5ad, compression="gzip")
        print(f"\n[write] {e.embed_h5ad}  ({embed.n_obs:,} x {embed.n_vars})")

    print("\ndone. Route A on one of them:")
    for e in missing:
        print(f"  python ../04_5_cell_first/cell_first_epi.py --embedding {e.name} "
              "--collection scie")


if __name__ == "__main__":
    main()
