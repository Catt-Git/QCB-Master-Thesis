#!/usr/bin/env python
"""
02_3 method plots: per-method UMAP panels for visual QC of the integrations.

Run after run_all.sh (02_2) and BEFORE the metrics (02_4_metrics): a cheap
visual check that a method actually mixed the cohorts without destroying the
cell-type structure. 

Five images per method, all rendered the same way so the two halves of every
comparison are directly comparable:

  1. {method}_cohort_integrated.png            cohort, integrated                (1 panel)
  2. {method}_cohort_int_vs_unint.png          cohort, integrated | unintegrated (2 panels)
  3. {method}_celltype_integrated.png          cell_type, integrated             (1 panel)
  4. {method}_celltype_int_vs_unint.png        cell_type, integrated | unint.    (2 panels)
  5. {method}_cohort_celltype_integrated.png   cohort + cell_type, both integr.  (2 panels)

The integrated embedding is derived from the output type:
    embed  -> neighbours on obsm['X_emb'], then UMAP
    full   -> PCA on the corrected .X, then neighbours + UMAP
    knn    -> UMAP straight on the corrected graph already in the object
The unintegrated UMAP is the reference's obsm['X_umap'] from 02_1 (the
unintegrated space); nothing is recomputed for it.

Colours are shared across every panel and method: cell_type inherits
uns['cell_type_colors'] from the reference, cohort gets a fixed large palette, so
a given cohort / cell type keeps its colour everywhere. Points are drawn in a
seeded shuffled order: the object is grouped by patient, so without shuffling the
last cohort drawn would cover the others and true mixing would read as separation.

Usage (one method):
If you want to runn all methods, refer to "plot_all.sh"
  export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
  python plot_methods_umaps.py -m harmony --type embed \
      -i $DATA_DIR/02_integration/harmony_unscaled.h5ad \
      -u $DATA_DIR/shiao_hvg_2k.h5ad \
      -o ../figures/harmony
"""

from __future__ import annotations

import argparse
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import scanpy as sc

RESULT_TYPES = ("full", "embed", "knn")


def parse_args():
    p = argparse.ArgumentParser(description="Per-method UMAP QC panels")
    p.add_argument("-i", "--integrated", required=True, help="integrated .h5ad")
    p.add_argument("-u", "--reference", required=True,
                   help="unintegrated reference .h5ad (must match the scaling variant)")
    p.add_argument("-m", "--method", required=True, help="method name, for titles/filenames")
    p.add_argument("--type", required=True, choices=RESULT_TYPES,
                   help="integration output type: full, embed or knn")
    p.add_argument("-o", "--outdir", required=True, help="directory to write the .png panels")
    p.add_argument("-b", "--batch-key", default="cohort")
    p.add_argument("-l", "--label-key", default="cell_type")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--dpi", type=int, default=150)
    return p.parse_args()


def integrated_umap(adata, type_, batch_key, seed):
    """Put a 2-D UMAP of the integrated space into obsm['X_umap'], by output type."""
    if type_ == "embed":
        assert "X_emb" in adata.obsm, "embed output has no obsm['X_emb']"
        sc.pp.neighbors(adata, use_rep="X_emb", random_state=seed)
    elif type_ == "full":
        sc.pp.pca(adata, n_comps=50, svd_solver="arpack", random_state=seed)
        sc.pp.neighbors(adata, use_rep="X_pca", random_state=seed)
    elif type_ == "knn":
        # bbknn already left a corrected neighbour graph in .obsp / .uns['neighbors'].
        assert "neighbors" in adata.uns or "connectivities" in adata.obsp, (
            "knn output carries no neighbour graph to lay out"
        )
    sc.tl.umap(adata, random_state=seed)
    return adata


def build_color_map(key, integrated, reference):
    """{category: hex} shared by both objects; cell_type reuses the reference palette."""
    cats = list(dict.fromkeys(
        list(reference.obs[key].astype("category").cat.categories)
        + list(integrated.obs[key].astype("category").cat.categories)
    ))
    uns_key = f"{key}_colors"
    ref_cats = list(reference.obs[key].astype("category").cat.categories)
    if uns_key in reference.uns and len(reference.uns[uns_key]) == len(ref_cats):
        cmap = {c: col for c, col in zip(ref_cats, reference.uns[uns_key])}
    else:
        cmap = {}
    # Fill anything missing from a large categorical palette (34 cohorts / 48 labels).
    palette = sc.pl.palettes.default_102
    missing = [c for c in cats if c not in cmap]
    for i, c in enumerate(missing):
        cmap[c] = palette[i % len(palette)]
    return cmap


def draw(ax, adata, key, cmap, seed, point_size, title, legend):
    """Scatter one UMAP into ax: seeded shuffle, shared colours, optional legend."""
    xy = np.asarray(adata.obsm["X_umap"])
    cats = list(adata.obs[key].astype("category").cat.categories)
    codes = adata.obs[key].astype("category").cat.codes.to_numpy()
    palette_arr = np.array([cmap[c] for c in cats])
    colours = palette_arr[codes]

    rng = np.random.default_rng(seed)
    order = rng.permutation(adata.n_obs)

    ax.scatter(xy[order, 0], xy[order, 1], c=colours[order], s=point_size,
               linewidths=0, rasterized=True)
    ax.set_title(title, fontsize=11)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("UMAP1", fontsize=8); ax.set_ylabel("UMAP2", fontsize=8)
    for spine in ax.spines.values():
        spine.set_visible(False)
    if legend:
        handles = [Line2D([0], [0], marker="o", linestyle="", markersize=5,
                          markerfacecolor=cmap[c], markeredgewidth=0, label=str(c))
                   for c in cats]
        ncol = 1 if len(cats) <= 16 else 2
        ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.01, 0.5),
                  frameon=False, fontsize=6, ncol=ncol, handletextpad=0.2,
                  columnspacing=0.6, borderaxespad=0.0)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    seed = args.seed

    print(f"[read] integrated  {args.integrated}", flush=True)
    integrated = sc.read_h5ad(args.integrated)
    print(f"[read] reference   {args.reference}", flush=True)
    reference = sc.read_h5ad(args.reference)

    for key in (args.batch_key, args.label_key):
        assert key in integrated.obs, f"integrated object missing obs[{key!r}]"
        assert key in reference.obs, f"reference object missing obs[{key!r}]"
    assert "X_umap" in reference.obsm, (
        "reference has no obsm['X_umap']; it should carry the unintegrated UMAP from 02_1"
    )
    # Same cells on both sides, so the comparison is honest.
    assert set(integrated.obs_names) == set(reference.obs_names), (
        "integrated and reference hold different cells"
    )

    print(f"[umap] computing integrated embedding (type={args.type})", flush=True)
    integrated = integrated_umap(integrated, args.type, args.batch_key, seed)

    cmap = {
        args.batch_key: build_color_map(args.batch_key, integrated, reference),
        args.label_key: build_color_map(args.label_key, integrated, reference),
    }

    n = integrated.n_obs
    ps = float(np.clip(12000.0 / n, 0.3, 8.0))
    m = args.method

    def save(fig, name):
        path = os.path.join(args.outdir, f"{m}_{name}.png")
        fig.savefig(path, dpi=args.dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"[write] {path}", flush=True)

    # 1. cohort, integrated (single panel)
    fig, ax = plt.subplots(figsize=(6.5, 6))
    draw(ax, integrated, args.batch_key, cmap[args.batch_key], seed, ps,
         f"{m} — cohort (integrated)", legend=True)
    save(fig, "cohort_integrated")

    # 2. cohort, integrated | unintegrated (side by side)
    fig, axes = plt.subplots(1, 2, figsize=(13, 6))
    draw(axes[0], integrated, args.batch_key, cmap[args.batch_key], seed, ps,
         f"{m} — cohort (integrated)", legend=False)
    draw(axes[1], reference, args.batch_key, cmap[args.batch_key], seed, ps,
         "cohort (unintegrated)", legend=True)
    save(fig, "cohort_int_vs_unint")

    # 3. cell_type, integrated (single panel)
    fig, ax = plt.subplots(figsize=(6.5, 6))
    draw(ax, integrated, args.label_key, cmap[args.label_key], seed, ps,
         f"{m} — cell_type (integrated)", legend=True)
    save(fig, "celltype_integrated")

    # 4. cell_type, integrated | unintegrated (side by side)
    fig, axes = plt.subplots(1, 2, figsize=(13, 6))
    draw(axes[0], integrated, args.label_key, cmap[args.label_key], seed, ps,
         f"{m} — cell_type (integrated)", legend=False)
    draw(axes[1], reference, args.label_key, cmap[args.label_key], seed, ps,
         "cell_type (unintegrated)", legend=True)
    save(fig, "celltype_int_vs_unint")

    # 5. cohort + cell_type, both integrated (side by side, no unintegrated)
    fig, axes = plt.subplots(1, 2, figsize=(13, 6))
    draw(axes[0], integrated, args.batch_key, cmap[args.batch_key], seed, ps,
         f"{m} — cohort (integrated)", legend=True)
    draw(axes[1], integrated, args.label_key, cmap[args.label_key], seed, ps,
         f"{m} — cell_type (integrated)", legend=True)
    save(fig, "cohort_celltype_integrated")

    print(f"[done] 5 panels written to {args.outdir}", flush=True)


if __name__ == "__main__":
    main()
