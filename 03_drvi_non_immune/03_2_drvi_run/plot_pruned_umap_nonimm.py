#!/usr/bin/env python
"""
03_2 extra: the DRVI non-immune UMAP recomputed on the **pruned** latent space.

The phase-03 counterpart of 02_3_plot_method_umap/plot_drvi_pruned_umap.py, kept
self-contained like every other script of this phase rather than imported across
phase boundaries.

Why this exists. The DRVI paper defines a latent dimension as *vanished* when its
maximum absolute value is below 1, and Supplemental Note 7 assumes the vanished
dimensions are pruned before anything else is evaluated. `run_drvi_nonimm.py`
does not prune: it builds its neighbour graph with `use_rep='X'` on all 64
dimensions, so every UMAP in figures/03_2_drvi_nonimm_64/ is an unpruned one.
This script is the hygiene check.

It touches nothing that already exists: it reads the latent object, writes one new
figure, and caches its layout as a .npy under DATA_DIR rather than back into any
.h5ad.

**Read the output against the control, not against the existing UMAPs.** UMAP's
approximate neighbour search and its spectral initialisation are sensitive enough
that the same data laid out through a slightly different code path lands in a
globally different arrangement, at the same seed and on the same cells in the same
order. `--no-prune` re-runs this identical code path on all 64 dimensions; that
control is the only image the pruned one is comparable to. The numbers the script
prints are the actual evidence - the figures illustrate it.

Usage:
  export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
  python plot_pruned_umap_nonimm.py --report-only   # the numbers, in seconds
  python plot_pruned_umap_nonimm.py                 # + the figure
  python plot_pruned_umap_nonimm.py --no-prune      # the control
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import h5py
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import scanpy as sc

# The paper's definition: vanished <=> max |z| < 1. run_drvi_nonimm.py called
# set_latent_dimension_stats with a threshold of 0.5, so the criterion is
# re-applied here rather than trusted.
VANISHED_THRESHOLD = 1.0

SEED = 0
N_NEIGHBORS = 15


def parse_args():
    p = argparse.ArgumentParser(description="non-immune DRVI UMAP, vanished dimensions pruned")
    p.add_argument("-n", "--n-latent", type=int, default=64,
                   help="latent size, i.e. which DRVI run [default: 64]")
    p.add_argument("-b", "--batch-key", default="cohort")
    p.add_argument("-l", "--label-key", default="cell_type")
    p.add_argument("--threshold", type=float, default=VANISHED_THRESHOLD,
                   help=f"max |z| below which a dimension is vanished [default: {VANISHED_THRESHOLD}]")
    p.add_argument("--seed", type=int, default=SEED)
    p.add_argument("--dpi", type=int, default=150)
    p.add_argument("--knn-sample", type=int, default=20000,
                   help="cells sampled for the kNN-preservation check (0 to skip)")
    p.add_argument("--recompute-umap", action="store_true",
                   help="ignore the cached .npy layout and lay it out again")
    p.add_argument("--report-only", action="store_true",
                   help="print the pruning report and stop, without any UMAP")
    p.add_argument("--no-prune", action="store_true",
                   help="the control: same code path and seed, NO pruning")
    return p.parse_args()


def reference_colors(path, keys):
    """{key: {category: hex}} lifted from the reference with h5py.

    Only the categories and the stored *_colors are needed, and the reference is
    a multi-GB file, so it is never opened as an AnnData. Categories without a
    stored colour fall back to scanpy's large categorical palette, in category
    order, exactly as sc.pl.umap would assign them.
    """
    out = {}
    with h5py.File(path, "r") as f:
        for key in keys:
            cats = [c.decode() if isinstance(c, bytes) else str(c)
                    for c in f["obs"][key]["categories"][:]]
            colour_key = f"{key}_colors"
            stored = []
            if colour_key in f["uns"]:
                stored = [c.decode() if isinstance(c, bytes) else str(c)
                          for c in f["uns"][colour_key][:]]
            if len(stored) == len(cats):
                out[key] = dict(zip(cats, stored))
            else:
                palette = sc.pl.palettes.default_102
                out[key] = {c: palette[i % len(palette)] for i, c in enumerate(cats)}
    return out


def draw(ax, adata, key, cmap, seed, point_size, title, legend):
    """Scatter one UMAP into ax: seeded shuffle, shared colours, optional legend.

    The point order is shuffled because the object is grouped by patient: without
    it the last cohort drawn would cover the others and true mixing would read as
    separation.
    """
    xy = np.asarray(adata.obsm["X_umap"])
    cats = list(adata.obs[key].astype("category").cat.categories)
    codes = adata.obs[key].astype("category").cat.codes.to_numpy()
    palette = sc.pl.palettes.default_102
    # One colour per category, then indexed by code: the same array backs both the
    # scatter and the legend, so a swatch cannot drift from its points.
    per_cat = np.array([cmap.get(c, palette[i % len(palette)]) for i, c in enumerate(cats)])
    colours = per_cat[codes]

    order = np.random.default_rng(seed).permutation(adata.n_obs)
    ax.scatter(xy[order, 0], xy[order, 1], c=colours[order], s=point_size,
               linewidths=0, rasterized=True)
    ax.set_title(title, fontsize=11)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("UMAP1", fontsize=8); ax.set_ylabel("UMAP2", fontsize=8)
    for spine in ax.spines.values():
        spine.set_visible(False)
    if legend:
        handles = [Line2D([0], [0], marker="o", linestyle="", markersize=5,
                          markerfacecolor=per_cat[i], markeredgewidth=0, label=str(c))
                   for i, c in enumerate(cats)]
        ncol = 1 if len(cats) <= 16 else 2
        ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.01, 0.5),
                  frameon=False, fontsize=6, ncol=ncol, handletextpad=0.2,
                  columnspacing=0.6, borderaxespad=0.0)


def pruning_report(embed, threshold, knn_sample, seed):
    """Print what the pruning removes, and return the boolean keep mask."""
    x = np.asarray(embed.X, dtype=np.float32)
    max_abs = np.abs(x).max(axis=0)
    keep = max_abs >= threshold
    n_van = int((~keep).sum())

    print(f"[prune] threshold: max |z| < {threshold} is vanished (the paper's definition)")
    print(f"[prune] {n_van} vanished / {embed.n_vars} dimensions -> {int(keep.sum())} kept")
    if "vanished" in embed.var:
        stored = embed.var["vanished"].to_numpy().astype(bool)
        same = "identical to" if np.array_equal(stored, ~keep) else "DIFFERENT from"
        print(f"[prune] the stored var['vanished'] ({int(stored.sum())} dims) is {same} this set")
    if n_van and keep.any():
        print(f"[prune] max |z|: vanished <= {max_abs[~keep].max():.4g}, "
              f"kept >= {max_abs[keep].min():.4g}")
        print(f"[prune] std:     vanished <= {x[:, ~keep].std(axis=0).max():.4g}, "
              f"kept >= {x[:, keep].std(axis=0).min():.4g}")

    var_all = x.var(axis=0)
    frac = var_all[~keep].sum() / var_all.sum() if n_van else 0.0
    print(f"[prune] the vanished dimensions carry {frac:.3e} of the total latent variance")

    if knn_sample and n_van:
        rng = np.random.default_rng(seed)
        idx = rng.choice(embed.n_obs, size=min(knn_sample, embed.n_obs), replace=False)
        a, b = x[idx], x[idx][:, keep]
        from sklearn.neighbors import NearestNeighbors
        from scipy.spatial.distance import pdist
        k = N_NEIGHBORS
        nn_a = NearestNeighbors(n_neighbors=k + 1).fit(a).kneighbors(a, return_distance=False)[:, 1:]
        nn_b = NearestNeighbors(n_neighbors=k + 1).fit(b).kneighbors(b, return_distance=False)[:, 1:]
        overlap = np.mean([len(set(u) & set(v)) / k for u, v in zip(nn_a, nn_b)])
        print(f"[prune] exact {k}-NN overlap full vs pruned on {len(idx):,} sampled cells: "
              f"{overlap:.6f} (1.0 = the graph is unchanged)")
        da, db = pdist(a[:4000]), pdist(b[:4000])
        rel = np.abs(da - db) / np.maximum(da, 1e-12)
        print(f"[prune] pairwise distances (4,000 cells): max rel dev {rel.max():.3e}, "
              f"mean {rel.mean():.3e}")
    return keep


def main():
    args = parse_args()

    data_dir = Path(os.environ.get("DATA_DIR", "")).expanduser()
    if not data_dir.is_dir():
        raise SystemExit("set DATA_DIR to the datasets root")

    run_id = f"drvi_nonimm_{args.n_latent}"
    nonimm_dir = data_dir / "03_nonimm"
    latent_h5ad = nonimm_dir / f"embed_{run_id}.h5ad"
    reference_h5ad = nonimm_dir / "shiao_nonimm_hvg_2k.h5ad"   # the DRVI input of this phase
    suffix = "allDims_control" if args.no_prune else "rmVanished"
    umap_npy = nonimm_dir / f"umap_{run_id}_{suffix}.npy"
    phase_dir = Path(__file__).resolve().parent.parent
    outdir = phase_dir / "figures" / f"03_2_{run_id}"

    print(f"[read] latent     {latent_h5ad}", flush=True)
    embed = sc.read_h5ad(latent_h5ad)
    for key in (args.batch_key, args.label_key):
        assert key in embed.obs, f"latent object missing obs[{key!r}]"

    keep = pruning_report(embed, args.threshold, args.knn_sample, args.seed)
    if args.report_only:
        return
    if args.no_prune:
        print("[control] --no-prune: keeping all dimensions", flush=True)
        keep = np.ones_like(keep)
    if not keep.any():
        raise SystemExit("every dimension is vanished at this threshold; nothing to plot")

    if umap_npy.exists() and not args.recompute_umap:
        xy = np.load(umap_npy)
        assert xy.shape == (embed.n_obs, 2), f"{umap_npy} holds {xy.shape}"
        print(f"[umap] reusing the cached layout {umap_npy}", flush=True)
    else:
        print(f"[umap] laying out {int(keep.sum())} dimensions ({embed.n_obs:,} cells)",
              flush=True)
        embed.obsm["X_pruned"] = np.asarray(embed.X[:, keep], dtype=np.float32)
        sc.pp.neighbors(embed, n_neighbors=N_NEIGHBORS, use_rep="X_pruned",
                        random_state=args.seed)
        sc.tl.umap(embed, random_state=args.seed)
        xy = np.asarray(embed.obsm["X_umap"], dtype=np.float32)
        np.save(umap_npy, xy)
        print(f"[write] {umap_npy}", flush=True)

    # Only the layout is swapped in; the object itself is never written back.
    embed.obsm["X_umap"] = xy

    cmap = reference_colors(reference_h5ad, (args.batch_key, args.label_key))

    ps = float(np.clip(12000.0 / embed.n_obs, 0.3, 8.0))
    state = "all dims (control)" if args.no_prune else "vanished pruned"
    tag = f"drvi{args.n_latent} non-immune - {int(keep.sum())}/{embed.n_vars} dims, {state}"

    outdir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    draw(axes[0], embed, args.batch_key, cmap[args.batch_key], args.seed, ps,
         f"{tag}\ncohort", legend=True)
    draw(axes[1], embed, args.label_key, cmap[args.label_key], args.seed, ps,
         f"{tag}\ncell_type", legend=True)
    fig.subplots_adjust(wspace=0.55)
    path = outdir / f"umap_cohort_celltype_{suffix}_{run_id}.png"
    fig.savefig(path, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"[write] {path}", flush=True)
    print("[done]", flush=True)


if __name__ == "__main__":
    main()
