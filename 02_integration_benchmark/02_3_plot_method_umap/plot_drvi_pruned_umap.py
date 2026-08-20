#!/usr/bin/env python
"""
02_3 extra: the DRVI-128 UMAP recomputed on the **pruned** latent space.

Why this exists. The DRVI paper defines a latent dimension as *vanished* when its
maximum absolute value is below 1, and Supplemental Note 7 assumes the vanished
dimensions are pruned before anything else is evaluated. Everything drawn and
scored in this phase so far - the five 02_3 panels, the 02_4 metrics - uses the
full 128-dimensional obsm['X_emb'], vanished dimensions included. This script is
the hygiene check: same cells, same neighbour parameters, same colours, same
seeded point order as plot_methods_umaps.py, with the vanished dimensions dropped
before the neighbour graph is built.

It touches nothing that already exists:
  * reads the latent space from 02_drvi/embed_<run_id>.h5ad (never the scored
    02_integration/<run_id>.h5ad, whose cached obsm['X_umap'] backs the five
    existing panels);
  * writes one new figure, <method>_cohort_celltype_integrated_rmVanished.png;
  * caches its layout as a .npy under DATA_DIR, not back into any .h5ad.

The colours are built by plot_methods_umaps.build_color_map from the same
unintegrated reference, so a cohort or a cell type keeps the colour it has in
every other panel of the phase and the two figures can be read side by side.

Besides the figure it prints a short pruning report: how many dimensions the
paper's threshold removes, how much of the total latent variance they carry, and
how much the k-nearest-neighbour structure moves when they are dropped - which is
what actually decides whether the 02_4 metrics would change.

Usage:
  export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
  python plot_drvi_pruned_umap.py                       # drvi_unscaled_128
  python plot_drvi_pruned_umap.py --n-latent 64         # the other DRVI run
  python plot_drvi_pruned_umap.py --recompute-umap      # ignore the cached layout
  python plot_drvi_pruned_umap.py --report-only         # the numbers, no UMAP
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import anndata as ad
import scanpy as sc

from plot_methods_umaps import build_color_map, draw

# The paper's definition: vanished <=> max |z| < 1. The stored var['vanished']
# was written by set_latent_dimension_stats with a different threshold (0.5), so
# the criterion is re-applied here rather than trusted.
VANISHED_THRESHOLD = 1.0

# Matched to plot_methods_umaps.py, so the only difference between the two
# figures is the pruning itself.
SEED = 0
N_NEIGHBORS = 15


def parse_args():
    p = argparse.ArgumentParser(description="DRVI UMAP with the vanished dimensions pruned")
    p.add_argument("-n", "--n-latent", type=int, default=128,
                   help="latent size, i.e. which DRVI run [default: 128]")
    p.add_argument("-m", "--method", default="drvi", help="method name, for titles/filenames")
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
                   help="the control: same code path, same seed, NO pruning. UMAP's "
                        "layout is stochastic enough that the pruned figure is only "
                        "readable next to this one, not next to the 02_3 panel.")
    return p.parse_args()


def reference_palette_stub(path, keys):
    """A minimal AnnData carrying only what build_color_map reads from the reference.

    The unscaled reference is 815 MB and the colours need two categorical columns
    and the matching uns entries, so the file is opened with h5py and only those
    elements are lifted out. The stub has one row per category, which is all
    .cat.categories needs.
    """
    obs, uns = {}, {}
    with h5py.File(path, "r") as f:
        n = 0
        for key in keys:
            cats = [c.decode() if isinstance(c, bytes) else str(c)
                    for c in f["obs"][key]["categories"][:]]
            obs[key] = cats
            n = max(n, len(cats))
            colour_key = f"{key}_colors"
            if colour_key in f["uns"]:
                uns[colour_key] = np.array(
                    [c.decode() if isinstance(c, bytes) else str(c)
                     for c in f["uns"][colour_key][:]]
                )
    # Pad the shorter column so every category survives as a category.
    frame = pd.DataFrame(index=[str(i) for i in range(n)])
    for key, cats in obs.items():
        values = list(cats) + [cats[0]] * (n - len(cats))
        frame[key] = pd.Categorical(values, categories=cats)
    stub = ad.AnnData(np.zeros((n, 1), dtype=np.float32), obs=frame)
    stub.uns.update(uns)
    return stub


def pruning_report(embed, threshold, batch_key, label_key, knn_sample, seed):
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
    print(f"[prune] the vanished dimensions carry {frac:.3e} of the total latent variance, "
          f"i.e. {frac * 100:.4f}% of every squared euclidean distance")

    if knn_sample and n_van:
        rng = np.random.default_rng(seed)
        idx = rng.choice(embed.n_obs, size=min(knn_sample, embed.n_obs), replace=False)
        a, b = x[idx], x[idx][:, keep]
        from sklearn.neighbors import NearestNeighbors
        k = N_NEIGHBORS
        nn_a = NearestNeighbors(n_neighbors=k + 1).fit(a).kneighbors(a, return_distance=False)[:, 1:]
        nn_b = NearestNeighbors(n_neighbors=k + 1).fit(b).kneighbors(b, return_distance=False)[:, 1:]
        overlap = np.mean([len(set(u) & set(v)) / k for u, v in zip(nn_a, nn_b)])
        print(f"[prune] {k}-NN overlap full vs pruned on {len(idx):,} sampled cells: "
              f"{overlap:.4f} (1.0 = the graph scib builds is unchanged)")
    return keep


def main():
    args = parse_args()

    data_dir = Path(os.environ.get("DATA_DIR", "")).expanduser()
    if not data_dir.is_dir():
        raise SystemExit("set DATA_DIR to the datasets root")

    run_id = f"drvi_unscaled_{args.n_latent}"
    latent_h5ad = data_dir / "02_drvi" / f"embed_{run_id}.h5ad"
    reference_h5ad = data_dir / "shiao_hvg_2k.h5ad"          # the grid's `reference` for this run
    suffix = "allDims_control" if args.no_prune else "rmVanished"
    umap_npy = data_dir / "02_drvi" / f"umap_{run_id}_{suffix}.npy"
    outdir = Path(__file__).resolve().parent.parent / "figures" / run_id

    print(f"[read] latent     {latent_h5ad}", flush=True)
    embed = sc.read_h5ad(latent_h5ad)
    for key in (args.batch_key, args.label_key):
        assert key in embed.obs, f"latent object missing obs[{key!r}]"

    keep = pruning_report(embed, args.threshold, args.batch_key, args.label_key,
                          args.knn_sample, args.seed)
    if args.report_only:
        return
    if args.no_prune:
        print("[control] --no-prune: keeping all dimensions", flush=True)
        keep = np.ones_like(keep)
    if not keep.any():
        raise SystemExit("every dimension is vanished at this threshold; nothing to plot")

    if umap_npy.exists() and not args.recompute_umap:
        xy = np.load(umap_npy)
        assert xy.shape == (embed.n_obs, 2), f"{umap_npy} holds {xy.shape}, expected ({embed.n_obs}, 2)"
        print(f"[umap] reusing the cached layout {umap_npy}", flush=True)
    else:
        print(f"[umap] laying out {int(keep.sum())} pruned dimensions "
              f"({embed.n_obs:,} cells; ~25 min)", flush=True)
        embed.obsm["X_pruned"] = np.asarray(embed.X[:, keep], dtype=np.float32)
        sc.pp.neighbors(embed, n_neighbors=N_NEIGHBORS, use_rep="X_pruned",
                        random_state=args.seed)
        sc.tl.umap(embed, random_state=args.seed)
        xy = np.asarray(embed.obsm["X_umap"], dtype=np.float32)
        np.save(umap_npy, xy)
        print(f"[write] {umap_npy}", flush=True)

    # Only the layout is swapped in; the object itself is never written back.
    embed.obsm["X_umap"] = xy

    reference = reference_palette_stub(reference_h5ad, (args.batch_key, args.label_key))
    cmap = {k: build_color_map(k, embed, reference) for k in (args.batch_key, args.label_key)}

    ps = float(np.clip(12000.0 / embed.n_obs, 0.3, 8.0))
    m = args.method
    state = "all dims (control)" if args.no_prune else "vanished pruned"
    tag = f"{m}{args.n_latent} - {int(keep.sum())}/{embed.n_vars} dims, {state}"

    outdir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(13, 6))
    draw(axes[0], embed, args.batch_key, cmap[args.batch_key], args.seed, ps,
         f"{tag}\ncohort (integrated)", legend=True)
    draw(axes[1], embed, args.label_key, cmap[args.label_key], args.seed, ps,
         f"{tag}\ncell_type (integrated)", legend=True)
    path = outdir / f"{m}_cohort_celltype_integrated_{suffix}.png"
    fig.savefig(path, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"[write] {path}", flush=True)
    print("[done]", flush=True)


if __name__ == "__main__":
    main()
