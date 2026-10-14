#!/usr/bin/env python3
"""05_11, Milo: does the COMPOSITION of the malignant states move with the treatment?

`treatment_programmes_tum.py` asks whether a named programme goes up or down in a biopsy.
That misses a change no single programme carries - a region of the latent space that fills
or empties between two timepoints. Milo (Dann et al. 2022) asks exactly that, without
clustering: it tiles the k-NN graph with overlapping neighbourhoods, counts each one's cells
per biopsy, and tests the counts with a negative-binomial GLM.

THE PIECES, AND WHERE EACH ONE IS DONE.

  graph       k-NN on 05_3's latent (the analysis dimensions, so PRUNE_VANISHED reaches it),
              built here with sklearn. DRVI was trained with `batch_key='cohort'`, so the
              graph is patient-corrected and the timepoint is still in it.
  nhoods      miloR's refined sampling: a random `--prop` of the cells, each replaced by the
              cell of its own neighbourhood closest to that neighbourhood's median position.
              A neighbourhood is that index cell and its k nearest neighbours.
  counts      cells per neighbourhood per biopsy (patient x timepoint, >= --min-cells).
  GLM         `milo_edger.R`, edgeR quasi-likelihood, ~ cohort + treatment: every contrast
              is within-patient. A second model, ~ cohort + response:timepoint, is written as
              exploratory - 4 to 8 patients a group.
  spatial FDR miloR's weighted BH, weights 1 / (distance of the index cell to its k-th
              neighbour), so the many overlapping neighbourhoods of a dense region do not
              count as many independent discoveries.

READING A NEIGHBOURHOOD BY THE PROGRAMMES. Every neighbourhood gets the mean of each named
programme's readout over its cells - the same oriented, z-scored composite the other script
tests - and its `top_programme`, the highest one if it clears 0.5 sd. The association table
is the Spearman of a contrast's logFC against each programme's neighbourhood mean: "the
neighbourhoods that grow after RT are the senescent ones". No p-value is attached to it, on
purpose: neighbourhoods overlap and share cells, so any p computed over them as if
independent would be too small.

WHAT IT DOES NOT SEE. Only the malignant cells are in the graph, so a neighbourhood's logFC
is its change as a SHARE of the biopsy's tumour cells, not in absolute numbers. A tumour that
shrinks uniformly moves nothing here; that is the right null for "did the states change", and
the cell count table of the other script is where the shrinkage itself is read.

Usage:
    export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
    N_LATENT=64 HVG_SET=nomt python3 milo_treatment_tum.py
    N_LATENT=64 HVG_SET=nomt python3 milo_treatment_tum.py --k 50 --prop 0.1
    RSCRIPT=/path/to/Rscript N_LATENT=64 HVG_SET=nomt python3 milo_treatment_tum.py
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr
from sklearn.neighbors import NearestNeighbors

import treatment_common as T
from treatment_common import C

R_SCRIPT = T.STEP_DIR / "milo_edger.R"

# A neighbourhood is annotated with a programme only when its mean on it clears this.
TOP_PROGRAMME_MIN = 0.5

SPATIAL_FDR_ALPHA = 0.10      # miloR's customary level for a DA neighbourhood


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--k", type=int, default=100,
                   help="neighbours per neighbourhood (default 100). Milo's rule of thumb is "
                        "a mean neighbourhood size of ~5 x the number of biopsies; with ~45 "
                        "biopsies of a few hundred cells each, a small k leaves most "
                        "neighbourhoods with 0-1 cells per biopsy. The run prints where it "
                        "lands against the rule")
    p.add_argument("--prop", type=float, default=0.05,
                   help="fraction of cells sampled as neighbourhood indices (default 0.05, "
                        "miloR's value for large data)")
    p.add_argument("--min-cells", type=int, default=T.MIN_CELLS_PER_SAMPLE,
                   help=f"malignant cells a biopsy needs to be counted "
                        f"(default {T.MIN_CELLS_PER_SAMPLE})")
    p.add_argument("--seed", type=int, default=C.SEED)
    p.add_argument("--rscript", default=os.environ.get("RSCRIPT", ""),
                   help="Rscript with edgeR installed (default: $RSCRIPT, else the Rscript "
                        "next to this Python, else the one on PATH)")
    p.add_argument("--reuse", action="store_true",
                   help="do not refit edgeR if its outputs are on disk AND were fitted on "
                        "exactly the counts built now (the neighbourhoods are rebuilt either "
                        "way: they are deterministic given --seed)")
    p.add_argument("--dpi", type=int, default=300)
    return p.parse_args()


def find_rscript(arg: str) -> str:
    """The Rscript of the active environment first: the one whose edgeR was checked."""
    if arg:
        return arg
    here = Path(sys.executable).with_name("Rscript")
    if here.exists():
        return str(here)
    found = shutil.which("Rscript")
    if not found:
        sys.exit("no Rscript found: pass --rscript or set $RSCRIPT to one with edgeR")
    return found


# --------------------------------------------------------------------------- #
# Neighbourhoods
# --------------------------------------------------------------------------- #


def build_nhoods(X: np.ndarray, k: int, prop: float, seed: int):
    """(index cells, membership nhoods x cells as CSR, k-th neighbour distance per index).

    Refined sampling as in miloR's `makeNhoods(refined=TRUE, refinement_scheme='reduced_dim')`:
    each sampled cell is replaced by the member of its own neighbourhood nearest to the
    neighbourhood's median position, which pulls the indices off the edges of the graph and
    leaves fewer, more representative neighbourhoods. Duplicates are dropped.
    """
    n = X.shape[0]
    nn = NearestNeighbors(n_neighbors=k + 1).fit(X)
    dist, idx = nn.kneighbors(X)                         # column 0 is the cell itself
    rng = np.random.default_rng(seed)
    sampled = rng.choice(n, size=max(1, int(round(prop * n))), replace=False)

    refined = np.empty_like(sampled)
    for j, c in enumerate(sampled):
        members = idx[c]
        med = np.median(X[members], axis=0)
        refined[j] = members[np.argmin(((X[members] - med) ** 2).sum(axis=1))]
    index = np.unique(refined)

    rows = np.repeat(np.arange(len(index)), k + 1)
    cols = idx[index].ravel()
    member = sparse.csr_matrix((np.ones_like(cols, dtype=np.float32), (rows, cols)),
                               shape=(len(index), n))
    kth = dist[index, k]
    print(f"[nhoods] {len(sampled)} sampled -> {len(index)} after refinement, "
          f"{k + 1} cells each")
    return index, member, kth


def spatial_fdr(p: np.ndarray, kth: np.ndarray) -> np.ndarray:
    """miloR's `graphSpatialFDR` with weighting 'k-distance': weighted BH, w = 1 / kth.

    A dense region has many overlapping neighbourhoods with small k-distances; down-weighting
    them stops one region from passing BH by sheer repetition.
    """
    out = np.full_like(p, np.nan, dtype=float)
    ok = ~np.isnan(p)
    if not ok.any():
        return out
    pv, w = p[ok], 1.0 / np.maximum(kth[ok], 1e-12)
    o = np.argsort(pv)
    pv, w = pv[o], w[o]
    adj = w.sum() * pv / np.cumsum(w)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    res = np.empty_like(adj)
    res[o] = np.minimum(adj, 1.0)
    out[np.where(ok)[0]] = res
    return out


# --------------------------------------------------------------------------- #
# The GLM, in R
# --------------------------------------------------------------------------- #


def run_edger(counts: pd.DataFrame, samples: pd.DataFrame, workdir: Path, rscript: str,
              reuse: bool) -> dict[str, pd.DataFrame]:
    workdir.mkdir(parents=True, exist_ok=True)
    done = sorted(workdir.glob("edger_*.csv"))
    saved = workdir / "nhood_counts.csv"
    if reuse and done and saved.exists():
        on_disk = pd.read_csv(saved, index_col=0)
        if not (on_disk.shape == counts.shape and (on_disk.columns == counts.columns).all()
                and (on_disk.to_numpy() == counts.to_numpy()).all()):
            print("[reuse] the counts on disk differ from these, refitting")
            reuse = False
    if not (reuse and done):
        for f in done:
            f.unlink()
        counts.to_csv(workdir / "nhood_counts.csv")
        samples.reset_index().to_csv(workdir / "samples.csv", index=False)
        print(f"[edgeR] {rscript} {R_SCRIPT.name} {workdir}")
        proc = subprocess.run([rscript, str(R_SCRIPT), str(workdir)], text=True,
                              capture_output=True)
        sys.stdout.write(proc.stderr)
        if proc.returncode != 0:
            sys.exit(f"milo_edger.R failed (exit {proc.returncode})")
    else:
        print(f"[reuse] {len(done)} edgeR outputs in {workdir}")
    out = {}
    for f in sorted(workdir.glob("edger_*.csv")):
        out[f.stem[len("edger_"):]] = pd.read_csv(f, index_col=0)
    return out


# --------------------------------------------------------------------------- #
# Figures
# --------------------------------------------------------------------------- #


def fig_umap(umap, index, nh: pd.DataFrame, contrasts, name, title, dpi) -> None:
    """Neighbourhoods at their index cell on the 05_3 UMAP; coloured only where DA."""
    ncol = len(contrasts)
    fig, axes = plt.subplots(1, ncol, figsize=(4.6 * ncol, 4.6), squeeze=False)
    lims = [float(np.nanmax(np.abs(nh[f"logFC_{c}"]))) for c in contrasts]
    lim = max([v for v in lims if np.isfinite(v)] or [1.0])
    for ax, c in zip(axes[0], contrasts):
        ax.scatter(umap[:, 0], umap[:, 1], s=0.3, c="#e6e5e1", lw=0, rasterized=True)
        sig = (nh[f"SpatialFDR_{c}"] < SPATIAL_FDR_ALPHA).to_numpy()
        xy = umap[index]
        ax.scatter(xy[~sig, 0], xy[~sig, 1], s=6, c="#c9c8c3", lw=0, rasterized=True)
        v = nh[f"logFC_{c}"].to_numpy()[sig]
        o = np.argsort(np.abs(v))
        sc = ax.scatter(xy[sig][o, 0], xy[sig][o, 1], s=14,
                        c=v[o], cmap=T.diverging_cmap(),
                        vmin=-lim, vmax=lim, lw=0.3, edgecolors="white", rasterized=True)
        ax.set_title(f"{c.replace('_', ' ')}   ({int(sig.sum())} DA of {len(nh)})",
                     fontsize=9, loc="left")
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_color("0.85")
    cb = fig.colorbar(sc, ax=axes[0].tolist(), fraction=0.02, pad=0.01)
    cb.set_label(f"logFC (SpatialFDR < {SPATIAL_FDR_ALPHA})", fontsize=7)
    cb.ax.tick_params(labelsize=6.5)
    cb.outline.set_visible(False)
    fig.suptitle(title, x=0.01, ha="left", fontsize=11, y=1.02)
    T.footnote(fig)
    C.savefig(name, T.STEP, T.CONSENSUS, fig=fig, dpi=dpi)
    plt.close(fig)


def fig_beeswarm(nh: pd.DataFrame, contrasts, name, title, dpi) -> None:
    """logFC per neighbourhood, one row per top programme: the miloR beeswarm, as a strip."""
    labels = (nh["top_programme"].value_counts().index.tolist())
    labels = [l for l in labels if l != "none"] + (["none"] if "none" in labels else [])
    ncol = len(contrasts)
    fig, axes = plt.subplots(1, ncol, figsize=(3.4 * ncol + 2.4, 0.32 * len(labels) + 1.6),
                             sharey=True, squeeze=False)
    rng = np.random.default_rng(0)
    for ax, c in zip(axes[0], contrasts):
        for i, lab in enumerate(labels):
            d = nh[nh["top_programme"] == lab]
            y = i + rng.uniform(-0.28, 0.28, len(d))
            sig = (d[f"SpatialFDR_{c}"] < SPATIAL_FDR_ALPHA).to_numpy()
            ax.scatter(d[f"logFC_{c}"][~sig], y[~sig], s=5, c="#c9c8c3", lw=0)
            v = d[f"logFC_{c}"][sig].to_numpy()
            ax.scatter(v, y[sig], s=8, lw=0,
                       c=np.where(v > 0, T.DIVERGING[4], T.DIVERGING[0]))
        ax.axvline(0, color="0.6", lw=0.8)
        ax.set_title(c.replace("_", " "), fontsize=9, loc="left")
        ax.set_xlabel("logFC", fontsize=7.5)
        ax.tick_params(labelsize=7)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[0][0].set_yticks(range(len(labels)), [l.replace("_", " ") for l in labels],
                          fontsize=7)
    axes[0][0].invert_yaxis()
    handles = [Line2D([], [], ls="", marker="o", color=T.DIVERGING[4], label="enriched, DA"),
               Line2D([], [], ls="", marker="o", color=T.DIVERGING[0], label="depleted, DA"),
               Line2D([], [], ls="", marker="o", color="#c9c8c3",
                      label=f"SpatialFDR >= {SPATIAL_FDR_ALPHA}")]
    fig.legend(handles=handles, loc="upper right", ncol=3, fontsize=7, frameon=False)
    fig.suptitle(title, x=0.01, ha="left", fontsize=11, y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    T.footnote(fig)
    C.savefig(name, T.STEP, T.CONSENSUS, fig=fig, dpi=dpi)
    plt.close(fig)


def fig_association(assoc: pd.DataFrame, name, title, dpi) -> None:
    M = assoc.pivot(index="programme", columns="test", values="rho")
    M = M.reindex(columns=list(dict.fromkeys(assoc["test"])))
    fig, ax = plt.subplots(figsize=(0.55 * M.shape[1] + 3.6, 0.3 * M.shape[0] + 1.8))
    im = ax.imshow(M.to_numpy(float), cmap=T.diverging_cmap(), vmin=-1, vmax=1,
                   aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M.iat[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=5.5,
                        color="white" if abs(v) > 0.6 else "0.15")
    ax.set_xticks(range(M.shape[1]), [c.replace("_", " ") for c in M.columns], fontsize=6.5,
                  rotation=60, ha="left")
    ax.xaxis.tick_top()
    ax.set_yticks(range(M.shape[0]), [p.replace("_", " ") for p in M.index], fontsize=7)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("Spearman, nhood logFC vs nhood programme mean", fontsize=7)
    cb.ax.tick_params(labelsize=6.5)
    cb.outline.set_visible(False)
    ax.set_title(title, fontsize=9, loc="left", pad=70)
    fig.tight_layout()
    T.footnote(fig)
    C.savefig(name, T.STEP, T.CONSENSUS, fig=fig, dpi=dpi)
    plt.close(fig)


# --------------------------------------------------------------------------- #


def main():
    args = parse_args()
    C.banner("05_11 - Milo: differential abundance of malignant states across the timepoints")
    print(f"run {C.RUN_ID}; k={args.k}, prop={args.prop}, min cells/biopsy={args.min_cells}")

    adata, L = T.load_embedding()
    prog = T.read_programmes()

    samples = T.sample_table(adata.obs, args.min_cells)
    T.report_samples(samples, args.min_cells)
    # A patient with one usable biopsy has no within-patient contrast: its cohort coefficient
    # would absorb it entirely. Out of the graph and the counts, and said so.
    kept = samples[samples["kept"]]
    n_tp = kept.groupby("cohort").size()
    single = n_tp.index[n_tp < 2]
    if len(single):
        print(f"patients with a single usable biopsy, left out: {', '.join(single)}")
    kept = kept[~kept["cohort"].isin(single)]

    sid = T.sample_ids(adata.obs)
    cell_mask = sid.isin(kept.index).to_numpy()
    dims = C.analysis_dimensions(adata)
    X = L.loc[cell_mask, dims].to_numpy(np.float32)
    print(f"{cell_mask.sum()} cells in {len(kept)} biopsies of {kept['cohort'].nunique()} "
          f"patients, graph on {len(dims)} dimensions")

    index, member, kth = build_nhoods(X, args.k, args.prop, args.seed)
    rule = 5 * len(kept)
    print(f"[nhoods] size {args.k + 1} vs Milo's rule of thumb ~{rule} (5 x {len(kept)} "
          f"biopsies){'  <- below it: expect sparse counts' if args.k + 1 < rule else ''}")

    cell_sample = pd.Categorical(sid[cell_mask].values, categories=kept.index)
    onehot = sparse.csr_matrix((np.ones(cell_mask.sum(), dtype=np.float32),
                                (np.arange(cell_mask.sum()), cell_sample.codes)),
                               shape=(cell_mask.sum(), len(kept)))
    counts = pd.DataFrame((member @ onehot).toarray().astype(int),
                          index=[f"nh{i}" for i in range(len(index))], columns=kept.index)

    # Same readouts as the other script, so a neighbourhood and a biopsy are read in one unit.
    directions = sorted({d for ds in prog["directions"] for d in ds}, key=C.dim_sort_key)
    P = T.programme_composites(T.oriented_directions(L, directions), prog)
    Pm = P.loc[cell_mask].to_numpy(float)
    nh_means = pd.DataFrame((member @ Pm) / (args.k + 1), index=counts.index,
                            columns=P.columns)

    workdir = (C.TUM_DIR / f"milo_{C.RUN_ID}{C.OUT_TAG}"
                         f"_k{args.k}_p{args.prop}_min{args.min_cells}_s{args.seed}")
    fits = run_edger(counts, kept, workdir, find_rscript(args.rscript), args.reuse)
    if not fits:
        sys.exit("edgeR wrote no contrast")

    nh = pd.DataFrame(index=counts.index)
    nh["index_cell"] = L.index[cell_mask][index]
    nh["size"] = args.k + 1
    nh["kth_distance"] = kth
    tp = kept["treatment"]
    for t in T.TREATMENTS:
        cols = tp.index[tp == t]
        nh[f"frac_{t}"] = counts[cols].sum(axis=1) / counts.sum(axis=1)
    top = nh_means.idxmax(axis=1)
    nh["top_programme"] = np.where(nh_means.max(axis=1) >= TOP_PROGRAMME_MIN, top, "none")
    for name, tab in fits.items():
        tab = tab.reindex(counts.index)
        nh[f"logFC_{name}"] = tab["logFC"]
        nh[f"PValue_{name}"] = tab["PValue"]
        nh[f"SpatialFDR_{name}"] = spatial_fdr(tab["PValue"].to_numpy(float), kth)
    nh = nh.join(nh_means.add_prefix("mean_"))
    C.write_table(nh, f"milo_nhoods_k{args.k}", T.CONSENSUS)

    tests = list(fits)
    summary = pd.DataFrame([{
        "test": t, "model": "all" if t.startswith("all_") else "by_response (exploratory)",
        "n_nhoods": int(nh[f"PValue_{t}"].notna().sum()),
        "n_DA": int((nh[f"SpatialFDR_{t}"] < SPATIAL_FDR_ALPHA).sum()),
        "n_enriched": int(((nh[f"SpatialFDR_{t}"] < SPATIAL_FDR_ALPHA)
                           & (nh[f"logFC_{t}"] > 0)).sum()),
        "n_depleted": int(((nh[f"SpatialFDR_{t}"] < SPATIAL_FDR_ALPHA)
                           & (nh[f"logFC_{t}"] < 0)).sum()),
    } for t in tests])
    C.write_table(summary, f"milo_summary_k{args.k}", T.CONSENSUS, index=False)
    C.banner(f"DA neighbourhoods at SpatialFDR < {SPATIAL_FDR_ALPHA}")
    print(summary.to_string(index=False))

    # Which programmes the moving neighbourhoods carry. Rho only: see the module docstring.
    assoc = pd.DataFrame([{
        "test": t, "programme": p,
        "rho": float(spearmanr(nh[f"logFC_{t}"], nh[f"mean_{p}"], nan_policy="omit")
                     .statistic),
        "mean_in_enriched": float(nh.loc[(nh[f"SpatialFDR_{t}"] < SPATIAL_FDR_ALPHA)
                                         & (nh[f"logFC_{t}"] > 0), f"mean_{p}"].mean()),
        "mean_in_depleted": float(nh.loc[(nh[f"SpatialFDR_{t}"] < SPATIAL_FDR_ALPHA)
                                         & (nh[f"logFC_{t}"] < 0), f"mean_{p}"].mean()),
    } for t in tests for p in P.columns])
    C.write_table(assoc, f"milo_programme_association_k{args.k}", T.CONSENSUS, index=False)

    umap = np.asarray(adata.obsm["X_umap"])[cell_mask] if "X_umap" in adata.obsm else None
    main_tests = [f"all_{c}" for c in T.CONTRAST_NAMES if f"all_{c}" in fits]
    if umap is not None and main_tests:
        fig_umap(umap, index, nh, main_tests, f"milo_nhood_umap_k{args.k}",
                 f"Milo neighbourhoods that change their share of the tumour, within patient (k={args.k})",
                 args.dpi)
    if main_tests:
        fig_beeswarm(nh, main_tests, f"milo_beeswarm_k{args.k}",
                     f"Neighbourhood logFC, grouped by the programme the neighbourhood is "
                     f"highest on (k={args.k})", args.dpi)
    fig_association(assoc, f"milo_programme_association_k{args.k}",
                    f"Which programmes the moving neighbourhoods carry, k={args.k} "
                    f"(columns after the first three: by response, exploratory)", args.dpi)
    print(f"\n[done] {C.fig_dir(T.STEP, T.CONSENSUS)}")


if __name__ == "__main__":
    main()
