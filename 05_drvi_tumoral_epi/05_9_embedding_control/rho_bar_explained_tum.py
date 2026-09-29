#!/usr/bin/env python3
"""What `ROUTE_A_RHO_MIN = 0.30` means, in one figure. An explainer, not a step: it measures
nothing the pipeline uses and decides nothing.

Two rows, two questions.

    top     WHAT A rho LOOKS LIKE IN CELLS. Three real dimension-direction x programme pairs,
            one well above the bar, one on it, one below. The cells are cut into ten deciles
            of the ORIENTED coordinate (multiplied by the sign of the direction, so the named
            side is always to the right, as in dr_umap_tum.py) and each decile shows the
            within-cohort z of the programme - the two vectors Route A correlates. A rho is
            the slope-without-units of that staircase.

    bottom  WHY 0.30. `A_rho` is a direction's BEST signature, a maximum over the collection's
            K gene sets, so its null is the same maximum over K size- and expression-matched
            random sets: the null route_a_null_tum.py prints as percentiles, drawn here as a
            distribution beside the real A_rho of every non-vanished direction. The bar sits
            where the null is spent and the real directions are not.

Reads only what is already on disk: 05_3's embedding, 05_6's per-cell scores, 05_8's pruned
convergence tables and 05_9's cached matched random scores. Pruned (PRUNE_VANISHED=1), the
setting the interpretation uses.

Usage:
    export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
    N_LATENT=64 HVG_SET=nomt python3 rho_bar_explained_tum.py
"""
from __future__ import annotations

import os
import sys
from types import SimpleNamespace

os.environ["PRUNE_VANISHED"] = "1"

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "utils"))
import signature_common as C  # noqa: E402
import sig_collections as SC  # noqa: E402
import cell_set as CS  # noqa: E402

STEP = "05_9_embedding_control"
CONSENSUS = SimpleNamespace(name="consensus", title="all three collections, read together")
BAR = 0.30
N_PSEUDO = 400          # as route_a_null_tum.py
SEED = 0

COLLS = {"scie": SC.SCIE, "gavish_tnbc": SC.GAVISH_TNBC, "emt": SC.EMT}

# (collection, direction, signature). The low one is picked from the table, not hand-chosen.
EXAMPLES = [("gavish_tnbc", "DR 3+", "MP1_CELL_CYCLE_G2_M"),
            ("scie", "DR 19-", "LIM_STEM")]
LOW_TARGET = 0.15

PURPLE, GREY, INK, MUTED = "#6A1FC2", "#BDBDBD", "#222222", "#6B6B6B"
NULL_C, REAL_C = "#9E9E9E", "#2F6FB0"


def ranked_z(X: np.ndarray) -> np.ndarray:
    R = np.apply_along_axis(rankdata, 0, X)
    return (R - R.mean(0)) / R.std(0)


def main():
    embed = ad.read_h5ad(C.get_embedding(C.DEFAULT_EMBEDDING).embed_h5ad)
    titles = embed.var["title"].to_numpy()
    live = ~embed.var["vanished"].astype(bool).to_numpy()
    X = np.asarray(embed.X)
    cells = embed.obs_names.astype(str)

    conv = {k: pd.read_csv(C.table_path("convergence", c), comment="#") for k, c in COLLS.items()}

    # The low example: the direction whose best signature is closest to LOW_TARGET, on gavish.
    g = conv["gavish_tnbc"]
    g = g.loc[~g["dimension_vanished"].astype(bool)]
    low = g.iloc[(g["A_rho"] - LOW_TARGET).abs().argsort().iloc[0]]
    examples = EXAMPLES + [("gavish_tnbc", low["dim_direction"], low["A_best_signature"])]

    fig = plt.figure(figsize=(16, 10.5))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1], hspace=0.55, wspace=0.28)

    # ------------------------------------------------------------ top row
    score_cache = {}
    for i, (cname, dd, sig) in enumerate(examples):
        if cname not in score_cache:
            s = pd.read_csv(C.scores_csv(COLLS[cname]), index_col=0)
            s.index = s.index.astype(str)
            score_cache[cname] = s.reindex(cells)
        z = score_cache[cname][f"z_{sig}"].to_numpy()
        dim, sign = dd[:-1], (1 if dd[-1] == "+" else -1)
        coord = sign * X[:, list(titles).index(dim)]
        rho = spearmanr(coord, z)[0]

        dec = pd.qcut(rankdata(coord), 10, labels=False)
        groups = [z[dec == d] for d in range(10)]

        ax = fig.add_subplot(gs[0, i])
        bp = ax.boxplot(groups, positions=range(1, 11), widths=0.6, showfliers=False,
                        patch_artist=True, medianprops={"color": INK, "lw": 1.4},
                        whiskerprops={"color": MUTED}, capprops={"color": MUTED})
        above = rho >= BAR
        for b in bp["boxes"]:
            b.set(facecolor=PURPLE if above else GREY, alpha=0.55 if above else 0.8,
                  edgecolor="white", lw=1)
        ax.axhline(0, color=MUTED, lw=0.8, ls=":")
        top = np.mean(groups[9] > 0) * 100
        bottom = np.mean(groups[0] > 0) * 100
        verdict = "above the bar" if above else "below the bar"
        ax.set_title(f"{sig.replace('_', ' ')}  on  {dd}\n"
                     f"rho = {rho:+.2f}   ({verdict})", fontsize=11, color=INK)
        ax.text(0.02, 0.97,
                f"cells above their cohort mean\n"
                f"lowest decile: {bottom:.0f}%   highest: {top:.0f}%",
                transform=ax.transAxes, va="top", fontsize=8.5, color=INK,
                bbox={"fc": "white", "ec": "none", "alpha": 0.85})
        ax.set_xticks(range(1, 11))
        ax.set_xticklabels([str(d) for d in range(1, 11)], fontsize=8)
        ax.set_xlabel(f"decile of the coordinate on {dd}\n(1 = furthest from the named side, "
                      f"10 = furthest into it)", fontsize=9, color=MUTED)
        if i == 0:
            ax.set_ylabel("within-cohort z of the programme", fontsize=9)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.grid(axis="y", alpha=0.3)
        ax.set_axisbelow(True)

    # ---------------------------------------------------------- bottom row
    Lr = ranked_z(X[:, live])
    rng = np.random.default_rng(SEED)
    for i, cname in enumerate(("scie", "gavish_tnbc", "emt")):
        coll = COLLS[cname]
        rnd = pd.read_csv(CS.tum_dir() / f"random_signature_scores_{coll.name}_{C.RUN_ID}.csv",
                          index_col=0)
        rnd.index = rnd.index.astype(str)
        rnd = rnd.reindex(cells)
        R = pd.DataFrame(Lr.T @ ranked_z(rnd.to_numpy()) / Lr.shape[0], columns=rnd.columns)
        base = sorted({c.split("__rnd")[0] for c in rnd.columns})
        draws = sorted({int(c.split("__rnd")[1]) for c in rnd.columns})
        vals = []
        for _ in range(N_PSEUDO):
            sub = R[[f"{b}__rnd{rng.choice(draws)}" for b in base]].to_numpy()
            vals += [sub.max(axis=1), (-sub).max(axis=1)]
        null = np.concatenate(vals)

        c = conv[cname]
        real = c.loc[~c["dimension_vanished"].astype(bool), "A_rho"].to_numpy()
        n = len(real)
        p_null = (null >= BAR).mean()
        n_real = int((real >= BAR).sum())

        ax = fig.add_subplot(gs[1, i])
        bins = np.linspace(min(null.min(), real.min(), 0), max(real.max(), null.max()) + 0.02, 45)
        ax.hist(null, bins=bins, density=True, color=NULL_C, alpha=0.6,
                label="random gene sets (null)")
        ax.hist(real, bins=bins, density=True, histtype="step", lw=2, color=REAL_C,
                label="real signatures")
        ax.set_ylim(0, ax.get_ylim()[1] * 1.12)       # headroom for the bar label
        ax.axvline(BAR, color=INK, lw=1.3, ls="--")
        ax.text(BAR + 0.008, 0.97, f"bar {BAR:.2f}", transform=ax.get_xaxis_transform(),
                fontsize=9, color=INK, va="top")
        ax.set_title(f"{cname}  ({len(base)} gene sets)\n"
                     f"past {BAR:.2f}:  by chance {p_null * 100:.1f}% (~{p_null * n:.0f} of {n})"
                     f"  |  real {n_real} of {n}", fontsize=10.5, color=INK)
        ax.set_xlabel("A_rho: best Spearman of a direction over the collection", fontsize=9)
        if i == 0:
            ax.set_ylabel("density", fontsize=9)
            handles, labels = ax.get_legend_handles_labels()
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)

    fig.text(0.07, 0.955, "A.  What a rho looks like in the cells", fontsize=13,
             weight="bold", color=INK)
    fig.text(0.07, 0.475, "B.  Why the bar is at 0.30: how often a random gene set of the same "
             "size gets there", fontsize=13, weight="bold", color=INK)
    fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False, fontsize=10,
               bbox_to_anchor=(0.5, 0.035))
    fig.suptitle(f"The Route A bar, rho >= {BAR:.2f}  -  {C.RUN_ID}, malignant cells, "
                 f"vanished dimensions pruned", fontsize=14, y=1.0)
    C.savefig("rho_bar_explained", STEP, CONSENSUS, fig=fig)


if __name__ == "__main__":
    main()
