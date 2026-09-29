#!/usr/bin/env python3
"""05_6, side figure: each DRVI direction against the SCIE readouts, as an AUROC.

Route A's heatmap correlates the whole coordinate of a dimension with each readout, so one
number stands for both sides of it. DRVI reads the two sides as possibly independent
programmes (Moinfar & Theis 2024), and this figure keeps them apart.

For every dimension and every side:

    active   = cells beyond the side's threshold (DR n+: value > ACTIVE_CUT, DR n-: < -ACTIVE_CUT)
    inactive = cells with |value| < INACTIVE_CUT, the same reference for both sides
    AUROC(active vs inactive) on each readout's within-cohort z-score

oriented so that positive always points at the SCIE target:

    stemness readouts:  AUROC - 0.5   (> 0: the active cells are MORE stem-like)
    immune readouts:    0.5 - AUROC   (> 0: the active cells are LESS immunogenic, i.e. evasive)

So a SCIE side is a row that is purple all the way across.

Inputs are exactly what 05_6 uses: the per-cell z-scores from its cache
(`C.scores_csv(coll)`, written by `cell_first_tum.py`, not recomputed here) and the DRVI
embedding at `C.EMBED_H5AD`, restricted to `C.analysis_dimensions`. Read-only: nothing in
05_6 / 05_7 / 05_8 reads what this writes. Run it with the chain's environment:

    PRUNE_VANISHED=1 OUT_TAG=_pruned_rho030 N_LATENT=64 python scie_side_heatmap_tum.py
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import anndata as ad  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402

UTILS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils")
sys.path.insert(0, UTILS_DIR)
import signature_common as C  # noqa: E402
import sig_collections as SC  # noqa: E402

# A cell is active on DR n+ if its value is > ACTIVE_CUT, on DR n- if it is < -ACTIVE_CUT.
# 1 is DRVI's vanished cutoff (max |value| < 1), reused per cell.
ACTIVE_CUT = 1.0
# The reference group: cells near zero on the dimension, |value| < INACTIVE_CUT.
INACTIVE_CUT = 0.5
# A side is kept only if at least this share of cells (in %) is active on it.
MIN_SIDE_PCT = 1.0
# Colour scale of the oriented AUROC.
VMAX = 0.3
N_TOP = 10

STEP = "05_6_cell_first"
COL_DROPPED = "0.55"


def side_table(L: pd.DataFrame, z: pd.DataFrame, readouts: list[str],
               axis_of: dict[str, str]) -> pd.DataFrame:
    rows = []
    n = len(L)
    for d in sorted(L.columns, key=lambda t: int(t.split()[1])):
        x = L[d].values
        inactive = np.abs(x) < INACTIVE_CUT
        for sign, active in (("+", x > ACTIVE_CUT), ("-", x < -ACTIVE_CUT)):
            row = {"side": f"{d}{sign}", "dimension": d, "direction": sign,
                   "n_active": int(active.sum()), "pct_active": 100 * active.sum() / n,
                   "n_inactive": int(inactive.sum())}
            for r in readouts:
                if active.sum() == 0 or inactive.sum() == 0:
                    row[r] = np.nan
                    continue
                keep = active | inactive
                auc = roc_auc_score(active[keep], z[f"z_{r}"].values[keep])
                row[r] = auc - 0.5 if axis_of[r] == "stemness" else 0.5 - auc
            rows.append(row)
    tab = pd.DataFrame(rows).set_index("side")
    tab["kept"] = tab["pct_active"] >= MIN_SIDE_PCT
    stem = [r for r in readouts if axis_of[r] == "stemness"]
    imm = [r for r in readouts if axis_of[r] == "immune"]
    tab["max_stemness"] = tab[stem].max(axis=1)
    tab["max_evasion"] = tab[imm].max(axis=1)
    tab["scie_score"] = tab[["max_stemness", "max_evasion"]].min(axis=1)
    return tab


def plot(tab: pd.DataFrame, cols: list[str], edges: list[int], coll, full: bool) -> None:
    t = tab if full else tab[tab["kept"]]
    m = t[cols].astype(float).rename(
        columns={c: f"{c} (evasion)" for c in cols if coll.axis_of[c] == "immune"})
    fig, ax = plt.subplots(figsize=(C.fig_span(len(cols), 1.0, 4.0),
                                    C.fig_span(len(m), 0.16, 3.0)))
    im = ax.imshow(m.values, cmap="PuOr", vmin=-VMAX, vmax=VMAX, aspect="auto",
                   interpolation="nearest")
    ax.set_facecolor(COL_DROPPED)                  # NaN (no active cell) shows as grey
    if full:
        # Dropped sides stay readable underneath, but greyed out.
        for i, kept in enumerate(t["kept"].values):
            if not kept:
                ax.axhspan(i - 0.5, i + 0.5, color=COL_DROPPED, alpha=0.75, lw=0)
    for pos in edges:
        ax.axvline(pos - 0.5, color="k", lw=1.5)
    # A thin line between dimensions, so the two sides of one read as a pair.
    for i in range(1, len(m)):
        if t["dimension"].iloc[i] != t["dimension"].iloc[i - 1]:
            ax.axhline(i - 0.5, color="white", lw=0.6)
    ax.set_xticks(range(len(m.columns)))
    ax.set_xticklabels(m.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(m)))
    ax.set_yticklabels([f"{s}  ({p:.1f}%)" for s, p in zip(m.index, t["pct_active"])],
                       fontsize=5)
    if full:
        for lab, kept in zip(ax.get_yticklabels(), t["kept"].values):
            if not kept:
                lab.set_color(COL_DROPPED)
    cb = fig.colorbar(im, ax=ax, shrink=0.3)
    cb.set_label("oriented AUROC (active vs inactive cells, within-cohort z)\n"
                 "purple = towards SCIE (more stem / less immune), orange = away", fontsize=8)
    n_drop = int((~tab["kept"]).sum())
    which = (f"all {len(tab)} sides; the {n_drop} with < {MIN_SIDE_PCT:g}% active cells in grey"
             if full else
             f"{len(m)} of {len(tab)} sides kept (>= {MIN_SIDE_PCT:g}% active cells)")
    ax.set_title(f"Route A, {coll.title}: DRVI directions x signatures, {C.RUN_ID}\n"
                 f"active: value beyond +/-{ACTIVE_CUT:g}; inactive: |value| < {INACTIVE_CUT:g}; "
                 f"stemness AUROC-0.5, immune 0.5-AUROC\n{which}; "
                 "a SCIE side is a row that is purple all the way across", fontsize=9)
    C.savefig("dim_side_auroc_heatmap" + ("_full" if full else ""), STEP, coll, fig)
    plt.close(fig)


def main() -> None:
    coll = SC.get("scie")
    C.banner(f"SCIE side AUROC | run {C.RUN_ID} | {C.EMBED_H5AD}")

    embed = ad.read_h5ad(C.EMBED_H5AD)
    dims = C.analysis_dimensions(embed)
    L = pd.DataFrame(np.asarray(embed.X), index=embed.obs_names,
                     columns=embed.var["title"].values)[dims]
    print(f"{embed.n_vars} dimensions, {len(dims)} used (PRUNE_VANISHED={C.PRUNE_VANISHED})")

    scores = C.scores_csv(coll)
    print(f"[read] {scores}")
    z = pd.read_csv(scores, index_col=0)
    assert (z.index == L.index).all(), "score cache and embedding are not in the same cell order"
    readouts = coll.order([c[2:] for c in z.columns if c.startswith("z_")], for_figure=True)
    edges = coll.block_edges(readouts)

    tab = side_table(L, z, readouts, coll.axis_of)
    C.write_table(tab, "dim_side_auroc", coll)

    plot(tab, readouts, edges, coll, full=False)
    plot(tab, readouts, edges, coll, full=True)

    top = tab[tab["kept"]].sort_values("scie_score", ascending=False).head(N_TOP)
    print(f"\ntop {N_TOP} kept sides by min(max stemness, max evasion):")
    print(top[["pct_active", "max_stemness", "max_evasion", "scie_score"]]
          .to_string(float_format="%+.3f"))


if __name__ == "__main__":
    main()
