#!/usr/bin/env python3
"""05_6, side analysis: is each DRVI dimension active on one side or on both?

In DRVI (Moinfar & Theis 2024) a dimension is read in its two directions, DR n+ and DR n-,
which can be independent programmes. This script asks, for every dimension of the run, how
many cells actually reach each side. Read-only: it reads the 05_3 embedding and writes one
table and one figure, and nothing in 05_6 / 05_7 / 05_8 reads either of them.

The activity threshold is |value| > 1, the same number the paper uses to call a dimension
vanished (max |value| < 1), reused here per cell and per side.

Run it with the same environment as the chain, so the run id and the output tag resolve to
the run being read:

    N_LATENT=64 OUT_TAG=_pruned_rho030 python dr_sidedness_tum.py

Outputs, both keyed on the run and on no collection (the per-dimension distribution does not
depend on which gene lists are being read):

    tables<OUT_TAG>/consensus/<run_id>/dr_sidedness_<run_id>.csv
    figures<OUT_TAG>/05_6_cell_first/dr_sidedness_violins_<run_id>.png

The figure is a grid of histograms with a log y axis rather than violins: the mass near zero
flattens a violin to a line and the tails, which are the whole question, disappear.
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import anndata as ad  # noqa: E402

UTILS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils")
sys.path.insert(0, UTILS_DIR)
import signature_common as C  # noqa: E402

# --------------------------------------------------------------------------- #
# Thresholds
# --------------------------------------------------------------------------- #

# A cell is active on the + side of a dimension if its value is > ACTIVE_CUT, on the - side
# if it is < -ACTIVE_CUT. 1 is DRVI's vanished cutoff (max |value| < 1), reused per cell.
ACTIVE_CUT = 1.0

# A side counts as present if at least this share of cells (in %) is active on it.
MIN_SIDE_PCT = 1.0

# A dimension is one-sided if its minor side is below MIN_SIDE_PCT AND the minor side holds
# less than this fraction of the major side's active cells (side_ratio = minor % / major %).
MAX_SIDE_RATIO = 0.1

PERCENTILES = (1, 5, 50, 95, 99)

STEP = "05_6_cell_first"
COL_HIST = "#2a78d6"
COL_CUT = "0.25"


def classify(pct_pos: float, pct_neg: float, ratio: float) -> str:
    major_pct = max(pct_pos, pct_neg)
    minor_pct = min(pct_pos, pct_neg)
    if major_pct < MIN_SIDE_PCT:
        return "inactive"
    if minor_pct >= MIN_SIDE_PCT:
        return "two-sided"
    if ratio < MAX_SIDE_RATIO:
        return "one-sided +" if pct_pos > pct_neg else "one-sided -"
    # Minor side below 1% but not small relative to the major one: only possible when the
    # major side is itself below 10%. Left unnamed rather than forced into a class.
    return "unclassified"


def sidedness_table(X: np.ndarray, var: pd.DataFrame) -> pd.DataFrame:
    n = X.shape[0]
    rows = []
    for j, (_, v) in enumerate(var.iterrows()):
        x = X[:, j]
        n_pos = int((x > ACTIVE_CUT).sum())
        n_neg = int((x < -ACTIVE_CUT).sum())
        pct_pos, pct_neg = 100 * n_pos / n, 100 * n_neg / n
        major = max(pct_pos, pct_neg)
        ratio = min(pct_pos, pct_neg) / major if major > 0 else np.nan
        row = {
            "dimension": v["title"],
            "vanished": bool(v["vanished"]),
            "vanished_positive_direction": bool(v["vanished_positive_direction"]),
            "vanished_negative_direction": bool(v["vanished_negative_direction"]),
            "n_pos": n_pos, "pct_pos": pct_pos,
            "n_neg": n_neg, "pct_neg": pct_neg,
            "min": float(x.min()), "max": float(x.max()),
        }
        for p, q in zip(PERCENTILES, np.percentile(x, PERCENTILES)):
            row[f"p{p:02d}"] = float(q)
        row["side_ratio"] = ratio
        row["sidedness"] = classify(pct_pos, pct_neg, ratio)
        rows.append(row)
    df = pd.DataFrame(rows)
    df["_n"] = df["dimension"].str.extract(r"(\d+)").astype(int)
    return df.sort_values("_n").drop(columns="_n").reset_index(drop=True)


def plot_grid(X: np.ndarray, var: pd.DataFrame, tab: pd.DataFrame, path) -> None:
    keep = tab.loc[~tab["vanished"], "dimension"].tolist()
    col_of = {t: i for i, t in enumerate(var["title"])}
    n_cols = 8
    n_rows = int(np.ceil(len(keep) / n_cols))
    lim = float(np.abs(X[:, [col_of[d] for d in keep]]).max())
    bins = np.linspace(-lim, lim, 81)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(2.3 * n_cols, 1.9 * n_rows),
                             sharex=True, squeeze=False)
    rec = tab.set_index("dimension")
    for ax, d in zip(axes.flat, keep):
        ax.hist(X[:, col_of[d]], bins=bins, color=COL_HIST, log=True)
        for c in (ACTIVE_CUT, -ACTIVE_CUT):
            ax.axvline(c, color=COL_CUT, ls="--", lw=0.8)
        r = rec.loc[d]
        ax.set_title(f"{d}  +: {r.pct_pos:.1f}%  -: {r.pct_neg:.1f}%", fontsize=7.5)
        ax.tick_params(labelsize=6)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    for ax in list(axes.flat)[len(keep):]:
        ax.set_visible(False)
    fig.suptitle(f"DRVI {C.RUN_ID}: per-dimension value distribution, "
                 f"{X.shape[0]:,} cells, {len(keep)} non-vanished dimensions "
                 f"(dashed: +/-{ACTIVE_CUT:g}; y log)", fontsize=10)
    fig.supxlabel("latent value", fontsize=8)
    fig.supylabel("cells (log)", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {path}")


def main() -> None:
    C.banner(f"DR sidedness | run {C.RUN_ID} | {C.EMBED_H5AD}")
    embed = ad.read_h5ad(C.EMBED_H5AD)
    X = np.asarray(embed.X, dtype=np.float64)   # cell order as stored
    var = embed.var

    tab = sidedness_table(X, var)

    table_dir = C.TABLE_DIR / "consensus" / C.RUN_ID
    table_dir.mkdir(parents=True, exist_ok=True)
    table_path = table_dir / f"dr_sidedness_{C.RUN_ID}.csv"
    with open(table_path, "w", encoding="utf-8") as fh:
        fh.write(f"# dr_sidedness | run {C.RUN_ID} | {X.shape[0]} cells | active: |value| > "
                 f"{ACTIVE_CUT:g} | side present: >= {MIN_SIDE_PCT:g}% cells | one-sided: "
                 f"minor side < {MIN_SIDE_PCT:g}% and side_ratio < {MAX_SIDE_RATIO:g}\n")
        tab.to_csv(fh, index=False)
    print(f"[table] {table_path}  ({tab.shape[0]} x {tab.shape[1]})")

    fig_dir = C.OUT_ROOT / f"figures{C.OUT_TAG}" / STEP
    fig_dir.mkdir(parents=True, exist_ok=True)
    plot_grid(X, var, tab, fig_dir / f"dr_sidedness_violins_{C.RUN_ID}.png")

    print("\nclass counts (all dimensions):")
    print(tab["sidedness"].value_counts().to_string())
    print("\nclass counts (non-vanished):")
    print(tab.loc[~tab["vanished"], "sidedness"].value_counts().to_string())


if __name__ == "__main__":
    main()
