#!/usr/bin/env python3
"""04_8: how much of "stemness" on this dataset is the cell cycle.

Every earlier step reports the cycle as a RISK - `confounders_*.csv` correlates each score
against `S_score` and `G2M_score`, 04_5 recomputes the target quadrant inside G1 alone, 04_7
carries a `cell_cycle` flag onto any row whose claimed signature is coupled to it. What none of
them does is put the four pieces of that story in one place and show that they are one story.
This step does, in four panels, following the confound from the gene lists to the latent space:

  A. THE MECHANISM. What fraction of each score's variance is carried by genes that are literally
     in the Regev cell-cycle list. This is a fact about the LISTS, not about the cells.

  B. THE CONSEQUENCE. R^2 of (S_score, G2M_score) on each within-stratum standardised score,
     against panel A on the x axis. If A causes B the points lie on a line through the origin,
     and the lists with no cycle genes sit at the origin.

  C. THE COMPOUNDING, and the panel that matters. AUROC of the cell cycle ALONE predicting
     membership of the target quadrant, as a function of how many stemness readouts have to
     agree. A majority vote is supposed to average an idiosyncratic confounder away. It does the
     opposite here, because the cycle is the one thing the readouts have in common, so
     intersecting them selects for it.

  D. THE LATENT SPACE. Cycle loading of each of the 64 dimensions against its strongest stemness
     association. The two dimensions that clear the 0.30 cycle bar are both convergent rows of
     04_7, which is the sharpest form the result takes.

CIRCULARITY, STATED. `S_score` and `G2M_score` were computed in 04_1 from the same Regev list
panel A counts against, so A and B are not independent: a list containing Regev genes correlates
with a score built from Regev genes partly by construction. That is why A is framed as literal
content - a checkable property of the list - and why CytoTRACE2 is the control that carries the
argument: it has no gene list at all, so it cannot contain a cycle gene, and it still lands at
R^2 0.07.

Usage:
    export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
    python cycle_confound_epi.py                      # the scie collection, the default
    python cycle_confound_epi.py --collection emt     # the same four panels on the EMT lists
"""

from __future__ import annotations

import argparse
import os
import sys

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import rankdata, spearmanr

UTILS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils")
sys.path.insert(0, UTILS_DIR)
import signature_common as C  # noqa: E402
import sig_collections as SC  # noqa: E402


# The quadrant cutoffs of 04_5. Re-declared rather than imported because `cell_first_epi.py` is
# a script with its own argparse; the run below asserts it reproduces 04_5's published vote
# distribution exactly, which is what actually keeps the two in step.
HIGH_Q, LOW_Q = 0.75, 0.25

# 04_7's bar, so panel D's threshold line is the one the convergence table was flagged with.
CYCLE_FLAG = 0.30

CYCLE_GENES_FILE = "regev_lab_cell_cycle_genes.txt"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    SC.add_argument(p)
    return p.parse_args()


def auroc(score: np.ndarray, label: np.ndarray) -> float:
    """Rank AUROC of a continuous score against a boolean label. Ties get midranks."""
    r = rankdata(score)
    n1 = int(label.sum())
    n0 = len(label) - n1
    if n1 == 0 or n0 == 0:
        return np.nan
    return float((r[label].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def within_stratum_z(df: pd.DataFrame, keys: list) -> pd.DataFrame:
    """The standardisation of 04_5: z-score inside each (cohort, cell_type) stratum."""
    return (df.groupby(keys, observed=True)
              .transform(lambda s: (s - s.mean()) / s.std(ddof=0) if s.std(ddof=0) > 0 else s * 0.0)
              .fillna(0.0))


def main():
    args = parse_args()
    coll = SC.get(args.collection)

    C.banner(f"04_8 - the cell cycle behind the stemness axis: {coll.title}")
    print(f"question    how much of what these lists call 'stemness' is proliferation?")

    scores_path = C.scores_csv(coll)
    if not scores_path.exists():
        sys.exit(f"[STOP] {scores_path} not found - run 04_5 first; this step reads its scores.")

    cyc_path = C.DATA_DIR / CYCLE_GENES_FILE
    if not cyc_path.exists():
        sys.exit(f"[STOP] {cyc_path} not found - it is the list 04_1 scored the cycle with.")
    cycle_genes = set(C.read_signature_file(cyc_path))
    print(f"cycle list  {cyc_path.name}, {len(cycle_genes)} genes (the same list 04_1 used)")

    adata = ad.read_h5ad(C.FULL_H5AD, backed="r")
    obs = adata.obs
    for k in ("S_score", "G2M_score", "phase"):
        if k not in obs.columns:
            sys.exit(f"[STOP] obs['{k}'] missing from {C.FULL_H5AD.name} - 04_1 did not run.")

    scores = pd.read_csv(scores_path, index_col=0, comment="#").loc[obs.index]
    readouts = [r for r in coll.axis_of if f"score_{r}" in scores.columns]
    print(f"readouts    {len(readouts)} scored, on {len(set(coll.axis_of[r] for r in readouts))} axes")

    S = obs["S_score"].to_numpy()
    G = obs["G2M_score"].to_numpy()
    phase = obs["phase"].to_numpy()
    keys = [obs["cohort"].values, obs["cell_type"].values]

    # ------------------------------------------------------------------ A + B
    gene_tbl = C.read_table("signature_gene_contribution", coll).reset_index()
    rows = []
    z_all = within_stratum_z(scores[[f"score_{r}" for r in readouts]], keys)
    z_all.columns = [c.replace("score_", "") for c in z_all.columns]
    X = np.column_stack([np.ones(len(obs)), S, G])
    for r in readouts:
        sub = gene_tbl[gene_tbl["signature"] == r]
        if len(sub):
            hit = sub[sub["gene"].isin(cycle_genes)]
            n_cycle, share = len(hit), float(hit["share_of_score_variance"].sum())
        else:
            # CytoTRACE2 and any derived readout: no gene list, so no cycle genes BY
            # CONSTRUCTION. That is what makes it the non-circular control of panel B.
            n_cycle, share = 0, np.nan
        y = z_all[r].to_numpy()
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        r2 = 1.0 - ((y - X @ beta) ** 2).sum() / ((y - y.mean()) ** 2).sum()
        rows.append({"readout": r, "axis": coll.axis_of[r],
                     "n_genes": len(sub), "n_cycle_genes": n_cycle,
                     "cycle_share_of_variance": share, "r2_cycle_on_score": r2,
                     "has_gene_list": bool(len(sub))})
    per_readout = pd.DataFrame(rows).set_index("readout")

    print("\nA/B - cycle gene content, and the variance it explains")
    print(per_readout[["axis", "n_genes", "n_cycle_genes", "cycle_share_of_variance",
                       "r2_cycle_on_score"]].to_string(float_format="%.4f"))

    # -------------------------------------------------------------------- C
    # One target per PLANE, exactly as 04_5's `define_target` builds it. Going through the
    # collection's planes rather than assuming a shape is what makes this work on both
    # collections: SCIE crosses many stemness readouts against one fixed immune axis, while EMT
    # pairs a different (x, y) per list version and shares no axis at all.
    z = z_all

    def in_region(v: pd.Series, rule: str) -> pd.Series:
        if rule == "high":
            return v >= v.quantile(HIGH_Q)
        if rule == "low":
            return v <= v.quantile(LOW_Q)
        raise ValueError(f"panel C does not handle the {rule!r} region rule")

    planes = [pl for pl in coll.planes(coll, list(z.columns))
              if pl.x in z.columns and pl.y in z.columns]
    if len(planes) < 2:
        sys.exit(f"[STOP] {coll.name} defines {len(planes)} plane(s); panel C needs at least two "
                 "target definitions to ask what agreeing between them does.")
    Q = pd.DataFrame({pl.label: in_region(z[pl.x], pl.x_rule) & in_region(z[pl.y], pl.y_rule)
                      for pl in planes})
    votes = Q.sum(axis=1)

    published = C.read_table("quadrant_vote_distribution", coll)["n_cells"]
    mine = votes.value_counts().sort_index()
    if not all(int(mine.get(k, 0)) == int(published.get(k, 0)) for k in published.index):
        print("[warn] the recomputed vote distribution does not match 04_5's published one; "
              "panel C is drawn on the recomputation and the two have drifted apart.")
    else:
        print(f"\n[check] vote distribution reproduces 04_5 exactly "
              f"({int((votes >= max(2, int(np.ceil(len(planes)/2)))).sum()):,} cells at the consensus bar)")

    cyc_score = S + G
    curve = []
    for k in range(1, len(planes) + 1):
        lab = (votes >= k).to_numpy()
        curve.append({"min_votes": k, "n_cells": int(lab.sum()),
                      "auroc_cycle_predicts_membership": auroc(cyc_score, lab),
                      "pct_G1": float((phase[lab] == "G1").mean() * 100),
                      "pct_G2M": float((phase[lab] == "G2M").mean() * 100),
                      "median_depth": float(obs["n_genes_by_counts"].to_numpy()[lab].mean())})
    curve = pd.DataFrame(curve).set_index("min_votes")
    # A single readout's own quadrant, as the reference the consensus is supposed to improve on.
    single = {c: auroc(cyc_score, Q[c].to_numpy()) for c in Q.columns}
    print("\nC - does agreement average the cycle away, or concentrate it?")
    print(curve[["n_cells", "auroc_cycle_predicts_membership", "pct_G1", "pct_G2M"]]
          .to_string(float_format="%.3f"))
    print("  single-readout quadrants for reference: "
          + ", ".join(f"{r} {v:.3f}" for r, v in sorted(single.items(), key=lambda t: -t[1])))

    # -------------------------------------------------------------------- D
    embed = ad.read_h5ad(C.EPI_DIR / f"embed_{C.RUN_ID}.h5ad")
    E = np.asarray(embed.X)
    # `var['title']` is the label, NOT the column position: DRVI stores the dimensions in its own
    # order, so column 0 is 'DR 32' on this run. Using i+1 here silently mislabels every
    # dimension and quietly corrupts the join with 04_5's table below.
    titles = list(embed.var["title"])
    dim = pd.DataFrame(
        {"rho_S": [spearmanr(E[:, i], S).statistic for i in range(E.shape[1])],
         "rho_G2M": [spearmanr(E[:, i], G).statistic for i in range(E.shape[1])]},
        index=titles)
    dim["cycle_loading"] = dim[["rho_S", "rho_G2M"]].abs().max(axis=1)

    sp = C.read_table("dim_signature_spearman", coll)
    sp.index = [i.strip() for i in sp.index]
    dim = dim.reindex(sp.index)
    # Which axis to put on the y of panel D is decided by the DATA, not hard-coded: the axis
    # whose readouts the cycle explains most. On scie that is stemness by an order of magnitude
    # (mean R2 0.092 against 0.009); on another collection it will name whichever axis is most
    # exposed, which is the one the panel exists to indict.
    listed_only = per_readout[per_readout["has_gene_list"]]
    focus_axis = listed_only.groupby("axis")["r2_cycle_on_score"].mean().idxmax()
    focus = [r for r in coll.by_axis(focus_axis) if r in sp.columns]
    rest = [r for r in sp.columns if r not in focus]
    print(f"\n    panel D focus axis: '{focus_axis}' "
          f"(mean R2 of the cycle on its scores is the highest of the collection)")
    dim["best_focus_axis"] = sp[focus].abs().max(axis=1) if focus else np.nan
    dim["best_other"] = sp[rest].abs().max(axis=1) if rest else np.nan

    conv = C.read_table("convergence", coll)
    conv_dims = set(conv[conv["verdict"] == "convergent"]["dimension"])
    dim["convergent"] = dim.index.isin(conv_dims)

    over = dim.index[dim["cycle_loading"] >= CYCLE_FLAG].tolist()
    print(f"\nD - the latent space: {len(over)} of {len(dim)} dimensions clear the "
          f"{CYCLE_FLAG} cycle bar")
    print(dim.loc[over, ["rho_S", "rho_G2M", "cycle_loading", "best_focus_axis",
                         "best_other", "convergent"]].to_string(float_format="%.3f"))
    rho_stem = spearmanr(dim["cycle_loading"], dim["best_focus_axis"]).statistic
    rho_other = spearmanr(dim["cycle_loading"], dim["best_other"]).statistic
    print(f"  across all {len(dim)} dimensions, Spearman(cycle loading, best {focus_axis} |rho|) "
          f"= {rho_stem:+.3f}")
    print(f"  the same against every readout off that axis             = {rho_other:+.3f}")

    C.write_table(per_readout.round(4), "cycle_confound_by_readout", coll)
    C.write_table(curve.round(4), "cycle_confound_by_vote", coll)
    C.write_table(dim.round(4), "cycle_confound_by_dimension", coll)

    # ---------------------------------------------------------------- figure
    axis_colors = dict(zip(coll.axes, sns.color_palette("deep", len(coll.axes))))
    STEM_C = sns.color_palette("deep")[1]
    CYC_C = "#8E44AD"
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10))

    # --- A: literal cycle-gene content
    ax = axes[0, 0]
    listed = per_readout[per_readout["has_gene_list"]]
    order_a = coll.order(list(listed.index))
    vals = listed.loc[order_a, "cycle_share_of_variance"] * 100
    ax.barh(range(len(order_a)), vals, color=[axis_colors[listed.loc[r, "axis"]] for r in order_a])
    for i, r in enumerate(order_a):
        ax.text(vals.iloc[i] + .18, i, f"{int(listed.loc[r,'n_cycle_genes'])}/{int(listed.loc[r,'n_genes'])}",
                va="center", fontsize=7.5, color="0.35")
    ax.set_yticks(range(len(order_a)))
    ax.set_yticklabels(order_a, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("% of the score's variance carried by Regev cell-cycle genes")
    ax.set_xlim(0, max(vals.max() * 1.28, 1))
    ax.set_title("A · The lists themselves contain cycle genes", fontsize=10.5, loc="left")

    # --- B: content vs variance explained
    ax = axes[0, 1]
    # Everything with no cycle genes piles up on x = 0 with near-identical y. Labelling those
    # points one by one produces an unreadable stack, so the cluster gets ONE annotation naming
    # its members and the per-readout numbers live in cycle_confound_by_readout_*.csv.
    at_zero = [r for r, row in per_readout.iterrows()
               if row["has_gene_list"] and not row["cycle_share_of_variance"] > 0]
    for r, row in per_readout.iterrows():
        listed_r = row["has_gene_list"]
        x = row["cycle_share_of_variance"] * 100 if listed_r else 0.0
        ax.scatter(x, row["r2_cycle_on_score"], s=64 if listed_r else 96,
                   marker="o" if listed_r else "D",
                   facecolor=axis_colors[row["axis"]] if listed_r else "none",
                   edgecolor=axis_colors[row["axis"]], linewidth=1.6, zorder=3)
        if r not in at_zero:
            ax.annotate(r, (x, row["r2_cycle_on_score"]), textcoords="offset points",
                        xytext=(9, 2), fontsize=7.4, color="0.3")
    if at_zero:
        y_hi = per_readout.loc[at_zero, "r2_cycle_on_score"].max()
        n_imm = sum(1 for r in at_zero if coll.axis_of[r] != "stemness")
        label = (f"{n_imm} {'list' if n_imm == 1 else 'lists'} off the stemness axis"
                 + (f"\n+ {', '.join(r for r in at_zero if coll.axis_of[r] == 'stemness')}"
                    if n_imm < len(at_zero) else "")
                 + f"\nno cycle genes at all, R\u00b2 \u2264 {y_hi:.3f}")
        ax.annotate(label, xy=(0.05, y_hi), xytext=(2.1, y_hi + .030),
                    fontsize=7.4, color="0.35", linespacing=1.5,
                    arrowprops=dict(arrowstyle="-", color="0.6", lw=.9,
                                    connectionstyle="arc3,rad=-0.2"))
    ax.set_xlabel("% of variance from cycle genes  (panel A)")
    ax.set_ylabel("R² of (S_score, G2M_score) on the score")
    ax.set_xlim(left=-0.6)
    ax.set_ylim(bottom=-0.008)
    ax.set_title("B · and that is what the cycle explains", fontsize=10.5, loc="left")
    ax.text(.98, .04, "hollow diamond = no gene list,\nso no cycle genes by construction",
            transform=ax.transAxes, ha="right", fontsize=7, color="0.4", linespacing=1.4)

    # --- C: the compounding
    ax = axes[1, 0]
    ax.plot(curve.index, curve["auroc_cycle_predicts_membership"], "-o", color=CYC_C,
            lw=2, ms=6, zorder=3, label="AUROC: cell cycle alone predicts membership")
    lo, hi = min(single.values()), max(single.values())
    ax.axhspan(lo, hi, color=STEM_C, alpha=.16, zorder=1)
    ax.text(len(planes) * .52, (lo + hi) / 2, "range over the single-readout quadrants",
            fontsize=7.2, color="0.35", va="center")
    ax.axhline(0.5, color="0.6", lw=.9, ls=":")
    ax.text(1.02, .503, "no better than chance", fontsize=7, color="0.45")
    for k, v in curve["auroc_cycle_predicts_membership"].items():
        ax.annotate(f"{v:.3f}\nn={curve.loc[k,'n_cells']:,}", (k, v), textcoords="offset points",
                    xytext=(0, 11), ha="center", fontsize=6.8, color="0.3", linespacing=1.3)
    ax2 = ax.twinx()
    ax2.plot(curve.index, curve["pct_G1"], "--s", color="0.45", lw=1.3, ms=4,
             label="% of the called cells in G1")
    ax2.set_ylabel("% in G1", fontsize=9, color="0.45")
    ax2.tick_params(axis="y", labelcolor="0.45")
    ax2.set_ylim(0, 100)
    ax.set_xlabel(f"number of target definitions required to agree "
                  f"(of {len(planes)}: {', '.join(Q.columns[:2])}, ...)"
                  if len(planes) > 2 else "number of target definitions required to agree")
    ax.set_ylabel("AUROC", color=CYC_C)
    ax.tick_params(axis="y", labelcolor=CYC_C)
    ax.set_xticks(list(curve.index))
    ax.set_ylim(.45, 1.0)
    ax.set_title("C · Consensus concentrates the cycle instead of averaging it",
                 fontsize=10.5, loc="left")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=7.4, loc="upper left", framealpha=.92)

    # --- D: the latent space
    ax = axes[1, 1]
    ax.scatter(dim["cycle_loading"], dim["best_other"], s=26, color="0.72",
               label=f"best readout off the {focus_axis} axis", zorder=2)
    ax.scatter(dim["cycle_loading"], dim["best_focus_axis"], s=34, color=STEM_C,
               label=f"best {focus_axis} readout", zorder=3)
    conv_only = dim[dim["convergent"]]
    ax.scatter(conv_only["cycle_loading"], conv_only["best_focus_axis"], s=110,
               facecolor="none", edgecolor=CYC_C, linewidth=1.5, zorder=4,
               label="convergent in 04_7")
    ax.axvline(CYCLE_FLAG, color="crimson", ls="--", lw=1)
    ax.text(CYCLE_FLAG + .012, ax.get_ylim()[1] * .04, f"04_7 cycle flag ({CYCLE_FLAG})",
            fontsize=7, color="crimson", rotation=90, va="bottom")
    for d in over:
        ax.annotate(d, (dim.loc[d, "cycle_loading"], dim.loc[d, "best_focus_axis"]),
                    textcoords="offset points", xytext=(-9, 9), fontsize=8, weight="bold",
                    color=CYC_C, ha="right")
    ax.set_xlabel("cycle loading of the dimension:  max(|ρ S_score|, |ρ G2M_score|)")
    ax.set_ylabel("strongest |ρ| with a readout")
    ax.set_title("D · The cycle dimensions are convergent stemness rows", fontsize=10.5, loc="left")
    ax.legend(fontsize=7.4, loc="upper right", framealpha=.92)
    ax.text(.02, .96, f"ρ(cycle loading, stemness) = {rho_stem:+.3f}\n"
                      f"ρ(cycle loading, other)    = {rho_other:+.3f}",
            transform=ax.transAxes, va="top", fontsize=7.6, family="monospace", color="0.3",
            linespacing=1.5)

    for a in axes.ravel():
        sns.despine(ax=a)
    fig.suptitle(f"{coll.title}: how much of \"stemness\" is the cell cycle", fontsize=12.5)
    fig.tight_layout(rect=[0, .035, 1, .98])
    fig.text(.5, .012,
             "S_score and G2M_score were computed in 04_1 from the same Regev list panel A counts "
             "against, so A and B are not independent measurements.\nCytoTRACE2 has no gene list "
             "and therefore no cycle genes by construction; it is the non-circular control.",
             ha="center", fontsize=7.4, style="italic", color="0.4", linespacing=1.5)
    C.savefig("cycle_behind_stemness", "04_8_cycle_confound", coll, fig, caveat=False)
    plt.close(fig)

    print("\ndone.")


if __name__ == "__main__":
    main()
