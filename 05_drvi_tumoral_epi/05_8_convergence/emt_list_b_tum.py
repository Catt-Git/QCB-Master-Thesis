#!/usr/bin/env python3
"""05_6 / 05_7 / 05_8 EMT figures redrawn on collaborator list B only.

The emt collection carries three list versions side by side (A, B, C). This redraws the
three readouts of the chain - Route A (cell-first), Route B (factor-first) and their
intersection (the two routes side by side and the convergence scatter) - keeping only the
primary triad, `EMT_B_EPITHELIAL`, `EMT_B_HYBRID`, `EMT_B_MESENCHYMAL`.

NOTHING IS RECOMPUTED UPSTREAM. The script reads the emt tables 05_6 and 05_7 already wrote
and drops the A and C columns:

  * Route A is a per-signature Spearman rho, so a column does not depend on the others;
  * Route B keeps the GLOBAL BH of 05_7 (every direction x every emt list + Hallmark). A BH
    over list B alone would have a smaller denominator and call slightly more pairs; keeping
    the original FDR makes these figures directly comparable with the full-collection ones;
  * the convergence verdict IS recomputed, because "the strongest signature on this side"
    now ranges over the three B lists only. Same rules and bars as convergence_tum.py. With
    one list per axis "same family" and "same signature" are the same thing.

The full-collection figures are not touched: everything goes under a new collection slug,
`emt_B`, i.e. figures/<step>/emt_B/<run_id>/ and tables/emt_B/<run_id>/.

Usage:
    export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
    N_LATENT=64 python emt_list_b_tum.py              # rho bar 0.20, as the reported run
    N_LATENT=64 python emt_list_b_tum.py --rho-min 0.30
"""
from __future__ import annotations

import argparse
import dataclasses
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

UTILS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils")
sys.path.insert(0, UTILS_DIR)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import signature_common as C  # noqa: E402
import sig_collections as SC  # noqa: E402
from convergence_tum import (FDR, CONVERGENT, BOTH_DIFFERENT, FACTOR_ONLY, CELL_ONLY,  # noqa: E402
                             NEITHER, VERDICT_LABEL, DEPTH_FLAG, CYCLE_FLAG)

VERSION = "B"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--rho-min", type=float, default=C.ROUTE_A_RHO_MIN,
                   help=f"Route A bar, as in convergence_tum.py (default {C.ROUTE_A_RHO_MIN})")
    return p.parse_args()


def main():
    args = parse_args()
    src = SC.EMT
    # The same collection restricted to list B, under its own slug so that every path it
    # writes to is new. The tables are still READ through `src`.
    coll = dataclasses.replace(
        src, name=f"emt_{VERSION}", title=f"{src.title}, list {VERSION} only",
        signatures=tuple(s for s in src.signatures if SC._emt_version(s.name) == VERSION),
        derived=tuple(d for d in src.derived if SC._emt_version(d.name) == VERSION))
    sigs = coll.order(coll.names, for_figure=True)
    C.banner(f"EMT, list {VERSION} only: {', '.join(sigs)}")

    rho = C.read_table("dim_signature_spearman", src)
    signed = C.read_table("dim_geneset_signed_significance", src)
    conf = C.read_table("confounders", src)
    eff = C.read_table("dim_target_effect_size", src)
    order = C.read_table("dimension_row_order", src)
    dims = order.index.tolist()
    thr = -np.log10(FDR)
    edges = coll.block_edges(sigs)
    pruning = "vanished PRUNED" if C.PRUNE_VANISHED else "nothing pruned"

    # ------------------------------------------------ Route A, cell-first
    fig, ax = plt.subplots(figsize=(C.fig_span(len(sigs), 1.0, 4.0),
                                    C.fig_span(len(dims), 0.24, 3.0)))
    sns.heatmap(rho.loc[dims, sigs].astype(float), cmap="vlag", center=0, vmin=-0.6, vmax=0.6,
                cbar_kws={"label": "Spearman rho (dimension vs within-stratum z-score)\n"
                                   "sign = direction: rho > 0 is DR n+, rho < 0 is DR n-",
                          "shrink": 0.4}, ax=ax)
    ax.set_title(f"Route A, {coll.title}: dimensions x signatures\n"
                 f"{len(dims)} dimensions of {C.RUN_ID}, {pruning}\n"
                 "rows are dimensions: the sign of rho IS the direction", fontsize=10)
    for pos in edges:
        ax.axvline(pos, color="k", lw=1.5)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=6)
    C.savefig("dim_signature_heatmap", "05_6_cell_first", coll, fig)
    plt.close(fig)

    # ------------------------------------------------ Route B, factor-first
    vmax = float(np.nanpercentile(signed.loc[dims, sigs].abs().values, 99)) or 1.0
    fig, ax = plt.subplots(figsize=(C.fig_span(len(sigs), 1.0, 4.0),
                                    C.fig_span(len(dims), 0.24, 3.0)))
    sns.heatmap(signed.loc[dims, sigs].astype(float), cmap="vlag", center=0, vmin=-vmax, vmax=vmax,
                cbar_kws={"label": "signed -log10 FDR  (+ = positive direction)", "shrink": 0.4},
                ax=ax)
    for pos in edges:
        ax.axvline(pos, color="k", lw=1.5)
    ax.set_title(f"Route B, {coll.title}: latent dimensions x gene sets, signed significance\n"
                 f"{len(dims)} dimensions of {C.RUN_ID}, {pruning}; FDR = the global BH of 05_7 "
                 f"over the full emt collection + Hallmark\n"
                 f"sign = the direction of the axis carrying the enrichment; "
                 f"|value| >= {thr:.2f} is FDR < {FDR}", fontsize=9)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=6)
    C.savefig("dim_geneset_signed_heatmap", "05_7_factor_first", coll, fig)
    plt.close(fig)

    # ------------------------------------------------ intersection: verdicts on B only
    rows = []
    for d in dims:
        for direction in ("+", "-"):
            s = 1.0 if direction == "+" else -1.0
            a_vals = rho.loc[d, sigs].astype(float) * s
            a_best, a_rho = a_vals.idxmax(), float(a_vals.max())
            b_vals = signed.loc[d, sigs].astype(float) * s
            b_best, b_neglog = b_vals.idxmax(), float(b_vals.max())
            auroc = float(eff.loc[d, "auroc_target_vs_rest"])

            a_hit, b_hit = a_rho >= args.rho_min, b_neglog >= thr
            same_family = (coll.axis_of[a_best] == coll.axis_of[b_best]) if (a_hit and b_hit) else False
            if a_hit and b_hit:
                verdict = CONVERGENT if same_family else BOTH_DIFFERENT
            elif b_hit:
                verdict = FACTOR_ONLY
            elif a_hit:
                verdict = CELL_ONLY
            else:
                verdict = NEITHER

            claimed = a_best if a_hit else (b_best if b_hit else a_best)
            flags = []
            if claimed in conf.index:
                if abs(float(conf.loc[claimed, "rho_n_genes_by_counts"])) >= DEPTH_FLAG:
                    flags.append("depth")
                if max(abs(float(conf.loc[claimed, "rho_S_score"])),
                       abs(float(conf.loc[claimed, "rho_G2M_score"]))) >= CYCLE_FLAG:
                    flags.append("cell_cycle")
            flags += coll.extra_flags(coll, claimed, float(rho.loc[d, claimed]))

            rows.append({
                "dimension": d, "direction": direction, "dim_direction": f"{d}{direction}",
                "A_best_signature": a_best, "A_rho": a_rho,
                "A_auroc_target_this_side": auroc if direction == "+" else 1 - auroc,
                "A_standardised_mean_difference": float(eff.loc[d, "standardised_mean_difference"]) * s,
                "B_best_signature": b_best, "B_neglog10_fdr": b_neglog,
                "B_fdr": float(10 ** (-b_neglog)) if b_neglog > 0 else 1.0,
                "A_significant": a_hit, "B_significant": b_hit,
                "same_family": same_family, "verdict": verdict,
                "confounder_flags": ",".join(flags) or "none",
            })
    conv = pd.DataFrame(rows).set_index("dim_direction")
    conv = conv.loc[sorted(conv.index, key=C.dim_sort_key)]
    C.write_table(conv, "convergence", coll)

    counts = conv["verdict"].value_counts()
    print(f"\nRoute A bar |rho| >= {args.rho_min}; Route B bar global FDR < {FDR}")
    for v in (CONVERGENT, BOTH_DIFFERENT, FACTOR_ONLY, CELL_ONLY, NEITHER):
        print(f"  {VERDICT_LABEL[v]:32s} {int(counts.get(v, 0)):3d}")
    hits = conv[conv["verdict"] != NEITHER].sort_values(["verdict", "A_rho"], ascending=[True, False])
    if len(hits):
        print(hits[["verdict", "A_best_signature", "A_rho", "B_best_signature", "B_fdr",
                    "confounder_flags"]].to_string(float_format="%.3g"))

    # side by side, same row order
    panel_w = C.fig_span(len(sigs), 0.9, 3.5, cap=C.MAX_FIG_IN / 2)
    fig, axes = plt.subplots(1, 2, figsize=(2 * panel_w, C.fig_span(len(dims), 0.24, 3.5)),
                             sharey=True)
    sns.heatmap(rho.loc[dims, sigs].astype(float), cmap="vlag", center=0, vmin=-0.6, vmax=0.6,
                cbar_kws={"label": "Spearman rho", "shrink": 0.4}, ax=axes[0])
    axes[0].set_title("Route A - cell-first\ndimension vs within-stratum z-score", fontsize=10)
    sns.heatmap(signed.loc[dims, sigs].astype(float), cmap="vlag", center=0, vmin=-vmax, vmax=vmax,
                cbar_kws={"label": "signed -log10 FDR", "shrink": 0.4}, ax=axes[1])
    axes[1].set_title("Route B - factor-first\ntop-gene ORA, HVG background", fontsize=10)
    for a in axes:
        for pos in edges:
            a.axvline(pos, color="k", lw=1.5)
        plt.setp(a.get_xticklabels(), rotation=45, ha="right", fontsize=8)
    plt.setp(axes[0].get_yticklabels(), fontsize=6)
    fig.suptitle(f"The two routes side by side, {coll.title}, same row order "
                 f"({len(dims)} dimensions of {C.RUN_ID}, {pruning})\n"
                 "convergence, not either panel alone, is the criterion for a cell state",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    C.savefig("routes_side_by_side", "05_8_convergence", coll, fig)
    plt.close(fig)

    # A strength vs B strength
    palette = {CONVERGENT: "#C44E52", BOTH_DIFFERENT: "#8172B3", FACTOR_ONLY: "#4C72B0",
               CELL_ONLY: "#DD8452", NEITHER: "0.8"}
    fig, ax = plt.subplots(figsize=(8, 6.5))
    for v, grp in conv.groupby("verdict"):
        ax.scatter(grp["A_rho"], grp["B_neglog10_fdr"], s=34, lw=0.4, edgecolor="w",
                   c=palette.get(v, "0.5"), label=f"{VERDICT_LABEL.get(v, v)} ({len(grp)})")
    ax.axvline(args.rho_min, color="k", ls="--", lw=0.9)
    ax.axhline(thr, color="k", ls="--", lw=0.9)
    for lbl, r in conv[conv["verdict"] == CONVERGENT].iterrows():
        ax.annotate(lbl, (r["A_rho"], r["B_neglog10_fdr"]), fontsize=6,
                    xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("Route A: strongest list-B association on this side (Spearman rho)")
    ax.set_ylabel("Route B: strongest list-B enrichment on this side (-log10 global FDR)")
    ax.set_title(f"Convergence of the two routes, {coll.title}\n"
                 "one point per dimension-direction; top-right quadrant = both routes; "
                 "only its same-family members are called cell states", fontsize=10)
    ax.legend(fontsize=7, loc="upper left", frameon=False)
    sns.despine(ax=ax)
    C.savefig("convergence_scatter", "05_8_convergence", coll, fig)
    plt.close(fig)

    print("\ndone.")


if __name__ == "__main__":
    main()
