#!/usr/bin/env python3
"""04_7, Route C: where the two routes agree. The main result of this stage.

Route A and Route B traverse the same mapping in opposite directions, and neither set they
map between is ground truth. What makes the pair worth running is that their failure modes
do not overlap:

  * a dimension can pass Route A by coincidence among heavily correlated per-cell scores;
  * a dimension can pass Route B by gene-set overlap with no cellular counterpart at all;
  * it is unlikely to pass BOTH for the wrong reason.

So agreement between the routes is the criterion for calling a dimension a genuine cell
state, and disagreement is informative rather than a failure:

  * B but not A  -> the axis carries the gene program but no coherent group of cells sits on
                    it: a candidate patient-specific or technical effect;
  * A but not B  -> the model separates the cells but does not encode the program cleanly on
                    a single axis, so the state is real and the axis is not its description;
  * A and B      -> convergent.

All three categories are reported separately below and NOTHING is promoted on a single route.

One row per dimension AND direction, for every dimension of the run: vanished ones are
not pruned anywhere in this stage (`PRUNE_VANISHED = False` in signature_common), so the
table below is 2 x n_latent rows and a dimension DRVI wrote off can still be read.

Route A is computed per dimension, so its direction is the SIGN of the Spearman
correlation: a signature correlating positively with DR 7 is a statement about DR 7+, and
the same signature correlating negatively is a statement about DR 7-. Route B is already
per direction. That is what makes the two joinable at all.

Usage:
    export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
    python convergence_epi.py
    python convergence_epi.py --rho-min 0.30      # a stricter cell-level bar
"""

from __future__ import annotations

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# signature_common lives in the phase's utils/, as in 02_2_integration.
UTILS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils")
sys.path.insert(0, UTILS_DIR)
import signature_common as C  # noqa: E402

FDR = 0.05

# A readout is flagged as confounded when its raw score correlates with a technical or
# cycle covariate above these. Both are conventions, not derived: chosen so the flag fires
# on the couplings the confounder table actually shows and stays quiet on the rest.
DEPTH_FLAG = 0.30
CYCLE_FLAG = 0.30


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--rho-min", type=float, default=0.20,
                   help="|Spearman rho| above which Route A counts as an association (default 0.20)")
    return p.parse_args()


def read_table(name: str) -> pd.DataFrame:
    return pd.read_csv(C.TABLE_DIR / f"{name}_{C.RUN_ID}.csv", comment="#", index_col=0)


def main():
    args = parse_args()
    C.banner("04_7 - Route C, convergence")

    rho = read_table("dim_signature_spearman")
    eff = read_table("dim_target_effect_size")
    signed = read_table("dim_geneset_signed_significance")
    conf = read_table("confounders")
    order = read_table("dimension_row_order")

    dims = order.index.tolist()
    # Reported, never used as a filter: a dimension DRVI flagged vanished is in the table
    # like any other, and the column is there so a hit on one can be spotted.
    vanished = (order["vanished"].astype(bool) if "vanished" in order.columns
                else pd.Series(False, index=order.index))
    thr = -np.log10(FDR)
    print(f"{len(dims)} dimensions (nothing pruned) x 2 directions = {2 * len(dims)} rows")
    print(f"Route A bar: |rho| >= {args.rho_min};  Route B bar: global FDR < {FDR}")

    # Only the signatures both routes actually carry.
    sigs = [s for s in rho.columns if s in signed.columns]
    independent = [s for s in rho.columns if s not in signed.columns]
    if independent:
        print(f"\nRoute A also carries {', '.join(independent)}, which has no Route B "
              "counterpart by construction (it is a per-cell predictor, not a gene set).")

    rows = []
    for d in dims:
        for direction in ("+", "-"):
            s = 1.0 if direction == "+" else -1.0

            # ---- Route A: the strongest signature association ON THIS SIDE of the axis
            a_vals = rho.loc[d, sigs].astype(float) * s
            a_best = a_vals.idxmax()
            a_rho = float(a_vals.max())

            # the same, over every Route A readout including any without a Route B counterpart
            all_vals = rho.loc[d].astype(float) * s
            a_best_any = all_vals.idxmax()
            a_rho_any = float(all_vals.max())

            # ---- effect size of the target quadrant on this axis, oriented to the side
            auroc = float(eff.loc[d, "auroc_target_vs_rest"])
            smd = float(eff.loc[d, "standardised_mean_difference"]) * s
            auroc_dir = auroc if direction == "+" else 1 - auroc

            # ---- Route B: the strongest enrichment on this side
            b_vals = signed.loc[d, sigs].astype(float) * s
            b_best = b_vals.idxmax()
            b_neglog = float(b_vals.max())
            b_fdr = float(10 ** (-b_neglog)) if b_neglog > 0 else 1.0

            a_hit = a_rho >= args.rho_min
            b_hit = b_neglog >= thr

            same_family = (C.SIG_AXIS.get(a_best) == C.SIG_AXIS.get(b_best)) if (a_hit and b_hit) else False
            same_signature = (a_best == b_best) if (a_hit and b_hit) else False

            if a_hit and b_hit:
                verdict = "convergent" if same_family else "both_routes_different_family"
            elif b_hit:
                verdict = "factor_only_candidate_patient_or_technical"
            elif a_hit:
                verdict = "cell_only_state_not_on_one_axis"
            else:
                verdict = "neither"

            # ---- confounder flags, carried over from A3 for the signature being claimed
            claimed = a_best if a_hit else (b_best if b_hit else a_best)
            depth_rho = float(conf.loc[claimed, "rho_n_genes_by_counts"]) if claimed in conf.index else np.nan
            cyc = max(abs(float(conf.loc[claimed, "rho_S_score"])),
                      abs(float(conf.loc[claimed, "rho_G2M_score"]))) if claimed in conf.index else np.nan
            flags = []
            if not np.isnan(depth_rho) and abs(depth_rho) >= DEPTH_FLAG:
                flags.append("depth")
            if not np.isnan(cyc) and cyc >= CYCLE_FLAG:
                flags.append("cell_cycle")
            if C.SIG_AXIS.get(claimed) == "immune" and a_rho < 0:
                flags.append("immune_low_is_absence_of_signal")

            rows.append({
                "dimension": d, "direction": direction, "dim_direction": f"{d}{direction}",
                "A_best_signature": a_best, "A_rho": a_rho,
                "A_best_any_readout": a_best_any, "A_rho_any": a_rho_any,
                "A_auroc_target_this_side": auroc_dir,
                "A_standardised_mean_difference": smd,
                "B_best_signature": b_best, "B_neglog10_fdr": b_neglog, "B_fdr": b_fdr,
                "A_significant": a_hit, "B_significant": b_hit,
                "same_signature": same_signature, "same_family": same_family,
                "A_family": C.SIG_AXIS.get(a_best), "B_family": C.SIG_AXIS.get(b_best),
                "verdict": verdict,
                "dimension_vanished": bool(vanished[d]),
                "confounder_flags": ",".join(flags) or "none",
            })

    conv = pd.DataFrame(rows).set_index("dim_direction")
    conv = conv.loc[sorted(conv.index, key=C.dim_sort_key)]

    C.banner("the three categories, reported separately")
    counts = conv["verdict"].value_counts()
    print(counts.to_string())
    print("\nNo dimension is promoted on a single route. A 'convergent' row is a candidate\n"
          "cell state; the two single-route categories are candidates for the OTHER thing\n"
          "each of them can be, and are listed here for that reason, not as weaker hits.")

    C.write_table(conv, "convergence")

    conv_rows = conv[conv["verdict"] == "convergent"].sort_values("A_rho", ascending=False)
    print(f"\nCONVERGENT ({len(conv_rows)}): both routes, same signature family")
    if len(conv_rows):
        print(conv_rows[["A_best_signature", "A_rho", "A_auroc_target_this_side",
                         "B_best_signature", "B_fdr", "same_signature",
                         "confounder_flags"]].to_string(float_format="%.3g"))

    b_only = conv[conv["verdict"] == "factor_only_candidate_patient_or_technical"]
    print(f"\nFACTOR-ONLY ({len(b_only)}): the gene program is on the axis, no coherent cell "
          "group is.\nCandidate patient-specific or technical effects - NOT cell states.")
    if len(b_only):
        print(b_only[["B_best_signature", "B_fdr", "A_best_signature", "A_rho",
                      "confounder_flags"]].head(20).to_string(float_format="%.3g"))

    a_only = conv[conv["verdict"] == "cell_only_state_not_on_one_axis"]
    print(f"\nCELL-ONLY ({len(a_only)}): the cells separate, the axis does not encode the "
          "program cleanly.\nThe state may be real; this single dimension is not its description.")
    if len(a_only):
        print(a_only[["A_best_signature", "A_rho", "A_auroc_target_this_side",
                      "B_best_signature", "B_neglog10_fdr", "confounder_flags"]]
              .head(20).to_string(float_format="%.3g"))

    flagged = conv[(conv["verdict"] == "convergent") & (conv["confounder_flags"] != "none")]
    print(f"\n{len(flagged)} of the {len(conv_rows)} convergent rows carry a confounder flag "
          "from A3 and cannot be read as clean.")

    # the project's actual target: immune-evasive AND stem-high on the same axis
    C.banner("the project's target: immune-evasive and stem-high on the same axis")
    tgt = []
    for d in dims:
        for direction in ("+", "-"):
            s = 1.0 if direction == "+" else -1.0
            stem = (rho.loc[d, [x for x in C.STEMNESS_SIGS if x in rho.columns]].astype(float) * s).max()
            imm = (rho.loc[d, C.PRIMARY_IMMUNE].astype(float) * s) * -1.0   # evasion = LOW immunogenicity
            if stem >= args.rho_min and imm >= args.rho_min:
                tgt.append({"dim_direction": f"{d}{direction}", "stem_rho": stem,
                            "immunogenic_low_rho": imm,
                            "auroc_target": conv.loc[f"{d}{direction}", "A_auroc_target_this_side"],
                            "verdict": conv.loc[f"{d}{direction}", "verdict"],
                            "flags": conv.loc[f"{d}{direction}", "confounder_flags"]})
    tgt_df = pd.DataFrame(tgt)
    if len(tgt_df):
        tgt_df = tgt_df.set_index("dim_direction").sort_values("auroc_target", ascending=False)
        print("axes on which stemness is high AND immunogenicity is low at the same time:")
        print(tgt_df.to_string(float_format="%.3f"))
        C.write_table(tgt_df, "target_axes")
    else:
        print(f"No single axis carries both at |rho| >= {args.rho_min}.\n"
              "That is a result: the target state is an INTERSECTION of two axes in the latent\n"
              "space rather than a direction of it, which is what the Route A quadrant already\n"
              "assumed by crossing two independent scores.")

    # ------------------------------------------------------------------ figures
    C.banner("figures")

    sigs_ord = [s for s in C.IMMUNE_SIGS + C.STEMNESS_SIGS if s in sigs]
    fig, axes = plt.subplots(1, 2, figsize=(2 * (0.9 * len(sigs_ord) + 3.5), 0.24 * len(dims) + 3.5),
                             sharey=True)
    sns.heatmap(rho.loc[dims, sigs_ord].astype(float), cmap="vlag", center=0, vmin=-0.6, vmax=0.6,
                cbar_kws={"label": "Spearman rho", "shrink": 0.4}, ax=axes[0])
    axes[0].set_title("Route A - cell-first\ndimension vs within-stratum z-score", fontsize=10)
    vmax = float(np.nanpercentile(signed.loc[dims, sigs_ord].abs().values, 99)) or 1.0
    sns.heatmap(signed.loc[dims, sigs_ord].astype(float), cmap="vlag", center=0,
                vmin=-vmax, vmax=vmax,
                cbar_kws={"label": "signed -log10 FDR", "shrink": 0.4}, ax=axes[1])
    axes[1].set_title("Route B - factor-first\ntop-gene ORA, HVG background", fontsize=10)
    n_imm = len([c for c in C.IMMUNE_SIGS if c in sigs_ord])
    for a in axes:
        a.axvline(n_imm, color="k", lw=1.5)
        plt.setp(a.get_xticklabels(), rotation=45, ha="right", fontsize=8)
    plt.setp(axes[0].get_yticklabels(), fontsize=6)
    fig.suptitle("The two routes side by side, same row order "
                 f"(all {len(dims)} dimensions of {C.RUN_ID}, nothing pruned)\n"
                 "convergence, not either panel alone, is the criterion for a cell state",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    C.savefig("routes_side_by_side", "04_7_convergence", fig)
    plt.close(fig)

    # A strength vs B strength, one point per dimension-direction
    palette = {"convergent": "#C44E52",
               "both_routes_different_family": "#8172B3",
               "factor_only_candidate_patient_or_technical": "#4C72B0",
               "cell_only_state_not_on_one_axis": "#DD8452",
               "neither": "0.8"}
    fig, ax = plt.subplots(figsize=(8, 6.5))
    for v, grp in conv.groupby("verdict"):
        ax.scatter(grp["A_rho"], grp["B_neglog10_fdr"], s=34, lw=0.4, edgecolor="w",
                   c=palette.get(v, "0.5"), label=f"{v} ({len(grp)})")
    ax.axvline(args.rho_min, color="k", ls="--", lw=0.9)
    ax.axhline(thr, color="k", ls="--", lw=0.9)
    for lbl, r in conv[conv["verdict"] == "convergent"].iterrows():
        ax.annotate(lbl, (r["A_rho"], r["B_neglog10_fdr"]), fontsize=6,
                    xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel(f"Route A: strongest signature association on this side (Spearman rho)")
    ax.set_ylabel(f"Route B: strongest enrichment on this side (-log10 global FDR)")
    ax.set_title("Convergence of the two routes, one point per dimension-direction\n"
                 "top-right quadrant = both routes; only its same-family members are "
                 "called cell states", fontsize=10)
    ax.legend(fontsize=7, loc="upper left", frameon=False)
    sns.despine(ax=ax)
    C.savefig("convergence_scatter", "04_7_convergence", fig)
    plt.close(fig)

    print("\ndone.")


if __name__ == "__main__":
    main()
