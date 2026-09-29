#!/usr/bin/env python3
"""05_11: do the EMT dimension calls survive a split of the cells by treatment?

The relevant dimension-directions are the ones the EMT analysis did not call `neither`:
  * list B only (05_8_convergence/emt_list_b_tum.py, tables/emt_B/): 14 rows, all single-route;
  * the full emt collection (convergence_tum.py, tables/emt/): the 3 two-route rows, which
    list B alone loses (DR 1+, DR 21+, DR 32-).

WHAT CAN CHANGE AND WHAT CANNOT. Route B, factor-first, reads the DRVI decoder's gene loadings
and never looks at a cell, so splitting the cells cannot move it: its value is carried over
unchanged. Route A, cell-first, is a Spearman rho over cells and CAN change, so it is
recomputed on each treatment arm (BASE, PD1, RTPD1) exactly as 05_6 computes it: raw
signature score -> z within `cohort` (now within cohort INSIDE the arm) -> Spearman against
the latent coordinate.

WHY THE PER-PATIENT LAYER. Treatment is not randomised across cells: it is a timepoint of a
patient, and the arms are made of different patients in different proportions (Patient16 alone
is 7,100 of the 25,487 PD1 cells and has 583 BASE cells). A pooled rho that moves between arms
can therefore be a change of patient mix rather than of biology. So each arm is also cut per
patient (strata of >= MIN_CELLS cells), and an arm difference is tested only between arms of
the SAME patient: a Wilcoxon signed-rank test on the paired per-patient rho, BH across every
row x arm pair. The pooled numbers describe; the paired test is the one that decides.

A side table does the same for the POSITION of the cells on each axis (per-patient mean of the
latent coordinate, BASE vs PD1 vs RTPD1, paired), because "the association changes" and "the
cells move along the axis" are different claims and the second one is the likelier one.

Usage:
    export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
    N_LATENT=64 python treatment_split_tum.py
"""
from __future__ import annotations

import dataclasses
import itertools
import os
import sys

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr, wilcoxon
from statsmodels.stats.multitest import multipletests

UTILS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils")
sys.path.insert(0, UTILS_DIR)
import signature_common as C  # noqa: E402
import sig_collections as SC  # noqa: E402

ARMS = ("BASE", "PD1", "RTPD1")
MIN_CELLS = 100      # smallest patient x arm stratum given a rho of its own
FDR = 0.05
STEP = "05_11_treatment_split"


def zscore_within(raw: pd.DataFrame, groups: pd.Series) -> pd.DataFrame:
    """05_6's standardise_within, on an arbitrary grouping: flat strata become 0."""
    g = raw.groupby(groups.values, observed=True)
    return g.transform(lambda s: (s - s.mean()) / s.std(ddof=0) if s.std(ddof=0) > 0 else s * 0.0).fillna(0.0)


def route_a(L: pd.DataFrame, z: pd.DataFrame) -> pd.DataFrame:
    """dims x readouts Spearman rho, as 05_6 computes it; NaN (no spread) -> 0."""
    rho = pd.DataFrame(index=L.columns, columns=z.columns, dtype=float)
    for r in z.columns:
        rr, _ = spearmanr(L.values, z[r].values)
        rho[r] = rr[:-1, -1]
    return rho.fillna(0.0)


def verdict(a_rho: float, b_neglog: float, a_best: str, b_best: str, coll, rho_min: float) -> str:
    a_hit, b_hit = a_rho >= rho_min, b_neglog >= -np.log10(FDR)
    if a_hit and b_hit:
        return "convergent" if coll.axis_of[a_best] == coll.axis_of[b_best] else "both_different_family"
    return "factor_only" if b_hit else ("cell_only" if a_hit else "neither")


def main():
    emt = SC.EMT
    emt_b = dataclasses.replace(
        emt, name="emt_B", title=f"{emt.title}, list B only",
        signatures=tuple(s for s in emt.signatures if SC._emt_version(s.name) == "B"),
        derived=tuple(d for d in emt.derived if SC._emt_version(d.name) == "B"))
    rho_min = C.ROUTE_A_RHO_MIN
    C.banner("05_11 - EMT dimension calls split by treatment")

    # --------------------------------------------------------- the relevant rows
    sources = {"list B": emt_b, "full emt": emt}
    rel = []
    for label, coll in sources.items():
        conv = C.read_table("convergence", coll)
        keep = conv[conv["verdict"] != "neither"]
        if label == "full emt":   # only what list B does not already carry
            keep = keep[keep["verdict"].isin(["convergent", "both_different_family"])]
        for dd, r in keep.iterrows():
            claimed = r["A_best_signature"] if r["A_significant"] else r["B_best_signature"]
            rel.append({"dim_direction": dd, "dimension": r["dimension"], "direction": r["direction"],
                        "source": label, "verdict_all": r["verdict"], "claimed": claimed})
    rel = pd.DataFrame(rel).set_index("dim_direction")
    rel["row"] = [f"{d} · {c.replace('EMT_', '')}" for d, c in zip(rel.index, rel["claimed"])]
    print(f"{len(rel)} relevant dimension-directions:")
    print(rel[["source", "verdict_all", "claimed"]].to_string())

    # ------------------------------------------------------------------ cells
    embed = ad.read_h5ad(C.EMBED_H5AD)
    obs = embed.obs[["cohort", "treatment"]].astype(str)
    dims = sorted(set(rel["dimension"]), key=lambda d: int(d[3:]))
    L = pd.DataFrame(np.asarray(embed.X), index=embed.obs_names,
                     columns=embed.var["title"].values)[dims]
    scores = pd.read_csv(C.scores_csv(emt), index_col=0)
    raw = scores[[c for c in scores.columns if c.startswith("score_")]]
    raw.columns = [c[len("score_"):] for c in raw.columns]
    assert set(raw.index) == set(L.index), "scores and embedding disagree on the cells"
    raw = raw.loc[L.index]
    print("\ncells per arm:", obs["treatment"].value_counts().to_dict())

    signed = C.read_table("dim_geneset_signed_significance", emt)

    # ---------------------------------------------------- pooled rho per arm
    # 'ALL' is recomputed through the same function, so it must reproduce 05_6: checked below.
    pooled = {}
    for arm in ("ALL",) + ARMS:
        m = np.ones(len(obs), bool) if arm == "ALL" else (obs["treatment"] == arm).values
        z = zscore_within(raw[m], obs.loc[m, "cohort"])
        pooled[arm] = route_a(L[m], z)
    ref = C.read_table("dim_signature_spearman", emt)
    # The nine lists only: the derived EMT_SCORE_* columns get their z in 05_6 as a contrast of
    # the parents' z, not by re-standardising the raw contrast, and no relevant row claims one.
    dev = (pooled["ALL"][emt.names] - ref.loc[dims, emt.names]).abs().max().max()
    print(f"ALL reproduces 05_6's dim_signature_spearman to within {dev:.1e}")
    assert dev < 1e-3

    rows = []
    for dd, r in rel.iterrows():
        s = 1.0 if r["direction"] == "+" else -1.0
        coll = sources[r["source"]]
        sigs = [n for n in coll.names if n in signed.columns]
        b_vals = signed.loc[r["dimension"], sigs].astype(float) * s
        out = {"dim_direction": dd, "source": r["source"], "claimed": r["claimed"],
               "B_neglog10_fdr": float(b_vals.max()), "B_best": b_vals.idxmax()}
        for arm, rho in pooled.items():
            a_vals = rho.loc[r["dimension"], sigs] * s
            out[f"rho_claimed_{arm}"] = float(rho.loc[r["dimension"], r["claimed"]] * s)
            out[f"A_best_{arm}"] = a_vals.idxmax()
            out[f"verdict_{arm}"] = verdict(float(a_vals.max()), out["B_neglog10_fdr"],
                                            a_vals.idxmax(), out["B_best"], coll, rho_min)
        out["n_cells"] = ", ".join(f"{a}={int((obs['treatment'] == a).sum())}" for a in ARMS)
        rows.append(out)
    pooled_tbl = pd.DataFrame(rows).set_index("dim_direction")
    C.write_table(pooled_tbl, "treatment_split_pooled", emt_b)

    C.banner("pooled, per arm: rho of the claimed signature, oriented to the side")
    show = pooled_tbl[[f"rho_claimed_{a}" for a in ("ALL",) + ARMS]
                      + [f"verdict_{a}" for a in ("ALL",) + ARMS]]
    show.columns = [c.replace("rho_claimed_", "rho_").replace("verdict_", "v_") for c in show.columns]
    print(show.to_string(float_format="%.3f"))
    changed = pooled_tbl[pooled_tbl[[f"verdict_{a}" for a in ARMS]].ne(pooled_tbl["verdict_ALL"], axis=0).any(axis=1)]
    print(f"\n{len(changed)} of {len(pooled_tbl)} rows change verdict in at least one arm")

    # -------------------------------------------- per patient x arm, paired
    strata = obs.groupby(["cohort", "treatment"], observed=True).size()
    strata = strata[strata >= MIN_CELLS]
    print(f"\n{len(strata)} patient x arm strata with >= {MIN_CELLS} cells")
    claimed_sigs = sorted(set(rel["claimed"]))
    per = []
    for (pat, arm), n in strata.items():
        m = ((obs["cohort"] == pat) & (obs["treatment"] == arm)).values
        zz = zscore_within(raw.loc[m, claimed_sigs], obs.loc[m, "cohort"])
        rho = route_a(L[m], zz)
        pos = L[m].mean()
        for dd, r in rel.iterrows():
            s = 1.0 if r["direction"] == "+" else -1.0
            per.append({"dim_direction": dd, "cohort": pat, "treatment": arm, "n_cells": int(n),
                        "rho_claimed": float(rho.loc[r["dimension"], r["claimed"]] * s),
                        "mean_latent_oriented": float(pos[r["dimension"]] * s)})
    per = pd.DataFrame(per)
    C.write_table(per, "treatment_split_per_patient", emt_b, index=False)

    tests = []
    for dd in rel.index:
        sub = per[per["dim_direction"] == dd]
        for quantity in ("rho_claimed", "mean_latent_oriented"):
            wide = sub.pivot(index="cohort", columns="treatment", values=quantity)
            for a1, a2 in itertools.combinations(ARMS, 2):
                if a1 not in wide or a2 not in wide:
                    continue
                pair = wide[[a1, a2]].dropna()
                delta = pair[a2] - pair[a1]
                p = float(wilcoxon(delta).pvalue) if len(pair) >= 5 and delta.abs().sum() > 0 else np.nan
                tests.append({"dim_direction": dd, "quantity": quantity, "contrast": f"{a2} - {a1}",
                              "n_patients": len(pair), "median_delta": float(delta.median()) if len(pair) else np.nan,
                              "n_up": int((delta > 0).sum()), "n_down": int((delta < 0).sum()),
                              "p_wilcoxon": p})
    tests = pd.DataFrame(tests)
    for q in tests["quantity"].unique():
        m = (tests["quantity"] == q) & tests["p_wilcoxon"].notna()
        tests.loc[m, "fdr_bh"] = multipletests(tests.loc[m, "p_wilcoxon"], method="fdr_bh")[1]
    tests["significant"] = tests["fdr_bh"] < FDR
    C.write_table(tests, "treatment_split_paired_tests", emt_b, index=False)

    C.banner("paired within-patient contrasts (Wilcoxon, BH per quantity)")
    for q in ("rho_claimed", "mean_latent_oriented"):
        t = tests[tests["quantity"] == q]
        print(f"\n{q}: {int(t['significant'].sum())} of {int(t['p_wilcoxon'].notna().sum())} "
              f"contrasts significant at FDR < {FDR}")
        print(t.sort_values("p_wilcoxon").head(12)[["dim_direction", "contrast", "n_patients",
              "median_delta", "n_up", "n_down", "p_wilcoxon", "fdr_bh"]].to_string(index=False, float_format="%.3g"))

    # ------------------------------------------------------------------ figures
    order = rel.index.tolist()
    labels = rel["row"].tolist()

    # 1. pooled rho per arm, with the verdict it gives
    mat = pooled_tbl.loc[order, [f"rho_claimed_{a}" for a in ("ALL",) + ARMS]]
    mat.columns = ["all cells"] + list(ARMS)
    mat.index = labels
    fig, ax = plt.subplots(figsize=(6.5, C.fig_span(len(order), 0.38, 2.5)))
    sns.heatmap(mat, cmap="vlag", center=0, vmin=-0.5, vmax=0.5, annot=True, fmt=".2f",
                annot_kws={"fontsize": 7}, cbar_kws={"label": "Spearman rho, oriented to the side",
                                                     "shrink": 0.5}, ax=ax)
    for i, dd in enumerate(order):
        for j, a in enumerate(("ALL",) + ARMS):
            if pooled_tbl.loc[dd, f"verdict_{a}"] != pooled_tbl.loc[dd, "verdict_ALL"]:
                ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=False, ec="k", lw=1.6))
    ax.axvline(1, color="k", lw=1.2)
    n_list_b = int((rel["source"] == "list B").sum())
    ax.axhline(n_list_b, color="k", lw=1.2, ls="--")
    ax.set_title("Route A per treatment arm, EMT dimension calls\n"
                 f"rho of the claimed signature; boxed = verdict differs from all cells "
                 f"(bar {rho_min:.2f}, Route B unchanged)\n"
                 f"above the dashed line: list B rows; below: full-collection convergent rows",
                 fontsize=9)
    ax.set_ylabel("")
    plt.setp(ax.get_yticklabels(), fontsize=7)
    C.savefig("pooled_rho_by_treatment", STEP, emt_b, fig)
    plt.close(fig)

    # 2. per patient, one line per patient across the arms
    ncol = 4
    nrow = int(np.ceil(len(order) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.1 * ncol, 2.5 * nrow), sharey=True)
    xs = {a: i for i, a in enumerate(ARMS)}
    pats = sorted(per["cohort"].unique())
    pal = dict(zip(pats, sns.color_palette("husl", len(pats))))
    for ax, dd, lab in zip(axes.flat, order, labels):
        sub = per[per["dim_direction"] == dd]
        for pat, g in sub.groupby("cohort"):
            g = g.sort_values("treatment", key=lambda s: s.map(xs))
            ax.plot(g["treatment"].map(xs), g["rho_claimed"], "-o", ms=3, lw=0.8,
                    color=pal[pat], alpha=0.8)
        for a in ARMS:
            ax.hlines(pooled_tbl.loc[dd, f"rho_claimed_{a}"], xs[a] - 0.25, xs[a] + 0.25,
                      color="k", lw=2)
        ax.axhline(rho_min, color="0.5", ls="--", lw=0.7)
        ax.axhline(0, color="0.8", lw=0.6)
        ax.set_xticks(range(len(ARMS)), ARMS, fontsize=7)
        ax.set_xlim(-0.4, len(ARMS) - 0.6)
        ax.set_title(lab, fontsize=8)
        sns.despine(ax=ax)
    for ax in axes.flat[len(order):]:
        ax.axis("off")
    for ax in axes[:, 0]:
        ax.set_ylabel("rho, oriented", fontsize=8)
    fig.suptitle(f"Per-patient Route A rho by treatment arm (patient x arm strata >= {MIN_CELLS} cells; "
                 "one colour per patient, lines join the same patient;\nblack bar = pooled arm, "
                 f"dashed = the {rho_min:.2f} bar)", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    C.savefig("per_patient_rho_by_treatment", STEP, emt_b, fig)
    plt.close(fig)

    print("\ndone.")


if __name__ == "__main__":
    main()
