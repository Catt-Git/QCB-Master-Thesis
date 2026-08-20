#!/usr/bin/env python3
"""04_5, Route A: cell-first. Prior knowledge defines the state, DRVI is the coordinate system.

The unit of analysis is the CELL. The question is whether the cells that prior knowledge
calls stem-like or immune-evasive occupy a distinct position along any latent dimension.
This is the only route that assigns cells to states - Route B can tell you a dimension is
enriched for stemness genes, and still not tell you which cells, how many, or in which
patients - so the deliverable of the project, an actual cell assignment, comes from here.

How it fails, and what this script does about it:

  * per-cell scores are noisy and correlate with depth, with the cycle, and with each
    other -> the confounder table, and the two named checks below;
  * the stemness lists are embryonic and proliferation-heavy while this compartment has a
    large Lumsec-prol population, so "stem-high" can just mean "cycling" -> the quadrant is
    recomputed inside G1 alone and the two cell sets compared;
  * immune evasion is defined by ABSENCE of signal, which shallow sequencing mimics
    perfectly -> the depth of the immunogenic-low group is compared against the rest;
  * absolute scores are not comparable across patients sequenced at different depths ->
    every readout is standardised within (cohort, cell_type) before it is used;
  * it is confirmatory by construction and can only find what was brought in from outside.
    That one has no fix inside Route A; it is why Route B exists.

Scoring is on the UNINTEGRATED, ALL-GENES object. A 150-gene signature reduced to whatever
survived HVG selection is no longer that signature, so `shiao_epi.h5ad` is read here and
never `shiao_epi_hvg_2k.h5ad`.

Every latent dimension is correlated, vanished ones included (`PRUNE_VANISHED = False` in
signature_common). Dropping dimensions before the correlations would decide ahead of the
analysis which axes are allowed to mean something; the vanished flag is reported next to
the results instead.

No stemness signature is primary: the quadrant is defined once per stemness readout and the
stability of the resulting cell set across those definitions is itself a reported result.

Usage:
    export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
    python cell_first_epi.py
    python cell_first_epi.py --high-q 0.80 --low-q 0.20   # a stricter quadrant
    python cell_first_epi.py --overwrite                  # re-score instead of reusing the csv
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
import scanpy as sc
import seaborn as sns
from scipy.stats import mannwhitneyu, spearmanr
from sklearn.metrics import roc_auc_score

# signature_common lives in the phase's utils/, as in 02_2_integration.
UTILS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils")
sys.path.insert(0, UTILS_DIR)
import signature_common as C  # noqa: E402

# Non-default scoring parameters, all on the record for the Methods.
N_BINS = 25                 # sc.tl.score_genes default is 25; stated because it is not optional
GROUPBY = ["cohort", "cell_type"]   # the standardisation strata


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--high-q", type=float, default=0.75,
                   help="quantile of the within-stratum z-score above which a cell is 'high' (default 0.75)")
    p.add_argument("--low-q", type=float, default=0.25,
                   help="quantile below which a cell is 'low' (default 0.25)")
    p.add_argument("--overwrite", action="store_true", help="re-score even if the per-cell csv exists")
    return p.parse_args()


# --------------------------------------------------------------------------- #
# Scoring
# --------------------------------------------------------------------------- #

def score_signatures(adata, sets: dict[str, list[str]]) -> tuple[pd.DataFrame, list[str]]:
    """One sc.tl.score_genes call per signature, on the all-genes log-normalised object.

    `ctrl_size` is matched to the signature length, as the scanpy documentation intends:
    the control set is drawn from the same expression bins, so a 25-gene signature scored
    against 50 control genes is not the same statistic as a 400-gene one.

    The control genes are sampled AT RANDOM, so `random_state` is not optional - without
    it the whole step is unreproducible.
    """
    scores, skipped = {}, []
    for name, genes in sets.items():
        if len(genes) < C.MIN_SIGNATURE_GENES:
            print(f"  {name:24s} SKIPPED, only {len(genes)} mapped genes "
                  f"(floor {C.MIN_SIGNATURE_GENES})")
            skipped.append(name)
            continue
        sc.tl.score_genes(
            adata, gene_list=genes, ctrl_size=len(genes), n_bins=N_BINS,
            random_state=C.SEED, use_raw=False, score_name=f"score_{name}",
        )
        v = adata.obs[f"score_{name}"].astype(float)
        scores[name] = v.values
        print(f"  {name:24s} {len(genes):5d} genes  mean {v.mean():+.4f}  sd {v.std():.4f}")
    return pd.DataFrame(scores, index=adata.obs_names), skipped


def standardise_within(df: pd.DataFrame, strata: pd.DataFrame) -> pd.DataFrame:
    """Z-score each column within (cohort, cell_type).

    Absolute scores are not comparable across patients sequenced at different depths, nor
    across cell types with different baseline expression of the signature. A stratum whose
    scores have no spread (a single cell, say) would divide by zero; those become 0, i.e.
    'exactly average for its stratum', which is the only honest value for a stratum of one.
    """
    g = df.groupby([strata[k].values for k in GROUPBY], observed=True)
    z = g.transform(lambda s: (s - s.mean()) / s.std(ddof=0) if s.std(ddof=0) > 0 else s * 0.0)
    return z.fillna(0.0).add_prefix("z_")


# --------------------------------------------------------------------------- #
# Quadrants
# --------------------------------------------------------------------------- #

def define_quadrant(z_stem: pd.Series, z_imm: pd.Series, high_q: float, low_q: float) -> pd.Series:
    """stem-high / immunogenic-low, on within-stratum z-scores.

    The cutoffs are quantiles of the standardised scores rather than fixed z values: the
    scores are not normal and a fixed z would give wildly different group sizes across
    readouts, which would make the stability comparison below meaningless.
    """
    return (z_stem >= z_stem.quantile(high_q)) & (z_imm <= z_imm.quantile(low_q))


def jaccard(a: pd.Series, b: pd.Series) -> float:
    inter = int((a & b).sum())
    union = int((a | b).sum())
    return inter / union if union else np.nan


def main():
    args = parse_args()
    C.banner("04_5 - Route A, cell-first")

    sets = C.read_gmt(C.GMT_PATH)
    print(f"{len(sets)} signatures from {C.GMT_PATH}")

    # ------------------------------------------------------------------ input
    adata = ad.read_h5ad(C.FULL_H5AD)
    print(f"{adata.n_obs:,} cells x {adata.n_vars:,} genes (all genes, unintegrated)")

    # `.X` must be log-normalised: score_genes on raw counts is a different statistic.
    xmax = float(adata.X.max())
    assert xmax < 50, (f".X looks like counts (max {xmax:.1f}), not log-normalised data. "
                       "Route A must score on the scran log1p matrix.")
    print(f".X max {xmax:.3f}, log1p in uns: {'log1p' in adata.uns}  -> log-normalised")

    for k in ("S_score", "G2M_score", "phase", "n_genes_by_counts", "pct_counts_mt"):
        assert k in adata.obs, f"missing obs key: {k}"
    print("cell cycle scores reused from 04_1, not recomputed and not regressed out")

    # ----------------------------------------------------------------- scores
    if C.SCORES_CSV.exists() and not args.overwrite:
        print(f"\n[read] {C.SCORES_CSV} (--overwrite to re-score)")
        allsc = pd.read_csv(C.SCORES_CSV, index_col=0)
        assert (allsc.index == adata.obs_names).all(), "score csv is not in the object's cell order"
        raw = allsc[[c for c in allsc.columns if not c.startswith("z_")]]
        raw.columns = [c.replace("score_", "") for c in raw.columns]
        skipped = [n for n in sets if n not in raw.columns and n != "CytoTRACE2"]
    else:
        print("\nscoring (sc.tl.score_genes, use_raw=False, ctrl_size=len(signature), "
              f"n_bins={N_BINS}, random_state={C.SEED})")
        raw, skipped = score_signatures(adata, sets)

    # ------------------------------------------------------------ CytoTRACE2
    ct2_col = None
    if C.CYTOTRACE_CSV.exists():
        ct2 = pd.read_csv(C.CYTOTRACE_CSV, index_col=0).reindex(adata.obs_names)
        if "CytoTRACE2_Score" in ct2:
            raw["CytoTRACE2"] = ct2["CytoTRACE2_Score"].values
            ct2_col = "CytoTRACE2"
            print(f"\n[read] {C.CYTOTRACE_CSV}: CytoTRACE2_Score joined as an "
                  "extra stemness readout, independent of the signature lists")
    else:
        print(f"\n[warn] {C.CYTOTRACE_CSV} not found - CytoTRACE2 is left out of the "
              "quadrant definitions.\n"
              "       The stemness axis then rests entirely on the lab's lists, with no "
              "evidence independent of them.\n"
              "       Run cytotrace2_epi.py (see its header; the package installs with "
              "`pip install cytotrace2-py`).")

    stem_readouts = [s for s in C.STEMNESS_SIGS if s in raw.columns] + ([ct2_col] if ct2_col else [])
    imm_readouts = [s for s in C.IMMUNE_SIGS if s in raw.columns]
    print(f"\nstemness readouts: {', '.join(stem_readouts)}")
    print(f"immune readouts  : {', '.join(imm_readouts)}   (primary: {C.PRIMARY_IMMUNE})")
    assert C.PRIMARY_IMMUNE in raw.columns, f"{C.PRIMARY_IMMUNE} was not scored"

    # ------------------------------------------------------- confounder table
    C.banner("A3 - confounders")
    conf_keys = ["n_genes_by_counts", "pct_counts_mt", "S_score", "G2M_score"]
    rows = []
    for name in raw.columns:
        r = {"readout": name, "axis": C.SIG_AXIS.get(name, "stemness (independent)")}
        for k in conf_keys:
            rho, p = spearmanr(raw[name].values, adata.obs[k].values)
            r[f"rho_{k}"] = rho
            r[f"p_{k}"] = p
        rows.append(r)
    conf = pd.DataFrame(rows).set_index("readout")
    print(conf[[f"rho_{k}" for k in conf_keys]].to_string(float_format="%+.3f"))
    C.write_table(conf, "confounders")

    cc_max = conf[["rho_S_score", "rho_G2M_score"]].abs().max(axis=1)
    print(f"\nstrongest cell-cycle coupling: {cc_max.idxmax()} (|rho| = {cc_max.max():.3f})")
    print(f"strongest depth coupling     : {conf['rho_n_genes_by_counts'].abs().idxmax()} "
          f"(|rho| = {conf['rho_n_genes_by_counts'].abs().max():.3f})")

    # -------------------------------------------------------- standardisation
    C.banner("A4 - standardisation within (cohort, cell_type)")
    z = standardise_within(raw, adata.obs)
    n_strata = adata.obs.groupby(GROUPBY, observed=True).size()
    print(f"{len(n_strata)} strata, smallest {n_strata.min()} cells, largest {n_strata.max():,}")
    print(f"singleton strata (z forced to 0): {int((n_strata == 1).sum())}")

    out = pd.concat([raw.add_prefix("score_"), z], axis=1)
    out.to_csv(C.SCORES_CSV)
    print(f"[write] {C.SCORES_CSV}")

    # --------------------------------------------------------- A5 - quadrants
    C.banner(f"A5 - the stemness x immunogenicity plane "
             f"(high >= q{args.high_q:.2f}, low <= q{args.low_q:.2f})")
    z_imm = z[f"z_{C.PRIMARY_IMMUNE}"]

    quads = {}
    for s in stem_readouts:
        quads[s] = define_quadrant(z[f"z_{s}"], z_imm, args.high_q, args.low_q)

    qdf = pd.DataFrame(quads)
    sizes = qdf.sum().rename("n_cells").to_frame()
    sizes["pct_of_compartment"] = 100 * sizes["n_cells"] / adata.n_obs
    print("\ntarget quadrant (stem-high / immunogenic-low), one definition per stemness readout")
    print(sizes.to_string(float_format="%.2f"))

    # stability across definitions
    stab = pd.DataFrame(index=stem_readouts, columns=stem_readouts, dtype=float)
    for a in stem_readouts:
        for b in stem_readouts:
            stab.loc[a, b] = jaccard(qdf[a], qdf[b])
    print("\nstability of the cell set across definitions (Jaccard)")
    print(stab.to_string(float_format="%.3f"))
    off = stab.where(~np.eye(len(stem_readouts), dtype=bool))
    print(f"median pairwise Jaccard: {np.nanmedian(off.values):.3f}  "
          f"(range {np.nanmin(off.values):.3f} - {np.nanmax(off.values):.3f})")
    C.write_table(stab.round(4), "quadrant_stability")

    n_defs = qdf.sum(axis=1)
    consensus = n_defs >= max(2, int(np.ceil(len(stem_readouts) / 2)))
    print(f"\ncalled by >=1 definition: {int((n_defs >= 1).sum()):,} cells; "
          f"by a majority: {int(consensus.sum()):,}; by all {len(stem_readouts)}: "
          f"{int((n_defs == len(stem_readouts)).sum()):,}")
    votes = n_defs.value_counts().sort_index().rename("n_cells").to_frame()
    votes.index.name = "n_definitions_calling_the_cell"
    C.write_table(votes, "quadrant_vote_distribution")

    # per-patient sizes: a state present in one patient is a patient effect until shown otherwise
    per_pat = qdf.copy()
    per_pat["cohort"] = adata.obs["cohort"].values
    pp = per_pat.groupby("cohort", observed=True).sum()
    pp["n_cells_in_cohort"] = adata.obs["cohort"].value_counts().reindex(pp.index).values
    pp["consensus"] = pd.Series(consensus.values, index=adata.obs["cohort"].values).groupby(level=0).sum()
    pp["pct_consensus"] = 100 * pp["consensus"] / pp["n_cells_in_cohort"]
    print("\nper-patient size of the consensus quadrant "
          "(a state in one patient is a patient effect until shown otherwise)")
    print(pp[["n_cells_in_cohort", "consensus", "pct_consensus"]]
          .sort_values("pct_consensus", ascending=False).to_string(float_format="%.2f"))
    print(f"\npatients with at least one consensus cell: "
          f"{int((pp['consensus'] > 0).sum())} / {len(pp)}")
    C.write_table(pp, "quadrant_per_patient")

    per_ct = pd.DataFrame({"cell_type": adata.obs["cell_type"].values,
                           "consensus": consensus.values}).groupby("cell_type", observed=True).agg(
        n_cells=("consensus", "size"), n_target=("consensus", "sum"))
    per_ct["pct"] = 100 * per_ct["n_target"] / per_ct["n_cells"]
    print("\nconsensus quadrant by cell type")
    print(per_ct.sort_values("pct", ascending=False).to_string(float_format="%.2f"))
    C.write_table(per_ct, "quadrant_per_cell_type")

    # ------------------------------------------- the two named risks (A3 cont.)
    C.banner("A3 - the two named risks")

    # 1. is stem-high just cycling?
    phase = adata.obs["phase"].astype(str).values
    comp = pd.crosstab(pd.Series(phase, name="phase"), consensus.values,
                       normalize="columns") * 100
    comp.columns = ["rest", "target"]
    print("phase composition, target quadrant vs the rest (%)")
    print(comp.to_string(float_format="%.2f"))

    g1 = phase == "G1"
    print(f"\nrecomputing the quadrant inside G1 alone ({g1.sum():,} cells), "
          "i.e. with the cycle held out")
    g1_quads = {}
    for s in stem_readouts:
        zs, zi = z[f"z_{s}"][g1], z_imm[g1]
        g1_quads[s] = define_quadrant(zs, zi, args.high_q, args.low_q)
    g1df = pd.DataFrame(g1_quads)
    g1_consensus = g1df.sum(axis=1) >= max(2, int(np.ceil(len(stem_readouts) / 2)))
    overlap = jaccard(consensus[g1], g1_consensus)
    print(f"consensus quadrant restricted to G1: {int(consensus[g1].sum()):,} cells")
    print(f"consensus quadrant recomputed within G1: {int(g1_consensus.sum()):,} cells")
    print(f"Jaccard between the two: {overlap:.3f}")
    print("A high Jaccard means the state is not an artefact of the cycle; a low one means\n"
          "'stem-high' was largely 'cycling' and the stemness axis does not survive the check.")

    cc_rows = [{"check": "phase_pct_G1_target", "value": float(comp.loc["G1", "target"]) if "G1" in comp.index else np.nan},
               {"check": "phase_pct_G1_rest", "value": float(comp.loc["G1", "rest"]) if "G1" in comp.index else np.nan},
               {"check": "n_target_all_phases", "value": float(consensus.sum())},
               {"check": "n_target_within_G1_recomputed", "value": float(g1_consensus.sum())},
               {"check": "jaccard_target_vs_G1_recomputed", "value": float(overlap)}]

    # 2. is immunogenic-low just shallow?
    imm_low = z_imm <= z_imm.quantile(args.low_q)
    depth = adata.obs["n_genes_by_counts"].values
    u, pu = mannwhitneyu(depth[imm_low.values], depth[~imm_low.values], alternative="two-sided")
    med_lo, med_hi = np.median(depth[imm_low.values]), np.median(depth[~imm_low.values])
    auc_depth = roc_auc_score(imm_low.values, -depth)
    print(f"\nimmunogenic-low ({int(imm_low.sum()):,} cells) vs the rest, n_genes_by_counts:")
    print(f"  median {med_lo:,.0f} vs {med_hi:,.0f}  (Mann-Whitney p = {pu:.3g})")
    print(f"  AUROC of 'shallower' predicting immunogenic-low: {auc_depth:.3f}")
    print("  0.5 means depth does not explain the group; well above it means the immune-evasive\n"
          "  group is largely the low-complexity group and the finding is technical.")
    cc_rows += [{"check": "median_depth_immunogenic_low", "value": float(med_lo)},
                {"check": "median_depth_rest", "value": float(med_hi)},
                {"check": "mannwhitney_p_depth", "value": float(pu)},
                {"check": "auroc_depth_predicts_immunogenic_low", "value": float(auc_depth)}]
    C.write_table(pd.DataFrame(cc_rows).set_index("check"), "confounder_checks", index=True)

    # ----------------------------------------------- A6 - dimensions x signatures
    C.banner("A6 - latent dimensions vs the standardised scores")
    embed = ad.read_h5ad(C.EMBED_H5AD)
    assert (embed.obs_names == adata.obs_names).all(), \
        "the embedding and the all-genes object do not hold the same cells in the same order"
    dims = C.analysis_dimensions(embed)
    n_van = C.n_vanished(embed)
    print(f"{embed.n_vars} dimensions, all {len(dims)} used "
          f"({n_van} flagged vanished in var['vanished'] and NOT pruned)")

    L = pd.DataFrame(np.asarray(embed.X), index=embed.obs_names, columns=embed.var["title"].values)[dims]

    readouts = list(raw.columns)
    rho = pd.DataFrame(index=dims, columns=readouts, dtype=float)
    for r in readouts:
        rr, _ = spearmanr(L.values, z[f"z_{r}"].values)
        rho[r] = rr[:-1, -1] if rr.ndim else np.nan
    # A dimension with no spread left in it has no rank correlation and comes back NaN.
    # 0 is the honest value - no association - and it keeps the row in the table, which is
    # the point of not pruning in the first place.
    rho = rho.fillna(0.0)
    C.write_table(rho.round(4), "dim_signature_spearman")
    print(f"\nstrongest |rho| per readout:")
    for r in readouts:
        print(f"  {r:24s} {rho[r].abs().idxmax():>8s}  rho = {rho.loc[rho[r].abs().idxmax(), r]:+.3f}")

    auroc = pd.Series({d: roc_auc_score(consensus.values, L[d].values) for d in dims},
                      name="auroc_target_vs_rest")
    # `sd == 0` is only reachable on a dead dimension, now that they are kept: 0 there, not inf.
    smd = pd.Series({d: ((L[d][consensus.values].mean() - L[d][~consensus.values].mean())
                         / L[d].std(ddof=0) if L[d].std(ddof=0) > 0 else 0.0)
                     for d in dims}, name="standardised_mean_difference")
    eff = pd.concat([auroc, smd], axis=1)
    eff["abs_auroc_from_half"] = (eff["auroc_target_vs_rest"] - 0.5).abs()
    eff = eff.sort_values("abs_auroc_from_half", ascending=False)
    print("\neffect size of the consensus target quadrant on each dimension (top 10)")
    print(eff.head(10).to_string(float_format="%+.3f"))
    C.write_table(eff, "dim_target_effect_size")

    # The row order every later heatmap uses, Route B included.
    order = pd.Series(dims, name="dimension").to_frame().assign(
        drvi_order=embed.var.set_index("title").loc[dims, "order"].values,
        vanished=embed.var.set_index("title").loc[dims, "vanished"].astype(bool).values)
    C.write_table(order.set_index("dimension"), "dimension_row_order")

    # -------------------------------------------------------------- figures
    C.banner("figures")

    # confounder heatmap
    fig, ax = plt.subplots(figsize=(6.5, 5))
    m = conf[[f"rho_{k}" for k in conf_keys]]
    m.columns = conf_keys
    sns.heatmap(m, cmap="vlag", center=0, vmin=-1, vmax=1, annot=True, fmt="+.2f",
                annot_kws={"size": 7}, cbar_kws={"label": "Spearman rho", "shrink": 0.7}, ax=ax)
    ax.set_title("Route A confounders: raw signature scores vs technical and cycle covariates",
                 fontsize=9)
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right", fontsize=8)
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=8)
    C.savefig("confounder_heatmap", "04_5_cell_first", fig)
    plt.close(fig)

    # the plane, one panel per stemness readout
    ncol = 4
    nrow = int(np.ceil(len(stem_readouts) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 3.8 * nrow), squeeze=False)
    rng = np.random.default_rng(C.SEED)
    idx = rng.choice(adata.n_obs, size=min(20000, adata.n_obs), replace=False)
    for a, s in zip(axes.ravel(), stem_readouts):
        x, y = z[f"z_{s}"].values[idx], z_imm.values[idx]
        tgt = qdf[s].values[idx]
        a.scatter(x[~tgt], y[~tgt], s=1.5, c="0.78", lw=0, rasterized=True)
        a.scatter(x[tgt], y[tgt], s=1.5, c="#C44E52", lw=0, rasterized=True)
        a.axvline(z[f"z_{s}"].quantile(args.high_q), color="k", ls="--", lw=0.8)
        a.axhline(z_imm.quantile(args.low_q), color="k", ls="--", lw=0.8)
        a.set_title(f"{s}\n{int(qdf[s].sum()):,} cells in quadrant", fontsize=9)
        a.set_xlabel(f"z {s} (within cohort x cell_type)", fontsize=8)
        a.set_ylabel(f"z {C.PRIMARY_IMMUNE}", fontsize=8)
        a.tick_params(labelsize=7)
        sns.despine(ax=a)
    for a in axes.ravel()[len(stem_readouts):]:
        a.axis("off")
    fig.suptitle("Route A: the stemness x immunogenicity plane, one definition per stemness readout\n"
                 f"target quadrant in red = stem-high (>=q{args.high_q:.2f}) and "
                 f"immunogenic-low (<=q{args.low_q:.2f}); {len(idx):,} cells shown",
                 fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    C.savefig("stemness_immunogenicity_plane", "04_5_cell_first", fig)
    plt.close(fig)

    # quadrant stability
    fig, ax = plt.subplots(figsize=(6, 5))
    sns.heatmap(stab.astype(float), cmap="rocket_r", vmin=0, vmax=1, annot=True, fmt=".2f",
                annot_kws={"size": 7}, square=True,
                cbar_kws={"label": "Jaccard of the called cell set", "shrink": 0.7}, ax=ax)
    ax.set_title("Stability of the target quadrant across stemness definitions", fontsize=10)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=8)
    C.savefig("quadrant_stability", "04_5_cell_first", fig)
    plt.close(fig)

    # dimensions x signatures
    col_order = [c for c in C.IMMUNE_SIGS if c in rho.columns] + \
                [c for c in C.STEMNESS_SIGS if c in rho.columns] + \
                ([ct2_col] if ct2_col else [])
    fig, ax = plt.subplots(figsize=(1.0 * len(col_order) + 4, 0.24 * len(dims) + 3))
    sns.heatmap(rho[col_order].astype(float), cmap="vlag", center=0, vmin=-0.6, vmax=0.6,
                cbar_kws={"label": "Spearman rho (dimension vs within-stratum z-score)",
                          "shrink": 0.4}, ax=ax)
    ax.set_title("Route A: latent dimensions x signatures\n"
                 f"all {len(dims)} dimensions of {C.RUN_ID}, nothing pruned "
                 f"({n_van} of them flagged vanished in var['vanished'])", fontsize=10)
    ax.axvline(len([c for c in C.IMMUNE_SIGS if c in rho.columns]), color="k", lw=1.5)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=6)
    C.savefig("dim_signature_heatmap", "04_5_cell_first", fig)
    plt.close(fig)

    if skipped:
        print(f"\n[warn] signatures skipped for having under {C.MIN_SIGNATURE_GENES} "
              f"mapped genes: {', '.join(skipped)}")
    print("\ndone.")


if __name__ == "__main__":
    main()
