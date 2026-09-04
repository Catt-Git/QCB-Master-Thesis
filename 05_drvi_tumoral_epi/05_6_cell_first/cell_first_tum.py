#!/usr/bin/env python3
"""05_6, Route A: cell-first. Prior knowledge defines the state, DRVI is the coordinate system.

The unit of analysis is the CELL. The question is whether the cells that prior knowledge calls
stem-like, immune-evasive, or hybrid-EMT occupy a distinct position along any latent dimension.
This is the only route that assigns cells to states - Route B can tell you a dimension is
enriched for stemness genes, and still not tell you which cells, how many, or in which
patients - so the deliverable of the project, an actual cell assignment, comes from here.

Which readouts, and what shape the target region has on the plane, come from the collection
(`utils/sig_collections.py`), not from this file:

  * `--collection scie` crosses each stemness readout with the primary immune one and takes
    the stem-HIGH / immunogenic-LOW corner;
  * `--collection emt` crosses the epithelial score with the mesenchymal score of the same list
    version and takes the corner where BOTH are high, because a partial-EMT cell is one
    co-expressing the two programmes rather than one sitting at either end of the axis. The
    hybrid gene lists are scored and reported but do not define the region - see
    `sig_collections.py` for the robustness argument behind that.

How it fails, and what this script does about it:

  * per-cell scores are noisy and correlate with depth, with the cycle, and with each
    other -> the confounder table, and the named checks below;
  * the stemness lists are embryonic and proliferation-heavy, and this subset is 51%
    `Lumsec-prol` (18,617 of 36,192, against 25% of 04's epithelium), so "stem-high" can just
    mean "cycling" -> the target set is recomputed inside G1 alone and the two cell sets
    compared. This check runs for every collection: an EMT score is as capable of tracking the
    cycle as a stemness one, and on this subset it is the check to read first;
  * immune evasion is defined by ABSENCE of signal, which shallow sequencing mimics
    perfectly -> the depth of the immunogenic-low group is compared against the rest
    (`risks = ("depth", ...)`);
  * a high mesenchymal score is as easily fibroblast ambient RNA or an epithelial-fibroblast
    doublet as it is a transition -> the doublet score of the mesenchymal-high group is
    compared against the rest (`risks = ("ambient", ...)`). The malignant subset NARROWS this
    risk: a cell had to be called aneuploid to be in the object, so an actual fibroblast is
    not among them. Ambient RNA survives - it is contamination, not identity - so the check
    stays;
  * absolute scores are not comparable across patients sequenced at different depths ->
    every readout is standardised within `cohort` before it is used. NOT within
    (cohort, cell_type), which is what 04 does: see GROUPBY below, it is the one deliberate
    departure of this step;
  * it is confirmatory by construction and can only find what was brought in from outside.
    That one has no fix inside Route A; it is why Route B exists.

THE COORDINATE SYSTEM IS AN ARGUMENT, NOT PART OF THE QUESTION. A1 - A5 below - scoring,
standardisation, the target region, the consensus and the confounder checks - are computed
from `shiao_tum.h5ad` alone and never see an embedding; only A6 does. `--embedding` therefore
reads the same Route A off another space, and today DRVI is the only one registered (04 also
registers Harmony and PCA, written by its 04_9 embedding control; that step is not planned
here). The mechanism is kept because it is what states, in code, that the latent space is not
part of the question. Route B and Route C are DRVI-only and are not parameterised.

Scoring is on the UNINTEGRATED, ALL-GENES object. A 150-gene signature reduced to whatever
survived HVG selection is no longer that signature, so `shiao_tum.h5ad` (36,192 x 24,779) is
read here and never `shiao_tum_hvg_2k.h5ad`.

Every latent dimension is correlated, vanished ones included (`PRUNE_VANISHED = False` in
signature_common). Dropping dimensions before the correlations would decide ahead of the
analysis which axes are allowed to mean something; the vanished flag is reported next to
the results instead.

No single definition of the target is primary: the region is defined once per plane - one per
stemness readout for `scie`, one per list version for `emt` - and the stability of the
resulting cell set across those definitions is itself a reported result.

Usage:
    export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
    python cell_first_tum.py                              # the scie collection, the default
    python cell_first_tum.py --collection emt             # the same procedure on the EMT lists
    python cell_first_tum.py --high-q 0.80 --low-q 0.20   # a stricter target region
    python cell_first_tum.py --overwrite                  # re-score instead of reusing the csv
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
import sig_collections as SC  # noqa: E402

# Non-default scoring parameters, all on the record for the Methods.
N_BINS = 25                 # sc.tl.score_genes default is 25; stated because it is not optional

# THE STANDARDISATION STRATA, AND THE ONE DELIBERATE DEPARTURE FROM 04.
#
# 04_5 standardises within (cohort, cell_type). This step uses `cohort` alone, and the reason
# is not that `cell_type` is constant here - that is true, and it would make the second key a
# no-op rather than a mistake. The reason is what the second key would be if it were not
# constant.
#
# In 04, `cell_type` was a LINEAGE: luminal against basal are different cells, and
# standardising within them asks "inside a cell type, does this dimension track stemness?" -
# a sensible question when the groups are different kinds of cell. Inside one malignant
# compartment the available groupings (`cell_type_01_4`, the leiden partition) are STATES of
# the same tumour, and state is the quantity this phase measures. Standardising within a
# state removes the contrast being looked for, by construction.
#
# The concrete danger, spelled out because it is not hypothetical: `Lumsec-prol` means
# PROLIFERATING and is 51% of this subset, and the cell cycle is a named risk of the `scie`
# collection. Standardising within that label would subtract the proliferation axis before
# anything is measured - which would not remove the confounder, it would HIDE it, and it is
# exactly what the G1 recomputation below exists to measure instead.
#
# The pre-CNV label and the leiden partition are therefore reported as COVARIATES of the
# target set (the `quadrant_per_*` tables), never regressed out of it.
GROUPBY = ["cohort"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    SC.add_argument(p)
    C.add_embedding_argument(p)
    C.add_cutoff_arguments(p)
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
    """Z-score each column within the strata of GROUPBY, i.e. within `cohort`.

    Absolute scores are not comparable across patients sequenced at different depths. A
    stratum whose scores have no spread (a single cell, say) would divide by zero; those
    become 0, i.e. 'exactly average for its stratum', which is the only honest value for a
    stratum of one. See GROUPBY for why the second key of 04 is not here.
    """
    g = df.groupby([strata[k].values for k in GROUPBY], observed=True)
    z = g.transform(lambda s: (s - s.mean()) / s.std(ddof=0) if s.std(ddof=0) > 0 else s * 0.0)
    return z.fillna(0.0).add_prefix("z_")


# --------------------------------------------------------------------------- #
# Quadrants
# --------------------------------------------------------------------------- #

# `in_region`, `define_target` and the majority rule live in signature_common: anything that
# re-derives the same consensus cell set has to do it from the same cached scores, and two
# copies of a quadrant rule are two chances to be compared on different cells.


def axis_label(coll, name: str) -> str:
    """The axis a readout sits on, marking the ones that do not come from a gene list.

    CytoTRACE2 is on the stemness axis but is not one of the lists, and that distinction is
    the whole reason it is in the table: it is the only stemness evidence in the stage that
    cannot be circular with the lists it is being compared to.
    """
    if name in coll.extra_readouts:
        return f"{coll.extra_readouts[name]} (independent)"
    return coll.axis_of.get(name, "")


def jaccard(a: pd.Series, b: pd.Series) -> float:
    inter = int((a & b).sum())
    union = int((a | b).sum())
    return inter / union if union else np.nan


def main():
    args = parse_args()
    coll = SC.get(args.collection)
    # Before anything is written: this sets the run id every table and figure name carries.
    emb = C.set_embedding(args.embedding)
    C.banner(f"05_6 - Route A, cell-first: {coll.title}")
    print(f"question  {coll.question}")
    print(f"space     {emb.title} ({emb.run_id}) - {emb.description}")

    # A1 - A5 describe the CELLS: they are computed from the all-genes object and cannot move
    # when the coordinate system does. The reference run owns them; a control run recomputes
    # them identically - which is the point, the target region must not be allowed to follow
    # the embedding - and does not write a second, byte-identical copy under its own run id.
    # Only A6 differs between runs, and that is exactly what is being compared.
    ref = C.EMBEDDINGS[C.DEFAULT_EMBEDDING]

    # A control run's figures would NOT go in this step's folder: this folder is the phase's
    # own run, and a second space's heatmap sitting next to the DRVI one invites the two to be
    # read as one result set. They would go under an embedding-control folder named for the
    # METHOD. With DRVI the only registered space, the first branch is the only one taken.
    fig_step = ("05_6_cell_first" if emb.is_reference
                else f"05_x_embedding_control/{emb.name}")

    def write_shared(df, name, *a, **k):
        """Write an embedding-independent table, or say who owns it and skip."""
        if not emb.is_reference:
            print(f"[skip] {name}: embedding-independent, owned by the {ref.run_id} run")
            return None
        return C.write_table(df, name, *a, **k)

    def savefig_shared(name, step, coll_, fig):
        if not emb.is_reference:
            print(f"[skip] {name}: embedding-independent, owned by the {ref.run_id} run")
            plt.close(fig)
            return None
        path = C.savefig(name, step, coll_, fig)
        plt.close(fig)
        return path

    sets = C.read_gmt(C.gmt_path(coll))
    print(f"{len(sets)} signatures from {C.gmt_path(coll)}")

    # ------------------------------------------------------------------ input
    adata = ad.read_h5ad(C.FULL_H5AD)
    print(f"{adata.n_obs:,} cells x {adata.n_vars:,} genes (all genes, unintegrated)")
    # The object has to be the malignant one: Route A's whole claim is that these scores are
    # read on aneuploid cells. A stale path pointing at 04's object would score 74,441
    # epithelial cells and report them under this phase's caveat.
    assert adata.obs[C.LABEL_KEY].astype(str).nunique() >= 1, "no cell_type column"
    assert (adata.obs["cnv_status"].astype(str) != "not_tested").all(), \
        "not_tested cells in the object: this is not the 05_2 subset"

    # `.X` must be log-normalised: score_genes on raw counts is a different statistic.
    xmax = float(adata.X.max())
    assert xmax < 50, (f".X looks like counts (max {xmax:.1f}), not log-normalised data. "
                       "Route A must score on the scran log1p matrix.")
    print(f".X max {xmax:.3f}, log1p in uns: {'log1p' in adata.uns}  -> log-normalised")

    for k in ("S_score", "G2M_score", "phase", "n_genes_by_counts", "pct_counts_mt"):
        assert k in adata.obs, f"missing obs key: {k}"
    print("cell cycle scores reused from 05_2, not recomputed and not regressed out")

    # ----------------------------------------------------------------- scores
    scores_csv = C.scores_csv(coll)
    if scores_csv.exists() and not args.overwrite:
        print(f"\n[read] {scores_csv} (--overwrite to re-score)")
        allsc = pd.read_csv(scores_csv, index_col=0)
        assert (allsc.index == adata.obs_names).all(), "score csv is not in the object's cell order"
        raw = allsc[[c for c in allsc.columns if not c.startswith("z_")]]
        raw.columns = [c.replace("score_", "") for c in raw.columns]
        # Derived readouts are rebuilt below from the z-scores, so they are dropped here
        # rather than read back: keeping both would let the two versions diverge.
        raw = raw.drop(columns=[d.name for d in coll.derived], errors="ignore")
        skipped = [n for n in sets if n not in raw.columns and n not in coll.extra_readouts]
    else:
        print("\nscoring (sc.tl.score_genes, use_raw=False, ctrl_size=len(signature), "
              f"n_bins={N_BINS}, random_state={C.SEED})")
        raw, skipped = score_signatures(adata, sets)

    # ------------------------------------------------------------ CytoTRACE2
    # Only the scie collection declares it: it is a stemness readout, and joining it into the
    # EMT run would put a stemness axis into a table that has none.
    ct2_col = None
    if "CytoTRACE2" in coll.extra_readouts and C.CYTOTRACE_CSV.exists():
        ct2 = pd.read_csv(C.CYTOTRACE_CSV, index_col=0).reindex(adata.obs_names)
        if "CytoTRACE2_Score" in ct2:
            raw["CytoTRACE2"] = ct2["CytoTRACE2_Score"].values
            ct2_col = "CytoTRACE2"
            print(f"\n[read] {C.CYTOTRACE_CSV}: CytoTRACE2_Score joined as an "
                  "extra stemness readout, independent of the signature lists")
    elif "CytoTRACE2" in coll.extra_readouts:
        print(f"\n[warn] {C.CYTOTRACE_CSV} not found - CytoTRACE2 is left out of the "
              "quadrant definitions.\n"
              "       The stemness axis then rests entirely on the lab's lists, with no "
              "evidence independent of them.\n"
              "       Run 05_5_cytotrace2/cytotrace2_tum.py from the `cytotrace2-py` env (see its header; "
              "the package pins numpy<2 and must not be installed into benchmark-py-r).")

    # -------------------------------------------------------- standardisation
    # Before the confounders, because the derived readouts are contrasts of z-scores and have
    # to exist by the time anything is correlated against a covariate.
    C.banner(f"A3 - standardisation within ({', '.join(GROUPBY)})")
    z = standardise_within(raw, adata.obs)
    n_strata = adata.obs.groupby(GROUPBY, observed=True).size()
    print(f"{len(n_strata)} strata, smallest {n_strata.min()} cells, largest {n_strata.max():,}")
    print(f"singleton strata (z forced to 0): {int((n_strata == 1).sum())}")

    # ----------------------------------------------------- derived readouts
    # z(plus) - z(minus), on the STANDARDISED scores: the two halves have to be on the same
    # scale before they are differenced, or the contrast just tracks whichever list is longer.
    for d in coll.derived:
        missing = [n for n in (d.plus, d.minus) if f"z_{n}" not in z.columns]
        if missing:
            print(f"  {d.name:24s} SKIPPED, missing {', '.join(missing)}")
            continue
        v = z[f"z_{d.plus}"] - z[f"z_{d.minus}"]
        z[f"z_{d.name}"] = v
        raw[d.name] = v          # so it is correlated and reported like any other readout
        print(f"  {d.name:24s} {d.description}")

    readout_axes = {n: axis_label(coll, n) for n in raw.columns}
    for ax_name in coll.axes:
        on_axis = [n for n in raw.columns if coll.axis_of.get(n) == ax_name]
        if on_axis:
            print(f"\n{ax_name:12s} readouts: {', '.join(on_axis)}")

    planes = coll.planes(coll, list(raw.columns))
    assert planes, (f"the {coll.name} collection defines no plane on the readouts that were "
                    "scored - nothing can be called a target region")
    for pl in planes:
        for k in (pl.x, pl.y):
            assert k in raw.columns, f"plane '{pl.label}' needs readout {k}, which was not scored"
    print(f"\n{len(planes)} definitions of the target region ({coll.target_label}):")
    for pl in planes:
        print(f"  {pl.label:24s} {pl.x} {pl.x_rule} x {pl.y} {pl.y_rule}")

    # ------------------------------------------------------- confounder table
    C.banner("A4 - confounders")
    conf_keys = ["n_genes_by_counts", "pct_counts_mt", "S_score", "G2M_score"]
    rows = []
    for name in raw.columns:
        r = {"readout": name, "axis": readout_axes[name]}
        for k in conf_keys:
            rho, p = spearmanr(raw[name].values, adata.obs[k].values)
            r[f"rho_{k}"] = rho
            r[f"p_{k}"] = p
        rows.append(r)
    conf = pd.DataFrame(rows).set_index("readout")
    print(conf[[f"rho_{k}" for k in conf_keys]].to_string(float_format="%+.3f"))
    write_shared(conf, "confounders", coll)

    cc_max = conf[["rho_S_score", "rho_G2M_score"]].abs().max(axis=1)
    print(f"\nstrongest cell-cycle coupling: {cc_max.idxmax()} (|rho| = {cc_max.max():.3f})")
    print(f"strongest depth coupling     : {conf['rho_n_genes_by_counts'].abs().idxmax()} "
          f"(|rho| = {conf['rho_n_genes_by_counts'].abs().max():.3f})")

    out = pd.concat([raw.add_prefix("score_"), z], axis=1)
    if emb.is_reference or not scores_csv.exists():
        out.to_csv(scores_csv)
        print(f"\n[write] {scores_csv}")
    else:
        # Named for the object and not for the run (see signature_common), so a control run
        # would be rewriting the reference run's file with identical content.
        print(f"\n[skip] {scores_csv} is embedding-independent and already there")

    # ----------------------------------------------------- A5 - target region
    C.banner(f"A5 - the target region: {coll.target_label}\n"
             f"(high >= q{args.high_q:.2f}, low <= q{args.low_q:.2f}, "
             f"mid = q{args.mid_lo_q:.2f} - q{args.mid_hi_q:.2f})")

    labels = [pl.label for pl in planes]
    cut = C.Cutoffs.from_args(args)
    quads = {pl.label: C.define_target(z, pl, cut) for pl in planes}

    qdf = pd.DataFrame(quads)
    sizes = qdf.sum().rename("n_cells").to_frame()
    sizes["pct_of_compartment"] = 100 * sizes["n_cells"] / adata.n_obs
    print(f"\ntarget region ({coll.target_label}), one definition per plane")
    print(sizes.to_string(float_format="%.2f"))

    # stability across definitions
    stab = pd.DataFrame(index=labels, columns=labels, dtype=float)
    for a in labels:
        for b in labels:
            stab.loc[a, b] = jaccard(qdf[a], qdf[b])
    print("\nstability of the cell set across definitions (Jaccard)")
    print(stab.to_string(float_format="%.3f"))
    off = stab.where(~np.eye(len(labels), dtype=bool))
    print(f"median pairwise Jaccard: {np.nanmedian(off.values):.3f}  "
          f"(range {np.nanmin(off.values):.3f} - {np.nanmax(off.values):.3f})")
    write_shared(stab.round(4), "quadrant_stability", coll)

    n_defs, consensus = C.consensus_vote(qdf)
    print(f"\ncalled by >=1 definition: {int((n_defs >= 1).sum()):,} cells; "
          f"by a majority: {int(consensus.sum()):,}; by all {len(labels)}: "
          f"{int((n_defs == len(labels)).sum()):,}")
    votes = n_defs.value_counts().sort_index().rename("n_cells").to_frame()
    votes.index.name = "n_definitions_calling_the_cell"
    write_shared(votes, "quadrant_vote_distribution", coll)

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
    write_shared(pp, "quadrant_per_patient", coll)

    # 04 breaks the target set down by `cell_type`. That column is the constant `malignant`
    # here, so the breakdown is by every grouping this subset actually has - the leiden
    # partition of 05_2 first, the pre-CNV CellTypist label second - plus `phase`, which is
    # not a grouping but the covariate the cycle risk is read on.
    #
    # These are COVARIATES of the target set, reported. They are deliberately not in GROUPBY:
    # standardising within a state is what would remove the contrast being measured. A target
    # set that turns out to be one leiden cluster, or all `Lumsec-prol`, is a finding to state
    # - possibly a negative one - not something to correct away here.
    for key in C.grouping_keys(adata.obs) + ["phase"]:
        per_g = pd.DataFrame({key: adata.obs[key].astype(str).values,
                              "consensus": consensus.values}).groupby(key, observed=True).agg(
            n_cells=("consensus", "size"), n_target=("consensus", "sum"))
        per_g["pct"] = 100 * per_g["n_target"] / per_g["n_cells"]
        print(f"\nconsensus quadrant by {key}")
        print(per_g.sort_values("pct", ascending=False).to_string(float_format="%.2f"))
        write_shared(per_g, f"quadrant_per_{key}", coll)

    # ------------------------------------------- the named risks (A4 cont.)
    C.banner(f"A4 - the named risks of this collection: {', '.join(coll.risks)}")

    # 1. is the target set just the cycling one? Runs for every collection.
    phase = adata.obs["phase"].astype(str).values
    comp = pd.crosstab(pd.Series(phase, name="phase"), consensus.values,
                       normalize="columns") * 100
    comp.columns = ["rest", "target"]
    print("phase composition, target quadrant vs the rest (%)")
    print(comp.to_string(float_format="%.2f"))

    g1 = phase == "G1"
    print(f"\nrecomputing the target region inside G1 alone ({g1.sum():,} cells), "
          "i.e. with the cycle held out")
    # The quantile cutoffs are recomputed WITHIN G1, not carried over: the point of the check
    # is what the definition would have called had the cycling cells never been there.
    z_g1 = z[g1]
    g1_quads = {pl.label: C.define_target(z_g1, pl, cut) for pl in planes}
    g1df = pd.DataFrame(g1_quads)
    _, g1_consensus = C.consensus_vote(g1df)
    overlap = jaccard(consensus[g1], g1_consensus)
    print(f"consensus target restricted to G1: {int(consensus[g1].sum()):,} cells")
    print(f"consensus target recomputed within G1: {int(g1_consensus.sum()):,} cells")
    print(f"Jaccard between the two: {overlap:.3f}")
    print("A high Jaccard means the state is not an artefact of the cycle; a low one means the\n"
          "target was largely 'cycling' and the readout does not survive the check.")

    cc_rows = [{"check": "phase_pct_G1_target", "value": float(comp.loc["G1", "target"]) if "G1" in comp.index else np.nan},
               {"check": "phase_pct_G1_rest", "value": float(comp.loc["G1", "rest"]) if "G1" in comp.index else np.nan},
               {"check": "n_target_all_phases", "value": float(consensus.sum())},
               {"check": "n_target_within_G1_recomputed", "value": float(g1_consensus.sum())},
               {"check": "jaccard_target_vs_G1_recomputed", "value": float(overlap)}]

    # 2. is the low end of the primary axis just shallow sequencing? (`risks` contains "depth")
    if "depth" in coll.risks and coll.depth_risk_readout:
        r = coll.depth_risk_readout
        low_grp = z[f"z_{r}"] <= z[f"z_{r}"].quantile(args.low_q)
        depth = adata.obs["n_genes_by_counts"].values
        u, pu = mannwhitneyu(depth[low_grp.values], depth[~low_grp.values], alternative="two-sided")
        med_lo, med_hi = np.median(depth[low_grp.values]), np.median(depth[~low_grp.values])
        auc_depth = roc_auc_score(low_grp.values, -depth)
        print(f"\n{r}-low ({int(low_grp.sum()):,} cells) vs the rest, n_genes_by_counts:")
        print(f"  median {med_lo:,.0f} vs {med_hi:,.0f}  (Mann-Whitney p = {pu:.3g})")
        print(f"  AUROC of 'shallower' predicting {r}-low: {auc_depth:.3f}")
        print("  0.5 means depth does not explain the group; well above it means the evasive\n"
              "  group is largely the low-complexity group and the finding is technical.")
        cc_rows += [{"check": "median_depth_immunogenic_low", "value": float(med_lo)},
                    {"check": "median_depth_rest", "value": float(med_hi)},
                    {"check": "mannwhitney_p_depth", "value": float(pu)},
                    {"check": "auroc_depth_predicts_immunogenic_low", "value": float(auc_depth)}]

    # 3. is the high end of the mesenchymal axis just ambient RNA or a doublet?
    #    (`risks` contains "ambient")
    #
    # This subset was defined by a CNV call, so an actual fibroblast is not in it - but
    # subsetting removes cells, not the fibroblast transcripts that leaked into the droplets of
    # the cells that remain. VIM / FN1 / SPARC / ACTA2 high is therefore still the expected
    # signature of contamination as much as of a transition, and without this check an EMT
    # result cannot be told apart from a soup result. What the malignant subset buys is the
    # other half of the risk: a high-mesenchymal cell here cannot simply BE a fibroblast.
    # `doublet_score` comes from 01_2 via 05_2 (Scrublet) and is not recomputed here.
    if "ambient" in coll.risks and coll.ambient_risk_axis:
        on_axis = [n for n in raw.columns if coll.axis_of.get(n) == coll.ambient_risk_axis]
        dbl = adata.obs["doublet_score"].astype(float).values
        print(f"\n{coll.ambient_risk_axis}-high vs the rest, doublet_score (Scrublet, via 05_2):")
        for r in on_axis:
            high = (z[f"z_{r}"] >= z[f"z_{r}"].quantile(args.high_q)).values
            auc = roc_auc_score(high, dbl)
            rho_d, p_d = spearmanr(raw[r].values, dbl)
            med_hi_d, med_lo_d = np.median(dbl[high]), np.median(dbl[~high])
            print(f"  {r:24s} AUROC {auc:.3f}   median {med_hi_d:.4f} vs {med_lo_d:.4f}   "
                  f"rho(score, doublet_score) {rho_d:+.3f}")
            cc_rows += [{"check": f"auroc_doublet_predicts_{r}_high", "value": float(auc)},
                        {"check": f"rho_doublet_score_{r}", "value": float(rho_d)}]
        # The predicted-doublet flag is a harder call than the score and is reported next to it.
        if "predicted_doublet" in adata.obs:
            pred = adata.obs["predicted_doublet"].astype(bool).values
            print(f"  predicted doublets in the consensus target: "
                  f"{int(pred[consensus.values].sum()):,} / {int(consensus.sum()):,} "
                  f"({100 * pred[consensus.values].mean():.2f}%) vs "
                  f"{100 * pred[~consensus.values].mean():.2f}% in the rest")
            cc_rows += [{"check": "pct_predicted_doublet_target",
                         "value": float(100 * pred[consensus.values].mean())},
                        {"check": "pct_predicted_doublet_rest",
                         "value": float(100 * pred[~consensus.values].mean())}]
        print("  AUROC near 0.5 means the mesenchymal signal is not the soup; well above it\n"
              "  means the EMT readout is measuring contamination and nothing here is a state.")

    write_shared(pd.DataFrame(cc_rows).set_index("check"), "confounder_checks", coll, index=True)

    # ----------------------------------------------- A6 - dimensions x signatures
    C.banner(f"A6 - the dimensions of {emb.title} vs the standardised scores")
    if not C.EMBED_H5AD.exists():
        raise SystemExit(
            f"{C.EMBED_H5AD} not found.\n"
            + ("Run 05_3_drvi_run/run_drvi_tum.py (or drvi_tum.ipynb) first."
               if emb.is_reference else
               "That embedding has no writer in this phase; see utils/signature_common.py."))
    embed = ad.read_h5ad(C.EMBED_H5AD)
    assert (embed.obs_names == adata.obs_names).all(), \
        "the embedding and the all-genes object do not hold the same cells in the same order"
    dims = C.analysis_dimensions(embed)
    n_van = C.n_vanished(embed)
    print(f"{embed.n_vars} dimensions, all {len(dims)} used")
    if emb.is_reference:
        print(f"({n_van} flagged vanished in var['vanished'] and NOT pruned)")
    else:
        # `vanished` is a DRVI notion. The column is there and all False; saying so is
        # clearer than a count of zero that could be read as a filter having run.
        print("no vanished dimensions: this space has no such notion, every dimension is used")

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
    C.write_table(rho.round(4), "dim_signature_spearman", coll)
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
    C.write_table(eff, "dim_target_effect_size", coll)

    # The row order every later heatmap uses, Route B included.
    # The order column is named for the embedding, so `drvi_order` is what the DRVI run
    # writes, unchanged, and a control run cannot pass its own ordering off as DRVI's.
    order = pd.Series(dims, name="dimension").to_frame().assign(
        **{f"{emb.name}_order": embed.var.set_index("title").loc[dims, "order"].values},
        vanished=embed.var.set_index("title").loc[dims, "vanished"].astype(bool).values)
    C.write_table(order.set_index("dimension"), "dimension_row_order", coll)

    # -------------------------------------------------------------- figures
    C.banner("figures")

    # confounder heatmap
    fig, ax = plt.subplots(figsize=(6.5, 5))
    m = conf[[f"rho_{k}" for k in conf_keys]]
    m.columns = conf_keys
    sns.heatmap(m, cmap="vlag", center=0, vmin=-1, vmax=1, annot=True, fmt="+.2f",
                annot_kws={"size": 7}, cbar_kws={"label": "Spearman rho", "shrink": 0.7}, ax=ax)
    ax.set_title(f"Route A confounders, {coll.title}\n"
                 "raw signature scores vs technical and cycle covariates", fontsize=9)
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right", fontsize=8)
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=8)
    savefig_shared("confounder_heatmap", "05_6_cell_first", coll, fig)

    # the plane, one panel per definition of the target region
    def cut_lines(a, v, rule, vertical: bool):
        """Draw the cutoff(s) of one rule. 'mid' has two, which is what makes it visible as a
        BAND rather than a corner - the reader has to be able to see that the EMT target is
        not an extreme of the axis."""
        draw = a.axvline if vertical else a.axhline
        qs = {"high": [args.high_q], "low": [args.low_q],
              "mid": [args.mid_lo_q, args.mid_hi_q]}[rule]
        for q in qs:
            draw(v.quantile(q), color="k", ls="--", lw=0.8)

    ncol = min(4, len(planes))
    nrow = int(np.ceil(len(planes) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 3.8 * nrow), squeeze=False)
    rng = np.random.default_rng(C.SEED)
    idx = rng.choice(adata.n_obs, size=min(20000, adata.n_obs), replace=False)
    for a, pl in zip(axes.ravel(), planes):
        zx, zy = z[f"z_{pl.x}"], z[f"z_{pl.y}"]
        x, y = zx.values[idx], zy.values[idx]
        tgt = qdf[pl.label].values[idx]
        a.scatter(x[~tgt], y[~tgt], s=1.5, c="0.78", lw=0, rasterized=True)
        a.scatter(x[tgt], y[tgt], s=1.5, c="#C44E52", lw=0, rasterized=True)
        cut_lines(a, zx, pl.x_rule, vertical=True)
        cut_lines(a, zy, pl.y_rule, vertical=False)
        a.set_title(f"{pl.label}\n{int(qdf[pl.label].sum()):,} cells in the target region", fontsize=9)
        a.set_xlabel(f"z {pl.x} ({pl.x_rule}) (within {' x '.join(GROUPBY)})", fontsize=8)
        a.set_ylabel(f"z {pl.y} ({pl.y_rule})", fontsize=8)
        a.tick_params(labelsize=7)
        sns.despine(ax=a)
    for a in axes.ravel()[len(planes):]:
        a.axis("off")
    fig.suptitle(f"Route A, {coll.title}: one definition per plane\n"
                 f"target region in red = {coll.target_label}; {len(idx):,} cells shown",
                 fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    savefig_shared(coll.plane_figure, "05_6_cell_first", coll, fig)

    # quadrant stability
    fig, ax = plt.subplots(figsize=(6, 5))
    sns.heatmap(stab.astype(float), cmap="rocket_r", vmin=0, vmax=1, annot=True, fmt=".2f",
                annot_kws={"size": 7}, square=True,
                cbar_kws={"label": "Jaccard of the called cell set", "shrink": 0.7}, ax=ax)
    ax.set_title(f"Stability of the target region across definitions, {coll.title}", fontsize=10)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=8)
    savefig_shared("quadrant_stability", "05_6_cell_first", coll, fig)

    # dimensions x signatures
    #
    # One row per DIMENSION, not per dimension-direction: Route A correlates the latent
    # coordinate itself, which has no direction of its own. The direction is carried by the
    # SIGN of rho, exactly as 05_8 reads it when it joins the two routes - rho > 0 on DR 7 is
    # a statement about `DR 7+`, rho < 0 about `DR 7-`. Route B's heatmap says so on its
    # colorbar and in its title, and without the same wording here the reader has a signed
    # colour scale with nothing telling them what the sign means. Hence both lines below.
    col_order = coll.order(list(rho.columns))
    fig, ax = plt.subplots(figsize=(1.0 * len(col_order) + 4, 0.24 * len(dims) + 3))
    sns.heatmap(rho[col_order].astype(float), cmap="vlag", center=0, vmin=-0.6, vmax=0.6,
                cbar_kws={"label": "Spearman rho (dimension vs within-stratum z-score)\n"
                                   "sign = direction: rho > 0 is DR n+, rho < 0 is DR n-",
                          "shrink": 0.4}, ax=ax)
    pruning = (f"nothing pruned ({n_van} of them flagged vanished in var['vanished'])"
               if emb.is_reference else "nothing pruned")
    ax.set_title(f"Route A, {coll.title}: {emb.title} dimensions x signatures\n"
                 f"all {len(dims)} dimensions of {emb.run_id}, {pruning}\n"
                 "rows are dimensions, not dimension-directions: the sign of rho IS the "
                 "direction, as 05_8 joins them", fontsize=10)
    for pos in coll.block_edges(col_order):
        ax.axvline(pos, color="k", lw=1.5)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=6)
    C.savefig("dim_signature_heatmap", fig_step, coll, fig)
    plt.close(fig)

    if skipped:
        print(f"\n[warn] signatures skipped for having under {C.MIN_SIGNATURE_GENES} "
              f"mapped genes: {', '.join(skipped)}")
    print("\ndone.")


if __name__ == "__main__":
    main()
