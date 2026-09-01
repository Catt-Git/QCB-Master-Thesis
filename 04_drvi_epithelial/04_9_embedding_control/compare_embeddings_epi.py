#!/usr/bin/env python3
"""04_9: Route A read off two coordinate systems, side by side.

Once `cell_first_epi.py` has been run with `--embedding drvi` and `--embedding harmony`, the
same cells, the same scores and the same target region exist in two spaces and the only thing
that differs is the axes. This step reads them back and asks what changes.

WHAT IS UNDER TEST, AND WHAT IS NOT
-----------------------------------
**DRVI is the hypothesis; Harmony is the reference level.** The phase chose DRVI for a latent
space whose axes are meant to be individually interpretable, and the README's justification
for that choice is about Route B ("without it any of the higher-scoring phase-02 methods
would have done"). Nothing in it is about Route A, which is the route that assigns the cells.
This step supplies the missing half.

The asymmetry is deliberate and has to be read that way: **Harmony has never claimed
axis-level interpretability**, it claims batch correction with biological signal preserved,
and it was benchmarked on exactly that in phase 02. A high `effective_n_dims` for Harmony is
therefore NOT a defect of Harmony, and this step must never be quoted as ranking the two
methods - phase 02 is where they are ranked, on what they both promise. What is being tested
here is whether DRVI delivers on ITS promise for these readouts, with Harmony as the level a
space that makes no such promise happens to reach.

TWO QUESTIONS, AND ONLY THE FIRST IS FAIR TO BOTH
-------------------------------------------------
1. **Does the space carry the state at all?** Multivariate, and method-agnostic: how well
   the consensus target region is predicted from ALL the dimensions of the space, cross-
   validated by patient. Neither method is asked for anything it does not offer.
2. **Is the state on ONE axis?** `max |rho|` per readout and the concentration of the
   association across dimensions. This is DRVI's own claim; Harmony is the floor.

The two together make the useful statement possible - "Harmony carries the same information
but on no single axis, DRVI puts it on one" - which neither number can make alone. If instead
the multivariate AUROCs match AND the concentrations match, Route A is not an argument for
DRVI and the whole weight of the choice sits on Route B. That is a result, not a failure.

READING `effective_n_dims`: IT NEEDS THE SPACE'S OWN REDUNDANCY NEXT TO IT
-------------------------------------------------------------------------
`effective_n_dims` mixes two things - how aligned a readout is with a single axis, and how
correlated the axes are with each other, because correlated axes SHARE an association. The
two spaces differ structurally on the second: a PCA is orthogonal by construction and its
`effective_rank_of_space` is exactly its dimension count, while a latent space need not be
and DRVI is not (64 dimensions behaving as ~31 independent ones on this object). That bias
runs AGAINST DRVI, not against Harmony, and `effective_rank_of_space` is reported in the same
table so the number is never read on its own.

One more asymmetry, small and stated rather than corrected: the maximum of |rho| over 64
draws is slightly inflated against the same maximum over 50. It cannot be removed without
retraining DRVI at 50, which would change the run the whole phase is built on.

THE CELLS ARE THE SAME CELLS. The consensus target region is embedding-independent, so it is
re-derived here from the per-cell scores 04_5 cached, with `signature_common`'s own quadrant
functions rather than a second copy of the rule, and checked against the vote distribution
04_5 reported. If that check fails the comparison stops: two spaces compared on two different
cell sets is not a comparison.

Outputs (`*` is `_<collection>_epi_embeddings`):

    tables/<coll>/embedding_comparison*.csv        one row per (space, readout)
    tables/<coll>/embedding_target_effect*.csv     the target region on each space, with the
                                                   multivariate AUROC and the space's own rank
    figures/04_9_embedding_control/<coll>/
        dim_signature_heatmap_side_by_side*.png    the heatmaps, one colour scale
        association_concentration*.png             sorted |rho| profile per readout
        max_abs_rho*.png                           the headline bars
        information_vs_alignment*.png              question 1 against question 2

Usage:
    export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
    cd 04_drvi_epithelial/04_9_embedding_control
    python compare_embeddings_epi.py                                  # drvi vs harmony
    python compare_embeddings_epi.py --collection emt
    python compare_embeddings_epi.py --embeddings drvi harmony pca    # with the null arm
    python compare_embeddings_epi.py --no-multivariate                # skip the CV fits

The four quantile flags are the same as 04_5's and must be given the same values: this step
re-derives the target region, so a region defined differently there and here would be two
different cell sets. The defaults match, and the vote check catches the rest.
"""

from __future__ import annotations

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import anndata as ad             # noqa: E402
import numpy as np               # noqa: E402
import pandas as pd              # noqa: E402
import seaborn as sns            # noqa: E402
from sklearn.linear_model import LogisticRegression      # noqa: E402
from sklearn.metrics import roc_auc_score                # noqa: E402
from sklearn.model_selection import GroupKFold           # noqa: E402
from sklearn.pipeline import make_pipeline               # noqa: E402
from sklearn.preprocessing import StandardScaler         # noqa: E402

UTILS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils")
sys.path.insert(0, UTILS_DIR)
import signature_common as C     # noqa: E402
import sig_collections as SC     # noqa: E402

# The run id every cross-run output carries, in place of a single embedding's. It is not a
# run: it says "this file is about several of them", which is why it is spelled out here.
COMPARISON_RUN_ID = "epi_embeddings"

# Grouped by patient, not shuffled: the question a multivariate readout has to answer is
# whether the state is recoverable in a patient the fit has never seen. A random split would
# put cells of the same patient on both sides and report the patient effect as a success.
CV_FOLDS = 5

PALETTE = {"drvi": "#4C72B0", "harmony": "#C44E52", "pca": "#8C8C8C"}

# Said on every figure that touches the concentration question, so the framing cannot be
# lost between the script and the slide.
FRAMING = ("DRVI is the hypothesis under test, Harmony the reference level: it never claimed "
           "axis-level interpretability,\nso a low concentration for Harmony is not a defect "
           "of Harmony. The two methods are ranked in phase 02, not here.")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    SC.add_argument(p)
    p.add_argument("--embeddings", nargs="+", default=["drvi", "harmony"],
                   choices=sorted(C.EMBEDDINGS),
                   help="the spaces to compare, reference first [default: drvi harmony]")
    p.add_argument("--rho-bar", type=float, default=0.20,
                   help="the Route A bar of 04_7, used for the counts [default: 0.20]")
    C.add_cutoff_arguments(p)
    p.add_argument("--no-multivariate", action="store_true",
                   help="skip the cross-validated multivariate AUROC (the only slow part)")
    return p.parse_args()


def load_tables(emb, coll):
    """The two embedding-dependent tables of one Route A run."""
    try:
        rho = C.read_table("dim_signature_spearman", coll, run_id=emb.run_id)
        eff = C.read_table("dim_target_effect_size", coll, run_id=emb.run_id)
    except FileNotFoundError as e:
        raise SystemExit(
            f"missing Route A output for --embedding {emb.name} ({emb.run_id}):\n  {e}\n"
            f"Run it first:  python ../04_5_cell_first/cell_first_epi.py "
            f"--embedding {emb.name} --collection {coll.name}") from None
    return rho.astype(float), eff.astype(float)


def rebuild_target(coll, cut) -> tuple[pd.Series, pd.DataFrame]:
    """The consensus target region, re-derived from the per-cell scores 04_5 cached.

    Embedding-independent by construction, which is the entire point: the two spaces are
    asked about the same cells. Re-derived rather than recomputed - the quadrant functions
    are `signature_common`'s, the same objects 04_5 calls - and then checked against the vote
    distribution 04_5 wrote, so a silent divergence is impossible rather than unlikely.
    """
    path = C.scores_csv(coll)
    if not path.exists():
        raise SystemExit(f"{path} not found - run 04_5 on the reference embedding first")
    allsc = pd.read_csv(path, index_col=0)
    z = allsc[[c for c in allsc.columns if c.startswith("z_")]]
    readouts = [c[len("score_"):] for c in allsc.columns if c.startswith("score_")]

    planes = coll.planes(coll, readouts)
    qdf = pd.DataFrame({pl.label: C.define_target(z, pl, cut) for pl in planes})
    n_defs, consensus = C.consensus_vote(qdf)

    got = n_defs.value_counts().sort_index()
    try:
        want = C.read_table("quadrant_vote_distribution", coll,
                            run_id=C.EMBEDDINGS[C.DEFAULT_EMBEDDING].run_id)["n_cells"]
    except FileNotFoundError:
        print("[warn] no quadrant_vote_distribution table to check the re-derived region "
              "against; proceeding on the re-derived one")
        return consensus, qdf
    want.index = want.index.astype(int)
    if not got.reindex(want.index).fillna(0).astype(int).equals(want.astype(int)):
        raise SystemExit(
            "the re-derived target region does not match the one 04_5 reported:\n"
            f"  here  {dict(got)}\n  04_5  {dict(want)}\n"
            "Give this step the same --high-q/--low-q/--mid-lo-q/--mid-hi-q 04_5 was run "
            "with, or re-run 04_5.")
    print(f"target region re-derived and checked against 04_5: {int(consensus.sum()):,} cells "
          f"in the consensus of {len(planes)} definitions")
    return consensus, qdf


def concentration(v: np.ndarray) -> tuple[float, float]:
    """Inverse Simpson on the shares of sum(v^2), and the share of the largest.

    `v` is one readout's rho across the dimensions of a space. Squaring first is what makes
    this a statement about the association rather than about its sign, and the shares are
    scale-free, so a space whose associations are all weak is not thereby called concentrated.
    """
    w = np.square(np.asarray(v, dtype=float))
    tot = w.sum()
    if tot <= 0:
        return float("nan"), float("nan")
    w = w / tot
    return float(1.0 / np.square(w).sum()), float(w.max())


def multivariate_auroc(X, y, groups, seed=C.SEED) -> tuple[float, float, int]:
    """Cross-validated AUROC of the target region from the WHOLE space, grouped by patient.

    The fair half of the comparison: it asks whether the space carries the state, not whether
    it isolates it on an axis, so it asks Harmony for nothing Harmony does not offer.
    Standardisation is inside the fold, so no fold sees the test cells' scale. L2 logistic
    regression at the default C - a linear readout, deliberately: a non-linear one would
    measure the classifier's capacity as much as the space's geometry, and the question is
    about the geometry.
    """
    n_splits = min(CV_FOLDS, len(np.unique(groups)))
    cv = GroupKFold(n_splits=n_splits)
    aucs = []
    for tr, te in cv.split(X, y, groups):
        if y[tr].sum() == 0 or y[te].sum() == 0:
            continue
        model = make_pipeline(
            StandardScaler(),
            LogisticRegression(max_iter=2000, random_state=seed))
        model.fit(X[tr], y[tr])
        aucs.append(roc_auc_score(y[te], model.predict_proba(X[te])[:, 1]))
    return float(np.mean(aucs)), float(np.std(aucs)), len(aucs)


def main():
    args = parse_args()
    coll = SC.get(args.collection)
    cut = C.Cutoffs.from_args(args)
    embs = [C.get_embedding(n) for n in args.embeddings]
    C.banner(f"04_9 - Route A across {len(embs)} coordinate systems: {coll.title}")
    print("DRVI is the hypothesis under test; Harmony is the reference level, not a rival:\n"
          "it never claimed axis-level interpretability. Phase 02 is where the methods are\n"
          "ranked, on what they both promise.\n")
    for e in embs:
        print(f"  {e.name:8s} {e.run_id:16s} {e.n_dims:3d} dims   {e.description}")

    loaded = {e.name: load_tables(e, coll) for e in embs}
    consensus, _ = rebuild_target(coll, cut)

    readout_sets = {n: list(r.columns) for n, (r, _) in loaded.items()}
    first = readout_sets[embs[0].name]
    for n, cols in readout_sets.items():
        if cols != first:
            print(f"[warn] {n} has a different readout set than {embs[0].name}: "
                  f"only in {n}: {sorted(set(cols) - set(first))}; "
                  f"only in {embs[0].name}: {sorted(set(first) - set(cols))}")
    readouts = coll.order([c for c in first if all(c in cols for cols in readout_sets.values())])
    print(f"\n{len(readouts)} readouts compared: {', '.join(readouts)}")

    # --------------------------------------------- what each space is like in itself
    C.banner("the spaces themselves: how many of their dimensions are independent")
    space_stats = {}
    for e in embs:
        embed = ad.read_h5ad(e.embed_h5ad)
        assert (embed.obs_names == consensus.index).all(), \
            f"{e.name}: the embedding is not in the cell order of the score table"
        X = np.asarray(embed.X, dtype=np.float64)
        pr, mean_r = C.effective_rank(X)
        space_stats[e.name] = {
            "n_dimensions": X.shape[1],
            "effective_rank_of_space": pr,
            "mean_abs_dim_correlation": mean_r,
            "_X": X,
            "_groups": embed.obs["cohort"].astype(str).values,
        }
        print(f"  {e.title:20s} {X.shape[1]:3d} dimensions, effective rank {pr:5.1f} "
              f"(mean |r| between dimensions {mean_r:.3f})")
    print("\nA PCA scores exactly its dimension count: its axes are orthogonal by "
          "construction.\nA latent space need not, and a space whose dimensions are "
          "redundant spreads any readout\nover more of them - which is why "
          "effective_n_dims below is never read without this column.")

    # ------------------------------------------------------- the comparison table
    rows = []
    for e in embs:
        rho, _ = loaded[e.name]
        for r in readouts:
            v = rho[r]
            best = v.abs().idxmax()
            n_above = int((v.abs() >= args.rho_bar).sum())
            eff_n, top_share = concentration(v.values)
            rows.append({
                "embedding": e.name,
                "space": e.title,
                "readout": r,
                "axis": coll.axis_of.get(r, coll.extra_readouts.get(r, "")),
                "n_dimensions": len(v),
                "max_abs_rho": float(v.abs().max()),
                "rho_at_best": float(v.loc[best]),
                "best_dimension": best,
                f"n_dims_abs_rho_ge_{args.rho_bar:g}": n_above,
                f"pct_dims_abs_rho_ge_{args.rho_bar:g}": 100.0 * n_above / len(v),
                "effective_n_dims": eff_n,
                "share_of_rho2_top_dimension": top_share,
                # Repeated on every row on purpose: effective_n_dims must not travel to a
                # slide without the redundancy of the space it was measured in.
                "effective_rank_of_space": space_stats[e.name]["effective_rank_of_space"],
                "mean_abs_dim_correlation": space_stats[e.name]["mean_abs_dim_correlation"],
            })
    comp = pd.DataFrame(rows)
    C.write_table(comp.set_index(["embedding", "readout"]), "embedding_comparison", coll,
                  run_id=COMPARISON_RUN_ID)

    wide = comp.pivot(index="readout", columns="embedding", values="max_abs_rho").loc[readouts]
    print("\nmax |rho| per readout  (question 2, and the bar 04_7 applies)")
    print(wide.to_string(float_format="%.3f"))
    ref = embs[0].name
    for e in embs[1:]:
        d = wide[e.name] - wide[ref]
        print(f"\n{e.name} - {ref}: median {d.median():+.3f}, "
              f"{int((d > 0).sum())}/{len(d)} readouts stronger on {e.name}")

    conc = comp.pivot(index="readout", columns="embedding",
                      values="effective_n_dims").loc[readouts]
    print("\neffective number of dimensions the association is spread over")
    print(conc.to_string(float_format="%.2f"))
    for e in embs:
        st = space_stats[e.name]
        print(f"  {e.name:8s} median {conc[e.name].median():5.2f} of {st['n_dimensions']} "
              f"dimensions, effective rank {st['effective_rank_of_space']:.1f}")

    # --------------------------------------- the target region: one axis vs the space
    C.banner(f"the consensus target region ({coll.target_label}) on each space")
    trows = []
    for e in embs:
        _, eff = loaded[e.name]
        a = eff["auroc_target_vs_rest"]
        dev = (a - 0.5).abs()
        eff_n, top_share = concentration(dev.values)
        row = {
            "embedding": e.name, "space": e.title,
            "n_dimensions": space_stats[e.name]["n_dimensions"],
            "effective_rank_of_space": space_stats[e.name]["effective_rank_of_space"],
            "mean_abs_dim_correlation": space_stats[e.name]["mean_abs_dim_correlation"],
            "best_dimension": dev.idxmax(),
            # ORIENTED, and it has to be: a latent dimension has no privileged sign - the
            # phase says so everywhere else, "the sign of rho IS the direction" - so a
            # dimension separating the target at AUROC 0.315 separates it exactly as well as
            # one at 0.685, in the other direction. Comparing the raw value across spaces
            # would report a strong anti-correlated axis as worse than chance; the signed
            # value is kept next to it because the direction is still worth reading.
            "auroc_best_single_dimension": float(0.5 + dev.max()),
            "auroc_best_single_dimension_signed": float(a.loc[dev.idxmax()]),
            "abs_auroc_from_half": float(dev.max()),
            "max_abs_smd": float(eff["standardised_mean_difference"].abs().max()),
            "n_dims_auroc_dev_ge_0.10": int((dev >= 0.10).sum()),
            "effective_n_dims": eff_n,
            "share_top_dimension": top_share,
        }
        if not args.no_multivariate:
            print(f"  {e.title}: fitting {CV_FOLDS}-fold grouped CV on {space_stats[e.name]['n_dimensions']} "
                  "dimensions ...", flush=True)
            m, sd, k = multivariate_auroc(space_stats[e.name]["_X"],
                                          consensus.values.astype(int),
                                          space_stats[e.name]["_groups"])
            row.update(auroc_multivariate_cv_mean=m, auroc_multivariate_cv_sd=sd,
                       cv_folds_used=k)
        trows.append(row)
    tgt = pd.DataFrame(trows).set_index("embedding")
    print()
    print(tgt.drop(columns=["space"]).to_string(float_format="%.3f"))
    C.write_table(tgt, "embedding_target_effect", coll, run_id=COMPARISON_RUN_ID)

    if not args.no_multivariate:
        print("\nquestion 1 (fair to both): does the space carry the state?  -> "
              "auroc_multivariate_cv_mean")
        print("question 2 (DRVI's own claim): is it on one axis?             -> "
              "auroc_best_single_dimension")
        print("The GAP between them is the state the space represents but does not isolate.\n"
          "The single-dimension figure is ORIENTED (0.5 + |AUROC - 0.5|): a latent axis has\n"
          "no privileged sign, so an anti-correlated dimension separates the target just as\n"
          "well. The signed value is in the table next to it.")
        for e in embs:
            r = tgt.loc[e.name]
            print(f"  {e.title:20s} multivariate {r['auroc_multivariate_cv_mean']:.3f} "
                  f"vs best single dimension {r['auroc_best_single_dimension']:.3f} "
                  f"(oriented; signed {r['auroc_best_single_dimension_signed']:.3f})  "
                  f"(gap {r['auroc_multivariate_cv_mean'] - r['auroc_best_single_dimension']:+.3f})")

    # -------------------------------------------------------------------- figures
    C.banner("figures")

    vmax = 0.6
    fig, axes = plt.subplots(1, len(embs),
                             figsize=(len(embs) * (0.85 * len(readouts) + 3.2), 12))
    axes = np.atleast_1d(axes)
    for a, e in zip(axes, embs):
        rho, _ = loaded[e.name]
        m = rho[readouts]
        m = m.loc[m.abs().max(axis=1).sort_values(ascending=False).index]
        sns.heatmap(m.astype(float), cmap="vlag", center=0, vmin=-vmax, vmax=vmax, ax=a,
                    cbar=(a is axes[-1]),
                    cbar_kws={"label": "Spearman rho (dimension vs within-stratum z-score)",
                              "shrink": 0.35})
        a.set_title(f"{e.title}\n{e.n_dims} dimensions, sorted by max |rho|", fontsize=10)
        plt.setp(a.get_xticklabels(), rotation=45, ha="right", fontsize=7)
        plt.setp(a.get_yticklabels(), fontsize=4)
        for pos in coll.block_edges(readouts):
            a.axvline(pos, color="k", lw=1.2)
    fig.suptitle(f"Route A, {coll.title}: the same scores read off {len(embs)} spaces\n"
                 "same cells, same signatures, same target region - only the axes differ\n"
                 + FRAMING, fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    C.savefig("dim_signature_heatmap_side_by_side", "04_9_embedding_control", coll, fig,
              run_id=COMPARISON_RUN_ID)
    plt.close(fig)

    ncol = 4
    nrow = int(np.ceil(len(readouts) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.4 * ncol, 2.9 * nrow), squeeze=False,
                             sharex=True, sharey=True)
    for a, r in zip(axes.ravel(), readouts):
        for e in embs:
            rho, _ = loaded[e.name]
            v = np.sort(rho[r].abs().values)[::-1]
            a.plot(np.arange(1, len(v) + 1), v, lw=1.4, label=e.title,
                   color=PALETTE.get(e.name))
        a.axhline(args.rho_bar, color="k", ls="--", lw=0.8)
        a.set_title(r, fontsize=8)
        a.set_xlim(0.5, 20.5)
        a.tick_params(labelsize=7)
        sns.despine(ax=a)
    for a in axes.ravel()[len(readouts):]:
        a.axis("off")
    # The x label goes on the lowest USED panel of each column: the last row is short
    # whenever the readouts do not fill the grid, and its empty cells are turned off above.
    for j in range(ncol):
        used = [i for i in range(nrow) if i * ncol + j < len(readouts)]
        if used:
            axes[used[-1], j].set_xlabel("dimension, ranked by |rho|", fontsize=8)
            axes[used[-1], j].tick_params(labelbottom=True)
    for a in axes[:, 0]:
        a.set_ylabel("|Spearman rho|", fontsize=8)
    axes.ravel()[0].legend(fontsize=7, frameon=False)
    ranks = "   ".join(f"{e.title}: effective rank {space_stats[e.name]['effective_rank_of_space']:.0f}"
                       f"/{space_stats[e.name]['n_dimensions']}" for e in embs)
    fig.suptitle(f"Route A, {coll.title}: how concentrated the association is\n"
                 f"top 20 dimensions of each space; dashed line = the 04_7 bar "
                 f"|rho| >= {args.rho_bar:g}\n{ranks}\n" + FRAMING, fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    C.savefig("association_concentration", "04_9_embedding_control", coll, fig,
              run_id=COMPARISON_RUN_ID)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(0.9 * len(readouts) + 3, 4.2))
    x = np.arange(len(readouts))
    w = 0.8 / len(embs)
    for i, e in enumerate(embs):
        ax.bar(x + i * w - 0.4 + w / 2, wide[e.name].values, width=w, label=e.title,
               color=PALETTE.get(e.name))
    ax.axhline(args.rho_bar, color="k", ls="--", lw=0.9)
    ax.text(len(readouts) - 0.5, args.rho_bar, f" 04_7 bar {args.rho_bar:g}", va="bottom",
            ha="right", fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels(readouts, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("max |Spearman rho| over the dimensions", fontsize=9)
    ax.set_title(f"Route A, {coll.title}: the strongest axis each space offers per readout",
                 fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    sns.despine(ax=ax)
    C.savefig("max_abs_rho", "04_9_embedding_control", coll, fig, run_id=COMPARISON_RUN_ID)
    plt.close(fig)

    # The figure the whole control is for: the fair question against the DRVI-specific one.
    if not args.no_multivariate:
        fig, ax = plt.subplots(figsize=(1.9 * len(embs) + 3.2, 4.4))
        x = np.arange(len(embs))
        w = 0.36
        mv = tgt["auroc_multivariate_cv_mean"].reindex([e.name for e in embs])
        sdv = tgt["auroc_multivariate_cv_sd"].reindex([e.name for e in embs])
        sd1 = tgt["auroc_best_single_dimension"].reindex([e.name for e in embs])
        ax.bar(x - w / 2, mv.values, width=w, yerr=sdv.values, capsize=3, color="#55A868",
               label="whole space, patient-grouped CV  (does it carry the state?)")
        ax.bar(x + w / 2, sd1.values, width=w, color="#DD8452",
               label="best single dimension, oriented  (is it on one axis?)")
        for i, e in enumerate(embs):
            gap = mv.iloc[i] - sd1.iloc[i]
            ax.annotate(f"gap {gap:+.3f}", (i, max(mv.iloc[i], sd1.iloc[i])),
                        textcoords="offset points", xytext=(0, 6), ha="center", fontsize=8)
        ax.axhline(0.5, color="k", ls="--", lw=0.9)
        ax.text(len(embs) - 0.5, 0.5, " chance", va="bottom", ha="right", fontsize=7)
        ax.set_xticks(x)
        ax.set_xticklabels([e.title for e in embs], fontsize=9)
        ax.set_ylabel(f"AUROC, {coll.target_label}", fontsize=9)
        ax.set_ylim(0.45, 1.0)
        ax.set_title(f"Route A, {coll.title}: information the space carries vs information "
                     "it isolates\nthe gap is the state the space represents but does not "
                     "put on an axis", fontsize=10)
        ax.legend(fontsize=8, frameon=False, loc="upper left")
        sns.despine(ax=ax)
        C.savefig("information_vs_alignment", "04_9_embedding_control", coll, fig,
                  run_id=COMPARISON_RUN_ID)
        plt.close(fig)

    print("\ndone.")


if __name__ == "__main__":
    main()
