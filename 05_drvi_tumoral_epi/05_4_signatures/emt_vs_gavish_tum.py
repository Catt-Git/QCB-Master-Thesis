#!/usr/bin/env python3
"""05_4: do the collaborator's EMT lists and Gavish's EMT metaprograms describe the same thing?

EMT enters this phase twice, from two unrelated directions, and the two never meet in any
other step:

  * the `emt` collection - the collaborator's nine lists, three axes (epithelial, hybrid,
    mesenchymal) x three curations (B primary, A its validated core, C Tomas' mouse list).
    A HYPOTHESIS: it names a state and 05_6 calls cells in it.
  * the `emt` axis of the `gavish` collection - MP12 - MP15, "EMT-I" to "EMT-IV", four of the
    41 recurrent pan-cancer metaprograms of Gavish et al. 2023 (MP16, the glioma mesenchymal
    programme, joins them when the run is `--all-metaprograms`). A VOCABULARY: it names
    dimensions and calls no cells.

That is exactly why the agreement is worth measuring. An EMT axis that shows up on the
collaborator's lists AND on Gavish's metaprograms is an axis that does not depend on whose
EMT list was used - the argument `sig_collections.py` makes for carrying the metaprograms at
all. But 05_8 currently reaches that conclusion through the CELLS (both collections
correlate against the same latent dimensions), and a shared correlation is only evidence of a
shared axis if the two are not simply the same genes twice. This step measures the genes.

WHAT COMES OUT, AND HOW TO READ IT. The headline number is small and easy to misread, so the
figure puts three reference bands next to it rather than the cross-block alone:

  * lab x lab, same axis      - the same biological axis, three curations of it. The CEILING:
                                what "these two lists mean the same thing" looks like here.
  * lab x lab, different axis - epithelial against mesenchymal, inside one curation. The
                                FLOOR: two lists that are deliberately about different things.
  * Gavish x Gavish, EMT only - MP12 - MP15 against each other. The band that matters most:
                                these four are all called EMT by their own authors, so
                                whatever overlap they show each other is what "two EMT gene
                                sets" is actually worth as a number.
  * lab x Gavish              - the question.

A Jaccard of 0.08 between two 30-to-50-gene lists reads as "no agreement" until it is put
beside the 0.00 - 0.10 that Gavish's own four EMT metaprograms show EACH OTHER. The
metaprograms are NMF programmes clustered across tumours, not curated marker panels: four
EMT programmes that barely share genes is their result, not a defect of ours. So the
comparison to make is against that band, not against 1.0, and the hypergeometric column of
the pairs table is the second half of the same point - six shared genes out of 25,133 is a
tiny Jaccard and an enormous enrichment at the same time, and those two facts are not in
conflict.

Genes are the MAPPED ones, as in `build_signatures_tum.py`: the collection as actually used
on this object, so the numbers here and the numbers in the 05_4 Jaccard matrix are on the
same footing. The enrichment universe is the object's full gene axis for the same reason.

Usage:
    export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
    python emt_vs_gavish_tum.py                        # the 4 TNBC-relevant EMT metaprograms
    python emt_vs_gavish_tum.py --all-metaprograms     # adds MP16, the glioma mesenchymal one
    CELL_SET=epi HVG_SET=nomt N_LATENT=64 python emt_vs_gavish_tum.py
"""

from __future__ import annotations

import argparse
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
from scipy.stats import hypergeom

UTILS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils")
sys.path.insert(0, UTILS_DIR)
import signature_common as C  # noqa: E402
import sig_collections as SC  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--all-metaprograms", action="store_true",
                   help="take the EMT axis of the full 41-metaprogram registry rather than of "
                        "the TNBC-relevant subset, which adds MP16 (glioma mesenchymal)")
    return p.parse_args()


def jaccard(a: set, b: set) -> float:
    return len(a & b) / len(a | b) if (a | b) else 0.0


def main():
    args = parse_args()
    gav = SC.GAVISH_VARIANTS[bool(args.all_metaprograms)]

    C.banner("05_4 - the collaborator's EMT lists against Gavish's EMT metaprograms")

    # ------------------------------------------------------------------ inputs
    # The EMT axis of the metaprogram collection, whichever variant was asked for. Selecting
    # on the axis rather than on the MP numbers is what makes `--all-metaprograms` pick up
    # MP16 without this script knowing which number it is.
    mp_names = [s.name for s in gav.signatures if gav.axis_of[s.name] == "emt"]
    if not mp_names:
        sys.exit(f"[stop] the {gav.name} collection declares no signature on the 'emt' axis")

    lab_raw = C.load_signatures(SC.EMT)
    gav_raw = {n: g for n, g in C.load_signatures(gav).items() if n in mp_names}

    print(f"\ncollaborator lists : {len(lab_raw)}  ({', '.join(lab_raw)})")
    print(f"Gavish EMT metaprograms : {len(gav_raw)}  ({', '.join(gav_raw)})")

    # ------------------------------------------------------------- the universe
    # Same object and same mapping as `build_signatures_tum.py`. A symbol that is not on this
    # gene axis was NOT MEASURED here, and counting it as a match between two lists would
    # claim an agreement on a gene neither collection could see.
    adata = ad.read_h5ad(C.FULL_H5AD, backed="r")
    universe = set(adata.var_names)
    n_universe = adata.n_vars
    adata.file.close()
    print(f"\nuniverse: {n_universe:,} genes of {C.FULL_H5AD.name}")

    lab = {k: set(g) & universe for k, g in lab_raw.items()}
    gavm = {k: set(g) & universe for k, g in gav_raw.items()}
    for k, g in {**lab, **gavm}.items():
        src = lab_raw.get(k, gav_raw.get(k))
        if len(g) < len(src):
            print(f"  {k:24s} {len(g):3d}/{len(src):3d} symbols mapped")

    # ------------------------------------------------------------ the matrices
    # Figure order: the collaborator's lists grouped by AXIS and, inside an axis, A / B / C -
    # `EMT.figure_key`, the same reordering the 05_4 Jaccard matrix uses, so the two figures
    # can be read side by side. The metaprograms follow, in MP number order.
    lab_order = SC.EMT.order(list(lab), for_figure=True)
    gav_order = [n for n in gav.order(list(gavm), for_figure=True)]
    order = lab_order + gav_order
    sets = {**lab, **gavm}

    J = pd.DataFrame(np.eye(len(order)), index=order, columns=order)
    O = pd.DataFrame(0, index=order, columns=order, dtype=int)
    for a, b in itertools.combinations(order, 2):
        J.loc[a, b] = J.loc[b, a] = jaccard(sets[a], sets[b])
        O.loc[a, b] = O.loc[b, a] = len(sets[a] & sets[b])
    for a in order:
        O.loc[a, a] = len(sets[a])

    C.write_table(J.round(4), "jaccard_emt_vs_gavish", SC.EMT)
    C.write_table(O, "shared_genes_emt_vs_gavish", SC.EMT)

    # ------------------------------------------------- the four reference bands
    # Each band is a list of (a, b) pairs; the figure and the printout read the same dict, so
    # a band cannot be drawn from one definition and quoted from another.
    bands: dict[str, list[tuple[str, str]]] = {
        "lab x lab\nsame axis": [(a, b) for a, b in itertools.combinations(lab_order, 2)
                                 if SC.EMT.axis_of[a] == SC.EMT.axis_of[b]],
        "lab x lab\ndifferent axis": [(a, b) for a, b in itertools.combinations(lab_order, 2)
                                      if SC.EMT.axis_of[a] != SC.EMT.axis_of[b]],
        "Gavish x Gavish\nEMT MPs only": list(itertools.combinations(gav_order, 2)),
        "lab x Gavish": [(a, b) for a in lab_order for b in gav_order],
    }

    print("\nJaccard by band (the cross-block is the question, the other three are its scale)")
    summary = []
    for label, pairs in bands.items():
        v = np.array([J.loc[a, b] for a, b in pairs])
        one = label.replace("\n", ", ")
        summary.append({"band": one, "n_pairs": len(v), "median": float(np.median(v)),
                        "min": float(v.min()), "max": float(v.max())})
        print(f"  {one:34s} n={len(v):3d}   median {np.median(v):.3f}   "
              f"range {v.min():.3f} - {v.max():.3f}")
    C.write_table(pd.DataFrame(summary).set_index("band").round(4),
                  "jaccard_bands_emt_vs_gavish", SC.EMT)

    cross = np.array([J.loc[a, b] for a, b in bands["lab x Gavish"]])
    within_g = np.array([J.loc[a, b] for a, b in bands["Gavish x Gavish\nEMT MPs only"]])
    print(f"\nthe comparison that matters: lab x Gavish median {np.median(cross):.3f} against "
          f"the {np.median(within_g):.3f} that Gavish's own EMT metaprograms show each other.")
    print("Read the cross-block against THAT number, not against 1.0.")

    # ------------------------------------------------------- the pairs, in full
    # One row per cross pair, with the genes spelled out: a Jaccard is a summary, and the
    # question "which genes actually agree?" is the one a reader asks next. The
    # hypergeometric column is the other half of the reading - see the header.
    rows = []
    for a, b in bands["lab x Gavish"]:
        A, B = sets[a], sets[b]
        inter = sorted(A & B)
        k, n_a, n_b = len(inter), len(A), len(B)
        # P(X >= k) for k genes shared between an m-gene and an n-gene draw from n_universe.
        p = float(hypergeom.sf(k - 1, n_universe, n_a, n_b)) if k else 1.0
        rows.append({"lab_list": a, "metaprogram": b, "n_lab": n_a, "n_mp": n_b,
                     "n_shared": k, "jaccard": J.loc[a, b],
                     "overlap_coefficient": k / min(n_a, n_b),
                     "expected_by_chance": n_a * n_b / n_universe,
                     "hypergeom_p": p, "shared_genes": ";".join(inter)})
    pairs = pd.DataFrame(rows)
    # BH over the cross-block only. The bands above are descriptive and are not tested: they
    # exist to scale the cross-block, and putting them in the same FDR would mix a reference
    # measurement into the family being corrected.
    o = pairs["hypergeom_p"].rank(method="first").astype(int)
    pairs["hypergeom_fdr"] = np.minimum.accumulate(
        (pairs["hypergeom_p"].sort_values(ascending=False).values
         * len(pairs) / np.arange(len(pairs), 0, -1)).clip(max=1.0)
    )[::-1][o - 1]
    pairs = pairs.sort_values("jaccard", ascending=False)
    C.write_table(pairs.round(6), "emt_vs_gavish_pairs", SC.EMT, index=False)

    print("\nthe ten strongest cross pairs:")
    print(pairs.head(10)[["lab_list", "metaprogram", "n_shared", "jaccard",
                          "expected_by_chance", "hypergeom_fdr", "shared_genes"]]
          .to_string(index=False, float_format="%.3g"))

    empty = pairs[pairs["n_shared"] == 0]
    if len(empty):
        print(f"\n{len(empty)} of the {len(pairs)} cross pairs share NO gene at all:")
        print("  " + ", ".join(f"{r.lab_list} x {r.metaprogram}" for r in empty.itertuples()))

    # ------------------------------------------------------------------ figure
    # THREE PANELS, because one is not enough and the reason is the colour scale. The full
    # matrix has to be drawn on a single scale to be honest, and on that scale the block this
    # step is about - a Jaccard of 0.00 - 0.10 - is indistinguishable from zero next to the
    # 0.76 of two curations of the same list. So the cross-block is drawn a second time on
    # its OWN scale, which is a zoom and is labelled as one, and the bands beneath it are
    # what says whether the zoomed numbers are large or small.
    n = len(order)
    fig = plt.figure(figsize=(16.8, 8.8))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.12, 1.0], height_ratios=[1.3, 1.0],
                          wspace=0.28, hspace=0.62)
    ax = fig.add_subplot(gs[:, 0])
    axz = fig.add_subplot(gs[0, 1])
    axb = fig.add_subplot(gs[1, 1])

    # ---- panel A: the whole matrix, one scale
    off = J.where(~np.eye(n, dtype=bool))
    labelled = J.copy()
    labelled.index = labelled.columns = [f"{m} ({len(sets[m])})" for m in order]
    sns.heatmap(labelled, cmap="rocket_r", vmin=0, vmax=float(np.nanmax(off.values)),
                annot=C.annotate_cells(n, n), fmt=".2f", annot_kws={"size": 7}, square=True,
                cbar_kws={"label": "Jaccard index", "shrink": 0.42, "pad": 0.02}, ax=ax)

    # Block edges: the three axes of the collaborator's collection, then the metaprograms as
    # a fourth block. The boundary between lab and Gavish is drawn heavier than the others -
    # it is the only one that separates two different KINDS of gene set.
    lab_edges = SC.EMT.block_edges(lab_order)
    for pos in lab_edges:
        ax.axhline(pos, color="black", lw=1.4)
        ax.axvline(pos, color="black", lw=1.4)
    ax.axhline(len(lab_order), color="black", lw=3)
    ax.axvline(len(lab_order), color="black", lw=3)
    blocks = list(zip([0] + lab_edges, lab_edges + [len(lab_order)]))
    for start_i, stop_i in blocks:
        ax.text((start_i + stop_i) / 2, -0.18, SC.EMT.axis_of[lab_order[start_i]], ha="center",
                va="bottom", fontsize=8, weight="bold")
    ax.text(len(lab_order) + len(gav_order) / 2, -0.18, "Gavish EMT MPs", ha="center",
            va="bottom", fontsize=8, weight="bold")
    ax.set_title("A. Gene overlap of the two EMT vocabularies, one scale\n"
                 "(mapped genes; the heavy line separates the collaborator's curated\n"
                 "lists from Gavish's metaprograms)", fontsize=10, pad=30)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=8)

    # ---- panel B: the cross-block alone, on its own scale
    # Rows are the collaborator's lists and columns the metaprograms, so the panel reads in
    # the direction the question is asked: for each of my lists, which metaprogram is it.
    sub = J.loc[lab_order, gav_order]
    # Two numbers per cell: the Jaccard above, the raw count of shared genes below it. The
    # count is not decoration - a Jaccard of 0.08 is 6 genes out of 81 between two lists of
    # this size and 2 genes out of 25 between two much shorter ones, and only one of those is
    # something to go and look at. `square=False` here and not in panel A: four columns drawn
    # square against nine rows leaves cells narrower than the numbers written in them.
    cell = np.array([[f"{sub.loc[a, b]:.2f}\n{O.loc[a, b]:d}" for b in gav_order]
                     for a in lab_order])
    sns.heatmap(sub, cmap="rocket_r", vmin=0, vmax=float(sub.values.max()),
                annot=cell, fmt="", annot_kws={"size": 7.5}, square=False, linewidths=0.8,
                linecolor="white",
                cbar_kws={"label": "Jaccard index", "shrink": 0.75, "pad": 0.02}, ax=axz)
    axz.set_yticks(np.arange(len(lab_order)) + 0.5)
    axz.set_yticklabels([f"{a} ({len(sets[a])})" for a in lab_order], rotation=0, fontsize=8)
    axz.set_xticks(np.arange(len(gav_order)) + 0.5)
    axz.set_xticklabels([f"{b} ({len(sets[b])})" for b in gav_order], rotation=20,
                        ha="right", fontsize=8)
    for pos in SC.EMT.block_edges(lab_order):
        axz.axhline(pos, color="black", lw=1.4)
    axz.set_title("B. The same cross-block, rescaled to itself "
                  f"(max {sub.values.max():.2f}): a ZOOM, not a second\n"
                  "measurement, and the only panel where these numbers are legible.\n"
                  "Each cell: Jaccard above, number of shared genes below",
                  fontsize=9.5, pad=12)

    # ---- panel C: every pair, by band
    # One dot per pair, the bands stacked. Grey for the three reference bands, one accent for
    # the band being asked about: the colour marks WHICH SERIES IS THE QUESTION, so it must
    # not be spent on telling the three references apart from each other.
    rng = np.random.default_rng(C.SEED)
    accent = sns.color_palette("rocket_r", 6)[4]
    labels = list(bands)
    for i, label in enumerate(labels):
        v = np.array([J.loc[a, b] for a, b in bands[label]])
        is_q = label == "lab x Gavish"
        axb.scatter(v, np.full(len(v), i) + rng.uniform(-0.15, 0.15, len(v)),
                    s=24, color=accent if is_q else "0.62",
                    edgecolor="white", linewidth=0.6, zorder=3, alpha=0.95)
        axb.plot([np.median(v)] * 2, [i - 0.3, i + 0.3],
                 color=accent if is_q else "0.3", lw=2.4, zorder=4, solid_capstyle="butt")
        axb.text(1.005, i, f"n={len(v)}", transform=axb.get_yaxis_transform(),
                 va="center", ha="left", fontsize=7, color="0.35")
    # The band the cross-block has to be read against, carried across the panel so the two
    # medians are compared by eye rather than by arithmetic.
    axb.axvline(np.median(within_g), color="0.3", ls=":", lw=1.2, zorder=1)
    axb.set_yticks(range(len(labels)))
    axb.set_yticklabels([l.replace("\n", " ") for l in labels], fontsize=8)
    axb.set_ylim(len(labels) - 0.5, -0.5)
    axb.set_xlim(-0.02, float(np.nanmax(off.values)) * 1.06)
    axb.set_xlabel("Jaccard index", fontsize=9)
    axb.set_title("C. Every pair, by band. Thick tick = median; the dotted line is the\n"
                  "median Gavish's own four EMT metaprograms reach against EACH OTHER,\n"
                  "which is the scale panel B has to be read on", fontsize=9.5, pad=12)
    axb.grid(axis="x", color="0.9", lw=0.7)
    axb.set_axisbelow(True)
    sns.despine(ax=axb, left=True)
    axb.tick_params(axis="y", length=0)

    C.savefig("jaccard_emt_vs_gavish", "05_4_signatures", SC.EMT, fig, caveat=False)
    plt.close(fig)

    print("\ndone.")


if __name__ == "__main__":
    main()
