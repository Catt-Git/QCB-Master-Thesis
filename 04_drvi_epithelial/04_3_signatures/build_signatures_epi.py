#!/usr/bin/env python3
"""04_3: ingest one collection's signature files into a .gmt, and characterise the collection.

Reads the plain-text files of the requested collection from `$DATA_DIR/signatures/` (one gene
symbol per line) and writes a single `.gmt` whose description field carries the provenance
string, so the identical collection feeds both routes and doubles as the Appendix table.

Which files, on which axes, is declared in `utils/sig_collections.py`: `--collection scie` is
the ten stemness/immunogenicity lists, `--collection emt` the nine EMT lists. The step is
the same either way, and so is every check below.

Two tables come out of it, both of which have to be read before any result of this stage
is believed:

  * coverage  - how much of each signature is actually measured in this object. These
                lists date from 2007-2012 and carry deprecated symbols; a gene that does
                not map is NOT MEASURED, which is not the same as not expressed. Below
                MIN_MAPPED_FRACTION the script stops rather than scoring a signature that
                is no longer the signature it is named after.
  * jaccard   - the pairwise overlap. The immune four are largely nested and the embryonic
                stemness lists overlap heavily, so they are not independent tests and the
                FDR of Route B must not be presented as if they were. This matrix is how
                the reader sees how much of the apparent agreement between two readouts is
                just shared genes.

The mapping is reported against two universes: all 26,371 genes of the epithelial object,
which is what Route A scores on, and the 2,000 HVGs DRVI was trained on, which is the
Route B. The second column is a warning, not a filter - see the README.

Usage:
    export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
    python build_signatures_epi.py                        # the scie collection, the default
    python build_signatures_epi.py --collection emt       # the same, on the EMT lists
    python build_signatures_epi.py --allow-low-coverage   # report, do not stop
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

# signature_common lives in the phase's utils/, as in 02_2_integration.
UTILS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils")
sys.path.insert(0, UTILS_DIR)
import signature_common as C  # noqa: E402
import sig_collections as SC  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    SC.add_argument(p)
    p.add_argument("--allow-low-coverage", action="store_true",
                   help="report signatures below the coverage floor instead of stopping")
    return p.parse_args()


def main():
    args = parse_args()
    coll = SC.get(args.collection)

    C.banner(f"04_3 - signature collection: {coll.title}")
    print(f"question    {coll.question}")
    print(f"signatures  {C.SIG_DIR}")
    print(f"object      {C.FULL_H5AD}")

    raw = C.load_signatures(coll)
    for name, genes in raw.items():
        role = "primary" if name in coll.primary_names() else "robustness"
        print(f"  {name:24s} {len(genes):5d} unique symbols  [{coll.axis_of[name]}, {role}]")

    # ---------------------------------------------------------------- universes
    adata = ad.read_h5ad(C.FULL_H5AD, backed="r")
    all_genes = set(adata.var_names)
    n_all = adata.n_vars
    adata.file.close()

    hvg = ad.read_h5ad(C.HVG_H5AD, backed="r")
    hvg_genes = set(hvg.var_names)
    n_hvg = hvg.n_vars
    hvg.file.close()

    print(f"\nuniverses: {n_all:,} genes (Route A) | {n_hvg:,} HVGs (Route B ORA background)")

    # ---------------------------------------------------------------- coverage
    rows = []
    mapped = {}
    for name, genes in raw.items():
        m = [g for g in genes if g in all_genes]
        mapped[name] = m
        in_hvg = [g for g in m if g in hvg_genes]
        rows.append({
            "signature": name,
            "axis": coll.axis_of[name],
            "primary": name in coll.primary_names(),
            "provenance": coll.provenance[name],
            "n_genes": len(genes),
            "n_mapped": len(m),
            "mapped_fraction": len(m) / len(genes),
            "n_in_hvg_background": len(in_hvg),
            "hvg_fraction_of_mapped": len(in_hvg) / len(m) if m else np.nan,
        })

    coverage = pd.DataFrame(rows).set_index("signature")
    coverage = coverage.loc[coll.names]

    print("\ncoverage")
    print(coverage[["axis", "primary", "n_genes", "n_mapped", "mapped_fraction",
                    "n_in_hvg_background"]].to_string(float_format="%.3f"))

    # ------------------------------------------------------------- stop checks
    low = coverage.index[coverage["mapped_fraction"] < C.MIN_MAPPED_FRACTION].tolist()
    tiny = coverage.index[coverage["n_mapped"] < C.MIN_SIGNATURE_GENES].tolist()

    if tiny:
        print(f"\n[warn] below {C.MIN_SIGNATURE_GENES} mapped genes, will be skipped "
              f"by Route A: {', '.join(tiny)}")
    if low:
        msg = (f"below the {C.MIN_MAPPED_FRACTION:.0%} mapping floor: {', '.join(low)}. "
               "Low coverage means NOT MEASURED, never not expressed - these lists carry "
               "deprecated symbols and need remapping before they mean anything.")
        if args.allow_low_coverage:
            print(f"\n[warn] {msg}")
        else:
            sys.exit(f"\n[STOP] {msg}\nRe-run with --allow-low-coverage to proceed anyway.")

    # Route B is powered by the HVG intersection, not the mapped one.
    weak_b = coverage.index[coverage["n_in_hvg_background"] < C.MIN_SIGNATURE_GENES].tolist()
    if weak_b:
        print(f"\n[warn] under {C.MIN_SIGNATURE_GENES} genes inside the 2,000-HVG ORA "
              f"background, so effectively untestable in Route B: {', '.join(weak_b)}")

    # ------------------------------------------------------------------- .gmt
    # The .gmt carries the MAPPED genes: it is the collection as actually used, so the
    # Appendix table and the tested sets cannot drift apart.
    C.write_gmt(mapped, C.gmt_path(coll), coll.provenance)
    print(f"\n[write] {C.gmt_path(coll)}")

    C.write_table(coverage, "coverage", coll)

    # ---------------------------------------------------------------- jaccard
    names = list(mapped)
    J = pd.DataFrame(np.eye(len(names)), index=names, columns=names)
    O = pd.DataFrame(0, index=names, columns=names, dtype=int)
    for i, a in enumerate(names):
        A = set(mapped[a])
        for b in names[i + 1:]:
            B = set(mapped[b])
            inter = len(A & B)
            j = inter / len(A | B) if (A | B) else 0.0
            J.loc[a, b] = J.loc[b, a] = j
            O.loc[a, b] = O.loc[b, a] = inter
        O.loc[a, a] = len(A)

    C.write_table(J.round(4), "jaccard", coll)
    C.write_table(O, "shared_genes", coll)

    off = J.where(~np.eye(len(names), dtype=bool))
    print(f"\npairwise Jaccard: max {np.nanmax(off.values):.3f}, "
          f"median {np.nanmedian(off.values):.3f}")
    pairs = (off.stack().rename("jaccard").reset_index()
             .rename(columns={"level_0": "a", "level_1": "b"}))
    pairs = pairs[pairs["a"] < pairs["b"]].sort_values("jaccard", ascending=False)
    print("\nmost overlapping pairs (these are NOT independent readouts):")
    print(pairs.head(8).to_string(index=False, float_format="%.3f"))

    # ----------------------------------------------------------------- figure
    order = coll.order(list(mapped))
    # The size of each set goes in the tick label: Jaccard is a ratio, and 0.05 between a
    # 25-gene list and a 1,698-gene one does not mean what the same number means between
    # two lists of equal size. The count is the MAPPED one, which is what J is built on.
    labelled = J.loc[order, order].copy()
    labelled.index = labelled.columns = [f"{n} ({len(mapped[n])})" for n in order]

    fig, ax = plt.subplots(figsize=(8.5, 7))
    sns.heatmap(labelled, cmap="rocket_r", vmin=0, vmax=float(np.nanmax(off.values)),
                annot=True, fmt=".2f", annot_kws={"size": 7}, square=True,
                cbar_kws={"label": "Jaccard index", "shrink": 0.7}, ax=ax)
    # `pad` leaves room for the two block labels below the title: at the default the
    # title sits where "immune" and "stemness" are drawn and the three overlap.
    ax.set_title(f"Pairwise overlap of the {coll.title} collection\n"
                 "(mapped genes; the block structure is why these are not independent tests)",
                 fontsize=10, pad=30)
    # One separator per axis boundary, and the axis name centred over its block. On the EMT
    # collection there are three blocks rather than two, so both are derived from the order
    # rather than hard-coded.
    edges = coll.block_edges(order)
    for pos in edges:
        ax.axhline(pos, color="black", lw=2)
        ax.axvline(pos, color="black", lw=2)
    for start, stop in zip([0] + edges, edges + [len(order)]):
        ax.text((start + stop) / 2, -0.15, coll.axis_of[order[start]], ha="center",
                va="bottom", fontsize=9, weight="bold")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=8)
    C.savefig("jaccard_signature_overlap", "04_3_signatures", coll, fig, caveat=False)
    plt.close(fig)

    # coverage barplot
    fig, ax = plt.subplots(figsize=(8, 4.5))
    cov = coverage.loc[order]
    axis_colors = dict(zip(coll.axes, sns.color_palette("deep", len(coll.axes))))
    colors = [axis_colors[a] for a in cov["axis"]]
    ax.bar(range(len(cov)), cov["mapped_fraction"], color=colors)
    ax.axhline(C.MIN_MAPPED_FRACTION, color="crimson", ls="--", lw=1,
               label=f"floor ({C.MIN_MAPPED_FRACTION:.0%})")
    for i, (n_m, n_g) in enumerate(zip(cov["n_mapped"], cov["n_genes"])):
        ax.text(i, cov["mapped_fraction"].iloc[i] + 0.015, f"{n_m}/{n_g}",
                ha="center", fontsize=7)
    ax.set_xticks(range(len(cov)))
    ax.set_xticklabels(cov.index, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("fraction of symbols mapped")
    ax.set_ylim(0, 1.18)   # headroom: at 1.08 the legend box sits on the LIM_STEM bar label
    handles = [plt.Rectangle((0, 0), 1, 1, color=axis_colors[a]) for a in coll.axes
               if a in set(cov["axis"])]
    labels = [a for a in coll.axes if a in set(cov["axis"])]
    ax.set_title(f"{coll.title}: signature coverage on the epithelial object "
                 f"({n_all:,} genes)", fontsize=10)
    ax.legend(handles + [ax.get_lines()[0]], labels + [f"floor ({C.MIN_MAPPED_FRACTION:.0%})"],
              fontsize=8)
    sns.despine(ax=ax)
    C.savefig("signature_coverage", "04_3_signatures", coll, fig, caveat=False)
    plt.close(fig)

    print("\ndone.")


if __name__ == "__main__":
    main()
