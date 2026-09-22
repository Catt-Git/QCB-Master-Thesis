#!/usr/bin/env python3
"""05_8, the picture of a name: the programme and its axis, on the same UMAP.

`dr_naming_tum.py` ends with a barplot - one bar per programme, one segment per latent
direction the two routes agreed on. That plot says THAT a direction carries a programme; it
cannot say WHERE. This script draws the where, once per bar segment:

    left   the integrated UMAP, every cell coloured by its score for THAT programme
    right  the same UMAP, every cell coloured by its coordinate on THAT dimension,
           oriented so that the named side of the axis is the purple side

The two panels are the two routes of 05_8 made visual, and the claim they test is the one the
bar makes: if `DR 50+` is MP7, the cells that are purple on the left are the cells that are
purple on the right. A bar whose two panels do not overlap is a correlation that lives in a
handful of cells, or in one cohort, and the barplot cannot show that.

WHICH PAIRS ARE DRAWN. Not decided here. The pair list is read straight out of
`programme_dimensions_consensus_<run>.csv`, which is the table `plot_programme_bars` is
handed - so what this script draws is, by construction, the bars of that figure and nothing
else: 32 programme-direction pairs on 25 dimensions on the run this phase reports. Change the
naming rule in `dr_naming_tum.py` and this follows without an edit. It therefore also reads
`tables_pruned/` by default, the same PRUNE_VANISHED=1 control that step reads, and
`--keep-vanished` sends both to `tables/` together.

WHAT IS PLOTTED, EXACTLY.

  left    the WITHIN-COHORT z of the programme's score (`z_<programme>` in 05_6's per-cell
          score table), not the raw score. That is the quantity Route A correlates - absolute
          scores are not comparable across patients sequenced at different depths - so this
          panel shows the number the rho was computed on rather than a cousin of it.
          `--score raw` draws `score_<programme>` instead, for the one question z cannot
          answer: whether the programme is high in absolute terms or only high for its cohort.

  right   the cell's coordinate on the dimension, taken POSITIONALLY from `embed.X` against
          `embed.var['title']` - the same idiom `cell_first_tum.py` builds its `L` with, and
          the reason it is not `X[:, n]`: DRVI titles the dimensions by
          `reconstruction_effect`, so `DR 50` is not latent column 50.
          The coordinate is MULTIPLIED BY THE SIGN OF THE DIRECTION. A name is a claim about
          one side of an axis - `DR 18-` is a statement about the negative side - and the
          orientation is what makes the two panels comparable by eye. The colourbar says so.

  colour  ONE SIDE ONLY, on both panels. Everything at or below zero is the same flat grey
          and everything above it climbs a single-hue ramp, clipped at the 99th percentile
          of the positive values so that a dozen extreme cells cannot flatten 42,000. The
          left panel's ramp is purple; the right panel's is RED for a `+` direction and BLUE
          for a `-` one, so the colour names the side of the axis without the reader having
          to look at the sign.

          A diverging map was the first thing tried here and it was the wrong instrument: it
          spends half its range on cells that do not have the programme, which are exactly
          the cells neither route is making a claim about. Grey is the honest colour for
          "this cell is not in it", and it leaves the eye nothing to do but compare the two
          coloured sets.

          The left panel's zero is the cohort mean (the score is a within-cohort z) and the
          right panel's is the origin of the latent axis. Neither is a threshold: they are
          where the quantity stops pointing the way the name points.

The UMAP is the one 05_3 computed on the 64-dimensional latent itself
(`sc.pp.neighbors(embed, use_rep='X')`), i.e. the integrated embedding - not a UMAP of the
HVG space with the batch effect still in it.

NOTHING IS MEASURED HERE. The rho and the FDR in every subtitle are read from
`convergence_<collection>_<run>.csv` and `dim_signature_spearman_<collection>_<run>.csv`. The
one number computed on the spot is the Spearman of the PLOTTED signature against the PLOTTED
axis, and it is computed only because the two can differ: the programme on a bar is Route B's
best enrichment, `A_rho` is Route A's best correlation, and `verdict == convergent` requires
the same FAMILY, not the same signature. Where they differ the subtitle names both.

Usage:
    export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
    N_LATENT=64 HVG_SET=nomt python3 dr_umap_tum.py
    N_LATENT=64 HVG_SET=nomt python3 dr_umap_tum.py --only "DR 50+,DR 4+,DR 9+"
    N_LATENT=64 HVG_SET=nomt python3 dr_umap_tum.py --score raw --cmap Purples
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from types import SimpleNamespace

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.ticker
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

# Same reason, and the same two lines, as dr_naming_tum.py: $PRUNE_VANISHED is what moves
# signature_common's TABLE_DIR and figure directory, and the module resolves both at import
# time, so it has to be set before the import. This step's default is the pruned side
# because its input is that step's output; `--keep-vanished` moves the pair together.
_KEEP_VANISHED = "--keep-vanished" in sys.argv
if not _KEEP_VANISHED:
    os.environ["PRUNE_VANISHED"] = "1"

UTILS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils")
sys.path.insert(0, UTILS_DIR)
import signature_common as C          # noqa: E402
import sig_collections as SC          # noqa: E402

# The naming bar lives in ONE place, the step that applies it. Imported rather than repeated
# so that moving it there moves the number this script reports against. The import is below
# the $PRUNE_VANISHED line above on purpose: that module resolves the same environment at
# import time and must see the same value this one does.
from dr_naming_tum import NAMING_RHO_MIN as NAMING_BAR   # noqa: E402

STEP = "05_8_convergence"

# The pseudo-collection dr_naming writes under: this figure belongs to the run and not to any
# one collection, and it has to land beside the barplot it explains.
CONSENSUS = SimpleNamespace(name="consensus",
                            title="all three collections, read together")

# One flat grey for everything at or below zero, and one saturated ramp above it. The grey is
# light enough to sit under the coloured cells and dark enough not to read as background.
GREY = "#dcdcdf"

RAMPS = {
    "score": (GREY, "#b88cf0", "#8b2fe0", "#5b00b0"),   # the programme: purple
    "+":     (GREY, "#f59b86", "#e5352a", "#a00a00"),   # a positive direction: red
    "-":     (GREY, "#8fc0ef", "#1f6fd0", "#0a3577"),   # a negative direction: blue
}

# Clip the colour scale here rather than at the extreme: a handful of cells three times
# further out than the rest would otherwise leave the whole cloud one shade of grey. The
# quantile is taken over the POSITIVE values only - the clipped half would drag it to zero.
CLIP_Q = 0.99

POINT_SIZE = 2.4


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--keep-vanished", action="store_true",
                   help="read tables/ instead of tables_pruned/, i.e. the pair list of a "
                        "dr_naming run that kept the vanished tail")
    p.add_argument("--only", default="",
                   help="comma-separated directions to draw, e.g. 'DR 50+,DR 4+,DR 9+' "
                        "(default: every pair in the barplot)")
    p.add_argument("--programme", default="",
                   help="comma-separated programme names to draw, e.g. 'MP7_STRESS_IN_VITRO'")
    p.add_argument("--score", choices=("z", "raw"), default="z",
                   help="the left panel: within-cohort z (default, what Route A correlates) "
                        "or the raw score_genes value")
    p.add_argument("--dpi", type=int, default=200,
                   help="200 by default and not the phase's 300: these are 42,096-point "
                        "scatters and there are 32 of them")
    return p.parse_args()


def dim_n(s: str) -> int:
    return int(re.search(r"\d+", str(s)).group())


def read_pairs(args) -> pd.DataFrame:
    """The bars of `programmes_named`, one row per bar segment: programme x direction.

    Read rather than recomputed, so that the set drawn here cannot drift from the set
    plotted there. `collections` is a comma-joined list because a programme could in
    principle be named out of two collections at once; each one gets its own row, which is
    also how the bar counts them.
    """
    prog = C.read_table("programme_dimensions", CONSENSUS)
    rows = []
    for programme, r in prog.iterrows():
        for direction in [d.strip() for d in str(r["dimensions"]).split(",")]:
            for coll in [c.strip() for c in str(r["collections"]).split(",")]:
                rows.append({"programme": programme, "dim_direction": direction,
                             "collection": coll, "dimension": direction[:-1].strip(),
                             "sign": 1.0 if direction.endswith("+") else -1.0,
                             "max_rho": float(r["max_rho"])})
    d = pd.DataFrame(rows)

    if args.only:
        want = {s.strip() for s in args.only.split(",")}
        unknown = want - set(d["dim_direction"])
        if unknown:
            print(f"[note] not on the barplot, ignored: {', '.join(sorted(unknown))}")
        d = d[d["dim_direction"].isin(want)]
    if args.programme:
        want = {s.strip() for s in args.programme.split(",")}
        d = d[d["programme"].isin(want)]
    if d.empty:
        sys.exit("no pair left to draw after --only/--programme")
    return d.sort_values(["dim_direction", "programme"],
                         key=lambda s: s.map(C.dim_sort_key) if s.name == "dim_direction" else s)


def collection_by_slug(slug: str):
    """A collection by the SLUG its outputs are written under - `gavish_tnbc` included.

    `sig_collections.get` takes `--collection` values, and `gavish_tnbc` deliberately is not
    one of them: it is the narrowing of `gavish` that `resolve` returns, and only the slug
    reaches disk. This step reads its slugs off a table, where that choice was already made
    and cannot be re-made here, so it looks the object up by what the paths were built from.
    """
    for c in (*SC.COLLECTIONS.values(), *SC.GAVISH_VARIANTS.values()):
        if c.name == slug:
            return c
    raise KeyError(f"no collection writes under the slug {slug!r}")


def load_collection(name: str, which: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Per-cell scores, the convergence table and the Route A rho matrix, for one collection.

    The score table is the heavy one and is keyed on the OBJECT rather than on the run - see
    `signature_common.scores_csv` - so it is the same file every embedding reads.
    """
    coll = collection_by_slug(name)
    path = C.scores_csv(coll)
    if not path.exists():
        sys.exit(f"missing {path}: run 05_6_cell_first/cell_first_tum.py --collection {name}")
    scores = pd.read_csv(path, index_col=0, comment="#")
    prefix = "z_" if which == "z" else "score_"
    keep = [c for c in scores.columns if c.startswith(prefix)]
    scores = scores[keep]
    scores.columns = [c[len(prefix):] for c in scores.columns]

    conv = C.read_table("convergence", coll)
    rho = C.read_table("dim_signature_spearman", coll)
    print(f"[read] {name:12s} {scores.shape[0]} cells x {scores.shape[1]} readouts "
          f"({prefix.rstrip('_')}), {len(conv)} convergence rows")
    return scores, conv, rho


def dimension_frame(adata) -> pd.DataFrame:
    """Cells x dimensions, columns TITLED - 'DR 1' ... - which is not the column order.

    The same construction `cell_first_tum.py` uses for its `L`: positional against
    `var['title']`, never `X[:, n]`. DRVI orders the dimensions by `reconstruction_effect`
    and titles them in that order, so latent column n and `DR n` are different things.
    """
    return pd.DataFrame(np.asarray(adata.X), index=adata.obs_names,
                        columns=adata.var["title"].values)


def top_limit(v: np.ndarray) -> float:
    """The high end of the colour scale: the 99th percentile of the POSITIVE values."""
    pos = v[v > 0]
    if pos.size == 0:
        return 1.0
    return float(np.quantile(pos, CLIP_Q)) or float(pos.max()) or 1.0


def panel(ax, umap: np.ndarray, v: np.ndarray, ramp: str, label: str, tick_sign: float = 1.0):
    """One UMAP, grey below zero, `ramp` above it, the strongest cells drawn last.

    `v` arrives already clipped at zero, so `np.argsort` puts the whole grey half first and
    the draw order follows the colour. That order is not cosmetic: at 42,096 points in a few
    square inches the cells plotted first are simply not visible, and a random order would
    bury exactly the cells the figure is about under the cells it is not about.

    `tick_sign` is -1 for a negative direction, where the plotted quantity is the coordinate
    NEGATED: the scale then runs 0, -0.5, -1.0 and the colourbar says what the axis says.
    """
    cmap = LinearSegmentedColormap.from_list(f"one_sided_{ramp}", RAMPS[ramp])
    order = np.argsort(v)
    s = ax.scatter(umap[order, 0], umap[order, 1], c=v[order], cmap=cmap,
                   vmin=0.0, vmax=top_limit(v), s=POINT_SIZE, lw=0, rasterized=True)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_color("0.85")
    ax.set_xlabel("UMAP 1", fontsize=8, color="0.35")
    ax.set_ylabel("UMAP 2", fontsize=8, color="0.35")
    cb = ax.figure.colorbar(s, ax=ax, fraction=0.046, pad=0.02, extend="max")
    cb.set_label(label, fontsize=8)
    cb.ax.tick_params(labelsize=7)
    if tick_sign < 0:
        cb.ax.yaxis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda t, _: f"{-t if t else 0.0:.1f}"))
    cb.outline.set_visible(False)
    return s


def draw_pair(row, umap, L, scores, conv, rho_mat, args) -> None:
    programme, direction = row["programme"], row["dim_direction"]
    dim, sign = row["dimension"], row["sign"]

    if programme not in scores.columns:
        print(f"[skip] {direction} {programme}: no {args.score} column in "
              f"{row['collection']}'s score table")
        return
    # Both panels are clipped at zero: below it there is no claim to colour. On the right
    # the coordinate is multiplied by the sign of the DIRECTION first, so that what survives
    # the clip is the named side of the axis - `DR 32-` keeps its negative half.
    y = np.clip(scores[programme].to_numpy(float), 0.0, None)
    x = np.clip((sign * L[dim]).to_numpy(float), 0.0, None)

    # Read, not measured: the two routes' own numbers for this direction.
    c = conv.loc[direction]
    # Route A's rho for the signature ACTUALLY PLOTTED, oriented to this side. Route A is
    # computed per dimension, so the sign of the correlation is what makes it a statement
    # about one side; `A_rho` in the row above is the maximum over the collection and can
    # belong to a different signature of the same family.
    rho_plotted = float(sign * rho_mat.loc[dim, programme]) if programme in rho_mat.columns \
        else float("nan")
    # The same quantity recomputed on the two vectors in front of us. It should reproduce
    # `rho_plotted`; it is here because a figure that quotes a number it did not draw is how
    # the wrong cells get plotted under the right caption.
    rho_drawn = float(spearmanr(sign * L[dim].to_numpy(float),
                                scores[programme].to_numpy(float)).statistic)

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 5.0))

    score_label = ("within-cohort z" if args.score == "z" else "score_genes") \
        + f", {programme.replace('_', ' ')}   (<= 0 grey)"
    panel(axes[0], umap, y, "score", score_label)
    axes[0].set_title(f"cells: {programme.replace('_', ' ')}",
                      fontsize=10, loc="left", color="0.2")

    side = "+" if sign > 0 else "-"
    axis_label = (f"DRVI coordinate on {dim}   "
                  f"({'<= 0' if sign > 0 else '>= 0'} grey)")
    panel(axes[1], umap, x, side, axis_label, tick_sign=sign)
    axes[1].set_title(f"axis: {direction}", fontsize=10, loc="left", color="0.2")

    flags = str(c["confounder_flags"])
    bits = [f"Route A rho {rho_plotted:+.3f} on this signature",
            f"Route B FDR {float(c['B_fdr']):.1e}",
            f"named out of {row['collection']}"]
    if str(c["A_best_signature"]) != programme:
        bits.append(f"Route A's own best here: {c['A_best_signature']} "
                    f"(rho {float(c['A_rho']):+.3f})")
    if flags and flags != "none":
        bits.append(f"flags: {flags}")

    fig.suptitle(f"{programme.replace('_', ' ')}  on  {direction}", fontsize=13, x=0.012,
                 ha="left", y=1.045)
    fig.text(0.012, 1.005, " · ".join(bits), fontsize=8, color="0.35", ha="left", va="top")
    coloured = ("red where the coordinate is positive" if sign > 0
                else "blue where the coordinate is negative")
    fig.text(0.012, 0.978,
             f"{C.RUN_ID}, UMAP of the {L.shape[1]}-dimensional DRVI latent (05_3). "
             f"Left: purple where the programme's score is above its cohort mean.\n"
             f"Right: {coloured}. Grey is everything on the other side of zero. "
             f"Spearman of the two full, unclipped vectors: {rho_drawn:+.3f}.",
             fontsize=7.5, color="0.45", ha="left", va="top", linespacing=1.5)

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    C.savefig(f"umap_{C.dim_slug(direction)}_{programme}", STEP, CONSENSUS, fig=fig,
              dpi=args.dpi)
    plt.close(fig)

    return {"dim_direction": direction, "programme": programme,
            "collection": row["collection"], "rho_this_signature": rho_plotted,
            "A_best_signature": str(c["A_best_signature"]), "A_rho": float(c["A_rho"]),
            "B_fdr": float(c["B_fdr"]), "confounder_flags": flags}


def main():
    args = parse_args()
    C.banner("05_8 - the programme and its axis on the same UMAP")
    print(f"run {C.RUN_ID}; reading {C.TABLE_DIR.name}/ "
          f"({'vanished tail NOT read' if not args.keep_vanished else 'vanished tail read'})")

    pairs = read_pairs(args)
    print(f"{len(pairs)} programme-direction pairs on "
          f"{pairs['dim_direction'].nunique()} directions "
          f"({pairs['dimension'].nunique()} dimensions), "
          f"from {', '.join(sorted(pairs['collection'].unique()))}")

    emb = C.get_embedding(C.DEFAULT_EMBEDDING)
    adata = ad.read_h5ad(emb.embed_h5ad)
    if "X_umap" not in adata.obsm:
        sys.exit(f"{emb.embed_h5ad.name} has no X_umap: it was written before 05_3's "
                 f"neighbors/umap step")
    umap = np.asarray(adata.obsm["X_umap"])
    L = dimension_frame(adata)
    print(f"[read] {emb.embed_h5ad.name}: {adata.n_obs} cells, {adata.n_vars} dimensions, "
          f"UMAP on the latent itself")

    loaded = {}
    for name in sorted(pairs["collection"].unique()):
        scores, conv, rho_mat = load_collection(name, args.score)
        missing = L.index.difference(scores.index)
        if len(missing):
            sys.exit(f"{name}: {len(missing)} cells of the embedding have no score row")
        loaded[name] = (scores.loc[L.index], conv, rho_mat)

    C.banner(f"{len(pairs)} figures")
    drawn = []
    for _, row in pairs.iterrows():
        scores, conv, rho_mat = loaded[row["collection"]]
        rec = draw_pair(row, umap, L, scores, conv, rho_mat, args)
        if rec is not None:
            drawn.append(rec)

    report(pd.DataFrame(drawn))
    print(f"\n[done] {C.fig_dir(STEP, CONSENSUS)}")


def report(d: pd.DataFrame) -> None:
    """What the barplot cannot show: the rho of the programme ON THE BAR.

    A bar's `max_rho` is Route A's best correlation on that direction, over the whole
    collection. The programme the bar is LABELLED with is Route B's best enrichment, and
    `verdict == convergent` asks only that the two sit on the same FAMILY. Where they are
    two different signatures of one family, the number next to the bar belongs to a
    signature the bar does not name, and the panel drawn here - the programme on the bar,
    against the axis - is the one that carries the bar's own claim.

    No table is written. This is a reading of tables 05_8 already owns, and the figures are
    the output; a second file holding the same numbers would be a second thing to keep in
    step with the naming rule.
    """
    if d.empty:
        return
    d = d.copy()
    d["same"] = d["A_best_signature"] == d["programme"]
    C.banner("the programme on the bar, against the axis")
    print(d.sort_values("rho_this_signature", ascending=False)
          [["dim_direction", "programme", "collection", "rho_this_signature",
            "A_best_signature", "A_rho", "confounder_flags"]]
          .to_string(index=False, float_format="%.3f"))
    split = d[~d["same"]]
    print(f"\n{len(split)} of {len(d)} bars are labelled with a signature that is not the "
          f"one Route A scored highest on that direction")
    weak = split[split["rho_this_signature"] < NAMING_BAR]
    if len(weak):
        print(f"of those, {len(weak)} where the labelled programme's OWN rho is BELOW the "
              f"naming bar ({NAMING_BAR}), i.e. the direction is over the bar on a sibling "
              f"signature of the same family and not on the one the bar is labelled with:")
        for r in weak.sort_values("rho_this_signature").itertuples():
            print(f"  {r.dim_direction:8s} {r.programme:28s} rho {r.rho_this_signature:+.3f}"
                  f"   (bar quotes {r.A_best_signature} at {r.A_rho:+.3f})")


if __name__ == "__main__":
    main()
