#!/usr/bin/env python3
"""05_8, the reading step: how many latent dimensions actually get a name.

Route C writes one convergence table PER COLLECTION. This script is the step after it and
the only one in the phase that reads all three at once, because the question it answers is
not about a collection - it is about the LATENT SPACE:

    of the dimensions this run actually carries, how many can be given a name, and with
    what confidence?

Nothing here is a new measurement. Every number comes out of the three convergence tables
05_8 already wrote; what this adds is three decisions, each of which is a decision and is
written down:

1. WHICH DIMENSIONS ARE ON THE TABLE. The vanished tail is NOT read. This is the opposite
   of the phase default and it is deliberate for this step alone: naming a dimension is a
   claim about the model's own representation, and a dimension DRVI collapsed onto the prior
   carries ~1e-05 of the latent variance - it has no representation to name. The phase keeps
   those rows everywhere else so that they can be dismissed on their merits rather than
   silently; this step is where they are dismissed. It therefore reads `tables_pruned/`,
   the PRUNE_VANISHED=1 control, which is that decision already computed.
   `--keep-vanished` reads the unpruned tables instead and reports what the tail would add.

2. WHAT COUNTS AS A NAME. A convergent row, and nothing else: both routes, same signature
   family. Single-route rows are not names, per 05_8's own rule.

3. WHERE THE BAR IS - AND THAT IT IS A BAR, NOT A LADDER. A name is a yes or a no: a
   direction on which both routes clear their thresholds carries the programme, and a
   direction on which they do not, does not. Whether Route A's rho came out at 0.31 or at
   0.84 is a statement about how tightly the axis tracks the score, not about whether the
   name holds, and it is reported as a number (`A_rho`, `max_rho`) rather than used to sort
   the result into grades. THIS REPLACED A FOUR-TIER LADDER, deliberately: the tiers were
   read as a second filter, which they never were, and the bar is the only decision here.

   A bar has to be read against the null of the statistic it is applied to. `A_rho` is a
   DIMENSION's best gene set - a maximum over the collection's K sets - and
   `05_9/route_a_null_tum.py` measures exactly that under size- and expression-matched
   random sets:

       collection     K   median genes   null p50   p95    p97.5   p99    0.20 passes
       scie          10        247         0.079   0.267   0.309   0.362   10.9% of directions
       emt            9         27         0.065   0.187   0.223   0.275    3.7%
       gavish_tnbc   22         49         0.085   0.253   0.275   0.353    9.5%

   `NAMING_RHO_MIN = 0.30` sits at the p97.5 of that null on `scie`, between p97.5 and p99
   on `gavish_tnbc`, and above p99 on `emt` - i.e. 1-3 directions in 112 reach it by chance
   on the two collections that carry the result. That is the whole justification for the
   number, and it is why the bar is 0.30 here and not the phase-wide `C.ROUTE_A_RHO_MIN`.

   NOT the null quoted in `sig_collections` (p50 0.225, p95 0.325, p99 0.372). That one is a
   maximum over the 64 DIMENSION-directions for one gene list, which is the right reference
   for "can this list find anything" - the question MP38 and MP11 were dropped on - and the
   wrong one for a per-direction bar. Keeping the two apart is the point of the other script.

   THE BAR IS APPLIED HERE, TO A TABLE COMPUTED AT A LOWER ONE, and that is exact rather
   than a shortcut. `convergence_*.csv` stores `A_rho`, `B_significant` and `same_family`,
   and `convergence_tum.py` builds its verdict out of those three and nothing else; raising
   the bar can only turn a `convergent` row into a non-convergent one, never the reverse, so
   re-deriving the verdict at 0.30 from the stored columns gives the file a re-run at 0.30
   would write, without re-running 05_6 and 05_7. `main` refuses to do it the other way
   round - a stored bar ABOVE this one cannot be undone, because `same_family` was never
   computed for the rows it dropped.

A fourth thing this step computes rather than reads: the correlation of each DIMENSION
(not of a signature score) with sequencing depth, cycle scores and mitochondrial fraction.
The confounder flags of 05_6 are about the readout being claimed; this column is about the
axis carrying it, and on this run it is what separates DR 1 and DR 2 - the two depth axes
diagnosed in 05 - from the dimensions whose name is biology.

    NOTE, the trap this file exists downstream of: `DR n` is NOT latent column n. DRVI
    titles the dimensions by `reconstruction_effect`, so the map is
    `embed.var.set_index('title')['original_dim_id']` and it is what the loop below uses.

Usage:
    export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
    N_LATENT=64 HVG_SET=nomt python3 dr_naming_tum.py
    N_LATENT=64 HVG_SET=nomt python3 dr_naming_tum.py --keep-vanished
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
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

# Read before signature_common is imported, because $PRUNE_VANISHED is what moves its
# TABLE_DIR and its figure directory, and the module resolves both at import time. This
# step's DEFAULT is the pruned side - the one place in the phase where it is - so the flag
# is set here rather than expected from the environment, and `--keep-vanished` is what
# leaves it alone. An explicit $PRUNE_VANISHED=0 in the environment still loses to this,
# which is correct: the outputs of this step must not land beside the unpruned run's.
_KEEP_VANISHED = "--keep-vanished" in sys.argv
if not _KEEP_VANISHED:
    os.environ["PRUNE_VANISHED"] = "1"

UTILS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils")
sys.path.insert(0, UTILS_DIR)
import signature_common as C  # noqa: E402

# The three collections whose convergence tables exist for this run. A name is a name
# whichever of them supplied it, and a dimension named by two of them independently is
# the strongest statement this phase can make about one axis.
COLLECTIONS = ("scie", "emt", "gavish_tnbc")

# The naming bar, read against the PER-DIRECTION null of 05_9 - see the module docstring,
# point 3, for why it is 0.30 and why it is not `C.ROUTE_A_RHO_MIN`. `--rho-min` moves it.
NAMING_RHO_MIN = 0.30

# A dimension whose own coordinate tracks depth this hard is a depth axis, whatever the
# gene lists say about it. Same 0.30 convention as 05_8's DEPTH_FLAG, applied one level
# down - to the axis instead of to the readout.
DIM_DEPTH_FLAG = 0.30

# One colour, slot 1 of the validated categorical palette. One colour because there is one
# category: a direction either carries the programme or is not on this plot.
COL_BAR = "#2a78d6"

STEP = "05_8_convergence"
# A pseudo-collection: this table belongs to the run, not to any one collection, and the
# folder layout is the phase's so that it sits where every other table of this step does.
CONSENSUS = SimpleNamespace(name="consensus",
                            title="all three collections, read together")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--keep-vanished", action="store_true",
                   help="read tables/ instead of tables_pruned/, i.e. put DR 57+ back on "
                        "the table (default: the vanished tail is not read, see docstring)")
    p.add_argument("--rho-min", type=float, default=NAMING_RHO_MIN,
                   help=f"Route A's bar for a NAME, applied here to the stored convergence "
                        f"table (default {NAMING_RHO_MIN}, the p97.5 of the per-direction "
                        f"null; must be >= the bar that table was written at)")
    return p.parse_args()


def dim_n(s: str) -> int:
    return int(re.search(r"\d+", str(s)).group())


def read_convergence(pruned: bool) -> pd.DataFrame:
    """The three convergence tables, stacked, with the collection as a column.

    `C.TABLE_DIR` already follows $OUT_TAG, so a variant run - a different Route A bar, say -
    is read from its own folder without this script knowing there is such a thing as a
    variant. Only the unpruned comparison at the end names a directory, because that one is
    a comparison against a specific other run.
    """
    root = C.TABLE_DIR if pruned else C.PHASE_DIR / "tables"
    out = []
    for name in COLLECTIONS:
        path = root / name / C.RUN_ID / f"convergence_{name}_{C.RUN_ID}.csv"
        if not path.exists():
            print(f"[skip] {path} does not exist")
            continue
        d = pd.read_csv(path, comment="#")
        d["collection"] = name
        out.append(d)
    if not out:
        sys.exit(f"no convergence table under {root}; run 05_8 first")
    d = pd.concat(out, ignore_index=True)
    d["dim_n"] = d["dimension"].map(dim_n)
    return d


def dimension_covariates(embed) -> pd.DataFrame:
    """Per DIMENSION: depth, cycle and mitochondrial correlation, plus DRVI's own ranking.

    Computed here rather than read because no table in the phase holds it: 05_6 correlates
    the SIGNATURE SCORES with these covariates, which answers the other question.
    """
    X = np.asarray(embed.X)
    obs = embed.obs
    log_umi = np.log10(obs["total_counts"].to_numpy(float) + 1)
    covs = {
        "rho_log_umi": log_umi,
        "rho_n_genes": obs["n_genes_by_counts"].to_numpy(float),
        "rho_S_score": obs["S_score"].to_numpy(float),
        "rho_G2M_score": obs["G2M_score"].to_numpy(float),
        "rho_pct_mt": obs["pct_counts_mt"].to_numpy(float),
    }
    v = embed.var
    rows = []
    for i in range(X.shape[1]):
        r = {"dimension": v["title"].iloc[i],
             "latent_column": int(v["original_dim_id"].iloc[i]),
             "drvi_order": int(v["order"].iloc[i]),
             "reconstruction_effect": float(v["reconstruction_effect"].iloc[i]),
             "vanished": bool(v["vanished"].iloc[i]),
             "vanished_positive_direction": bool(v["vanished_positive_direction"].iloc[i]),
             "vanished_negative_direction": bool(v["vanished_negative_direction"].iloc[i])}
        r.update({k: float(spearmanr(X[:, i], y).statistic) for k, y in covs.items()})
        rows.append(r)
    d = pd.DataFrame(rows)
    d["latent_variance_share"] = d["reconstruction_effect"] / d["reconstruction_effect"].sum()
    d["dim_n"] = d["dimension"].map(dim_n)
    return d.set_index("dim_n").sort_index()


def main():
    args = parse_args()
    C.banner("05_8 - how many dimensions get a name")
    pruned = not args.keep_vanished
    print(f"run {C.RUN_ID}; reading {(C.TABLE_DIR if pruned else C.PHASE_DIR / 'tables').name}/ "
          f"({'vanished tail NOT read' if pruned else 'vanished tail read'})")
    if args.rho_min < C.ROUTE_A_RHO_MIN:
        sys.exit(f"--rho-min {args.rho_min} is BELOW the bar the convergence tables were "
                 f"written at ({C.ROUTE_A_RHO_MIN}). The rows it would add have no "
                 f"`same_family` - it was never computed for them - so they cannot be "
                 f"recovered here. Re-run 05_8 with ROUTE_A_RHO_MIN={args.rho_min}.")
    print(f"naming bar: Route A rho >= {args.rho_min} AND Route B FDR < 0.05, same family. "
          f"A name is a yes or a no; the magnitude is reported, not graded.")
    print(f"{args.rho_min} against the PER-DIRECTION null of 05_9: p97.5 on scie, ~p98 on "
          f"gavish_tnbc, above p99 on emt - NOT the per-list null of sig_collections")
    print(f"the convergence tables were written at rho >= {C.ROUTE_A_RHO_MIN}, so the bar "
          f"is re-derived here from A_rho / B_significant / same_family: exact, because "
          f"raising a threshold can only remove rows")

    conv = read_convergence(pruned)
    embed = C.get_embedding(C.DEFAULT_EMBEDDING)
    adata = ad.read_h5ad(embed.embed_h5ad)
    cov = dimension_covariates(adata)

    dims = sorted(conv["dim_n"].unique())
    print(f"\n{len(dims)} dimensions x 2 directions = {2 * len(dims)} rows, "
          f"x {conv['collection'].nunique()} collections = {len(conv)} verdicts")

    # The verdict AT THIS BAR, rebuilt from the three columns `convergence_tum.py` builds
    # its own out of. NOT `verdict == "convergent"`: that column is the stored bar's answer,
    # and the whole point of this step is that the bar is applied here.
    hits = conv[(conv["A_rho"] >= args.rho_min)
                & conv["B_significant"] & conv["same_family"]].copy()
    for c in ("rho_log_umi", "rho_S_score", "rho_G2M_score", "rho_pct_mt",
              "latent_variance_share", "drvi_order", "latent_column"):
        hits[c] = hits["dim_n"].map(cov[c])
    hits["axis_is_depth"] = hits["rho_log_umi"].abs() >= DIM_DEPTH_FLAG

    C.banner("the bar")
    stored = conv[conv["verdict"] == "convergent"]
    print(f"convergent at the stored bar ({C.ROUTE_A_RHO_MIN}): {len(stored):3d} rows on "
          f"{stored['dim_n'].nunique():2d} dimensions")
    print(f"convergent at this bar      ({args.rho_min}): {len(hits):3d} rows on "
          f"{hits['dim_n'].nunique():2d} dimensions")

    named = hits.copy()
    clean = named[~named.axis_is_depth]
    print(f"\nnamed:                          {named['dim_n'].nunique()} dimensions, "
          f"{named['dim_direction'].nunique()} directions")
    print(f"  of which the axis is depth:   "
          f"{named[named.axis_is_depth]['dim_n'].nunique()} "
          f"({', '.join(sorted(named[named.axis_is_depth]['dimension'].unique(), key=dim_n))})")
    print(f"named and not a depth axis:     {clean['dim_n'].nunique()} dimensions, "
          f"{clean['B_best_signature'].nunique()} distinct programmes")
    two = clean.groupby("dim_direction")["collection"].nunique()
    print(f"  named by two collections:     {(two > 1).sum()} directions "
          f"({', '.join(sorted(two[two > 1].index, key=lambda s: dim_n(s)))})")

    # ---------------------------------------------------------------- per-dimension table
    rows = []
    for n in dims:
        sub = hits[hits.dim_n == n].sort_values("A_rho", ascending=False)
        best = sub.iloc[0] if len(sub) else None
        # Every row of `sub` is over the bar now, so `strong` and `sub` are one set. The
        # name is kept because the calls below are about the AXIS, not about the row.
        strong = sub
        c = cov.loc[n]
        if len(strong) and abs(c["rho_log_umi"]) >= DIM_DEPTH_FLAG:
            call = "depth_axis"
        elif len(strong):
            call = "named"
        else:
            call = "unnamed"
        rows.append({
            "dimension": f"DR {n}",
            "drvi_order": int(c["drvi_order"]),
            "latent_column": int(c["latent_column"]),
            "latent_variance_share": c["latent_variance_share"],
            "call": call,
            "programme": best["B_best_signature"] if best is not None else "",
            "direction": best["dim_direction"] if best is not None else "",
            "A_rho": best["A_rho"] if best is not None else np.nan,
            "B_fdr": best["B_fdr"] if best is not None else np.nan,
            "same_signature": bool(best["same_signature"]) if best is not None else False,
            "named_by": ", ".join(sorted(strong["collection"].unique())),
            "n_convergent_rows": len(sub),
            "other_programmes": ", ".join(
                f"{r.dim_direction} {r.B_best_signature} (rho {r.A_rho:.2f})"
                for r in sub.iloc[1:].itertuples()),
            "confounder_flags": best["confounder_flags"] if best is not None else "",
            "rho_log_umi": c["rho_log_umi"], "rho_S_score": c["rho_S_score"],
            "rho_G2M_score": c["rho_G2M_score"], "rho_pct_mt": c["rho_pct_mt"],
        }
        )
    per_dim = pd.DataFrame(rows).set_index("dimension")
    C.write_table(per_dim, "dr_naming", CONSENSUS)

    C.banner("per dimension, DRVI order")
    show = ["drvi_order", "latent_variance_share", "call", "direction", "programme",
            "A_rho", "named_by", "rho_log_umi"]
    print(per_dim.sort_values("drvi_order")[show].to_string(float_format="%.3f"))

    # ---------------------------------------------------------------- programme table
    prog = (clean.groupby("B_best_signature")
            .agg(n_directions=("dim_direction", "nunique"),
                 max_rho=("A_rho", "max"), min_fdr=("B_fdr", "min"),
                 dimensions=("dim_direction",
                             lambda s: ", ".join(sorted(set(s), key=lambda x: dim_n(x)))),
                 collections=("collection", lambda s: ", ".join(sorted(set(s)))))
            .sort_values(["n_directions", "max_rho"], ascending=False))
    prog.index.name = "programme"
    C.write_table(prog, "programme_dimensions", CONSENSUS)
    C.banner(f"{len(prog)} programmes on {clean['dim_n'].nunique()} dimensions")
    print(prog.to_string(float_format="%.3g"))

    # ---------------------------------------------------------------- figure
    plot_programme_bars(prog, clean["dim_n"].nunique(), len(dims), args.rho_min)

    if pruned:
        C.banner("what not reading the vanished tail costs")
        un = read_convergence(pruned=False)
        tail = un[(un.dim_n > len(dims)) & (un.verdict == "convergent")].copy()
        strong = tail[(tail["A_rho"] >= args.rho_min) & tail["B_significant"]
                      & tail["same_family"]]
        print(f"{len(tail)} convergent rows on DR {len(dims) + 1}+ at the stored bar, of "
              f"which {len(strong)} clear this one:")
        if len(strong):
            print(strong[["dim_direction", "collection", "B_best_signature", "A_rho",
                          "B_fdr", "confounder_flags"]]
                  .sort_values("A_rho", ascending=False).to_string(index=False,
                                                                   float_format="%.3g"))
        lost = set(strong["B_best_signature"]) - set(prog.index)
        print(f"\nprogrammes the tail would add that the kept dimensions do not carry: "
              f"{', '.join(sorted(lost)) if lost else 'none'}")


def plot_programme_bars(prog: pd.DataFrame, n_dims_named: int, n_dims: int,
                        rho_min: float) -> None:
    """One horizontal bar per programme, one segment per direction that carries it.

    ONE COLOUR, on purpose. The bar was split into two rungs of a null-based ladder and that
    split was read as a grading of the names, which it never was: a direction over the bar
    carries the programme, and 0.31 and 0.84 are the same verdict told with different
    margins. The margin is still on the plot - `max_rho` is printed at the end of each bar -
    it is simply no longer the thing the geometry encodes.
    """
    p = prog.iloc[::-1]                       # largest at the top of a horizontal axis
    y = np.arange(len(p))
    fig, ax = plt.subplots(figsize=(10.5, 0.38 * len(p) + 1.9))
    ax.barh(y, p["n_directions"], height=0.46, color=COL_BAR)
    for i, (n, row) in enumerate(p.iterrows()):
        ax.text(row["n_directions"] + 0.10, i,
                f"{row['dimensions']}   (max rho {row['max_rho']:.2f})",
                va="center", ha="left", fontsize=8, color="0.30")
    ax.set_yticks(y)
    ax.set_yticklabels([n.replace("_", " ") for n in p.index], fontsize=9)
    ax.set_xlabel("latent directions on which the two routes converge", fontsize=9)
    ax.set_xlim(0, float(p["n_directions"].max()) + 2.6)
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    ax.grid(axis="x", color="0.88", lw=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color("0.7")
    ax.tick_params(length=0)
    ax.set_title(f"Programmes named on {C.RUN_ID}", fontsize=12, loc="left", pad=52)
    ax.text(0, 1.004,
            f"Route A rho >= {rho_min} and Route B FDR < 0.05 on the same signature family: "
            f"a name is a yes or a no,\nand both routes agree on "
            f"{n_dims_named} of the {n_dims} non-vanished dimensions. A bar counts one entry "
            f"per collection\nthat named the direction, so the two interferon rows share "
            f"DR 18-, DR 23- and DR 35-.",
            transform=ax.transAxes, fontsize=8.5, color="0.35", va="bottom",
            linespacing=1.45)
    fig.tight_layout()
    C.savefig("programmes_named", STEP, CONSENSUS, fig=fig)
    plt.close(fig)


if __name__ == "__main__":
    main()
