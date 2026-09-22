#!/usr/bin/env python3
"""05_10: the standard pipeline against the factor pipeline, on the same cells and the same
gene sets. How many programmes does each one RECOVER?

05_9 asked whether a programme is PRESENT in a coordinate system and answered yes for
Harmony - on most readouts more cleanly than for DRVI. That answer is real and it is the
reason this step exists: presence is not the question the phase turned on. The question is
whether a pipeline, run on these cells with no prior naming of an axis, ARRIVES at the
programme. So this step does not compare two spaces. It compares two ways of turning cells
into a gene list, and then runs the identical ORA on both.

    standard   space -> k-NN graph -> Leiden -> one-vs-rest DE -> top N genes -> ORA
    factor     DRVI decoder -> top N genes per dimension-direction -> ORA

Everything downstream of "top N genes" is one code path. The gene sets, the background, the
list depth, the hypergeometric test and the Benjamini-Hochberg treatment are the same object
in both arms, so the only thing that differs is WHERE THE GENE LIST CAME FROM. That is the
whole design, and it is what lets the difference be attributed.

THE ARMS, AND WHY THERE ARE THREE OF THEM

  harmony_leiden   05_3b's corrected components, clustered      the pipeline anyone would run
  drvi_leiden      05_3's latent space, clustered the same way  the control that separates
                                                                the SPACE from the UNIT
  drvi_axes        05_3's decoder, read per dimension-direction the phase's own route (05_7)
  pca_leiden       the uncorrected PCA, clustered               optional, what integration bought

`drvi_leiden` is the arm that makes this an experiment rather than a demonstration. If it
lands with `harmony_leiden` and both land far below `drvi_axes`, the loss is a property of
CLUSTERING and not of Harmony - which is the honest reading and the stronger claim, because
it does not require Harmony to be the worse method at anything it promises. If instead
`drvi_leiden` sits above `harmony_leiden`, the space mattered too, and the gap between them
is how much.

THE TWO OUTCOMES, AND WHY THE OBVIOUS ONE IS USELESS

The first thing anyone measures is "did some gene list of this arm enrich for this
signature at FDR < 0.05". On these cells that number SATURATES: it is 17/17 for the
metaprogram collection in both arms, because a 200-gene marker list of a 2,000-gene
background overlaps something in a 22-set catalogue almost always. That saturation is not a
nuisance to be tuned away, it is the first result of the step, and it is reported as
`recall_any` precisely so it cannot be quoted as agreement between the two pipelines.

The outcome that discriminates is the one that matches what naming a dimension MEANS in
05_8:

    a unit of analysis - a cluster or a dimension-direction - is IDENTIFIED BY the
    signature it enriches for most strongly, and a programme counts as NAMED by an arm when
    some unit of that arm has it as its best hit.

`recall_named` is the fraction of the testable catalogue that ends up as some unit's
identity, and it is the headline. The gap between the two numbers is the whole phenomenon:
an arm whose `recall_any` is 1.00 and whose `recall_named` is 0.59 did not find seventeen
programmes, it found ten coarse identities each of which is significant for a dozen
overlapping lists. `mean_hits_per_list` is the same statement from the other side, and
`gene_list_hits` is the table it is counted in.

Neither number is a ranking of integration methods, and neither is charged to Harmony: the
`drvi_leiden` arm runs the identical clustering on DRVI's own space, so whatever the
clustering costs is visible without Harmony in the picture at all.

THE TARGET LIST, AND WHY IT IS NOT "THE PROGRAMMES DRVI NAMED"

Recall needs a denominator that neither pipeline chose. Using 05_8's 19 named programmes
would be circular: they are the output of one of the arms. The denominator here is instead a
property of the CATALOGUE and the OBJECT - every signature of the collection with at least
`C.MIN_SIGNATURE_GENES` genes inside the ORA background, i.e. every set that either arm
could in principle enrich for. It is computed before either arm runs, it is identical for
all of them, and it is printed.

That denominator is generous to the standard pipeline in one direction and harsh in another,
and both are stated rather than corrected: generous because a signature no cell in this
object expresses is in the denominator and neither arm can find it, harsh because a set that
is genuinely absent counts against both equally. `--presence-filter` narrows it to the
signatures 05_9 found present in the Harmony space, when that table is on disk; the headline
number is the unfiltered one.

THE BUDGET OBJECTION, WHICH IS REAL

The factor arm brings 2 x n_latent gene lists to the test and Leiden at resolution 0.2 brings
eight. More lists is more chances to hit, and also a bigger BH denominator. Neither effect is
argued away here; the resolution scan is what answers it empirically. It runs up to
`--resolutions 3.0`, which puts the cluster count in the same range as the direction count,
and every row of the output table carries `n_gene_lists` and `n_tests`, so recall can be read
against the budget that bought it. A standard pipeline that is still below the factor arm at
matched list count has not been starved of tests.

THE NULL, WHICH IS NOT OPTIONAL

A gene list of 200 genes drawn from a 2,000-gene background overlaps a 50-gene signature by
chance often enough that "some cluster was significant for something" is not evidence. The
floor is a SIZE-MATCHED RANDOM PARTITION: the cluster labels of the reference resolution
shuffled among the cells, so the number of groups and every group size are preserved and
nothing else is. It goes through the identical DE and ORA. A random partition recovering
programmes at the same rate as a real one would mean the DE is reading list length rather
than biology, and the step says so.

WHAT THIS STEP DOES NOT CLAIM

It does not rank Harmony and DRVI as integration methods - that is phase 02, on what they
both promise - and it does not say the standard pipeline is wrong. A cluster is the right
unit for a question about discrete populations, and this object has none: it is one
compartment of one lineage, and its structure is continuous by construction. What the step
measures is what that costs, in programmes, on these cells.

Usage:
    export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
    cd 05_drvi_tumoral_epi/05_10_pipeline_recall

    N_LATENT=64 python3 pipeline_recall_tum.py --collection gavish
    N_LATENT=64 python3 pipeline_recall_tum.py --collection scie --arms harmony_leiden drvi_axes
    N_LATENT=64 python3 pipeline_recall_tum.py --resolutions 0.4 1.0 2.0 --no-null

`N_LATENT=64` is not optional in practice: 05_3b only ever ran at 64, so at any other value
the Harmony embedding this step needs is not on disk. `CELL_SET`, `HVG_SET` and
`PRUNE_VANISHED` behave exactly as in 05_4 - 05_9.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
import warnings

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import gseapy as gp
from statsmodels.stats.multitest import multipletests

UTILS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils")
sys.path.insert(0, UTILS_DIR)
import signature_common as C  # noqa: E402
import sig_collections as SC  # noqa: E402

CS = C.CS

STEP = "05_10_pipeline_recall"
FDR = 0.05                       # as 05_7, and applied the same way
K_NEIGHBOURS = 15                # as 05_3, 05_3b and 05_9: the graph is the phase's graph
N_TOP_GENES = 200                # as 05_7: the list depth is shared by both arms
DEFAULT_RESOLUTIONS = (0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0)
REF_RESOLUTION = 1.0

HEADER_NOTE = ("recall of a signature collection by two pipelines on ONE set of cells: "
               "Leiden clusters vs DRVI decoder directions, identical ORA downstream")

# One arm, one colour, fixed so dropping an arm never recolours the others.
ARM_COLOUR = {
    "harmony_leiden": "#2f6f9f",
    "drvi_leiden": "#7fb3d5",
    "pca_leiden": "#7a7a7a",
    "drvi_axes": "#c25e00",
    "drvi_axes_matched": "#e8993f",
    "null_partition": "#b0b0b0",
}
ARM_LABEL = {
    "harmony_leiden": "Harmony + Leiden + DE",
    "drvi_leiden": "DRVI + Leiden + DE",
    "pca_leiden": "uncorrected PCA + Leiden + DE",
    "drvi_axes": "DRVI decoder directions",
    "drvi_axes_matched": "DRVI decoder, list budget matched",
    "null_partition": "size-matched random partition",
}
LEIDEN_ARMS = ("harmony_leiden", "drvi_leiden", "pca_leiden")
ARM_SPACE = {"harmony_leiden": "harmony", "drvi_leiden": "drvi", "pca_leiden": "pca"}
INK = "#222222"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    SC.add_argument(p)
    p.add_argument("--arms", nargs="+",
                   default=["harmony_leiden", "drvi_leiden", "drvi_axes"],
                   choices=sorted(ARM_LABEL),
                   help="which pipelines to run (default: the three that make the comparison "
                        "an experiment; add pca_leiden for the uncorrected space)")
    p.add_argument("--resolutions", nargs="+", type=float, default=list(DEFAULT_RESOLUTIONS),
                   help="the Leiden resolution scan (default %(default)s). The top of the "
                        "range is there to match the factor arm's list budget")
    p.add_argument("--ref-resolution", type=float, default=REF_RESOLUTION,
                   help="the resolution the per-programme and per-cluster tables are written "
                        "at (default %(default)s); must be in --resolutions")
    p.add_argument("--k", type=int, default=K_NEIGHBOURS,
                   help=f"neighbours in the graph (default {K_NEIGHBOURS}, as 05_3/05_3b)")
    p.add_argument("--n-top-genes", type=int, default=N_TOP_GENES,
                   help=f"genes per gene list, BOTH arms (default {N_TOP_GENES}, as 05_7)")
    p.add_argument("--de-fdr", type=float, default=0.05,
                   help="BH-adjusted p a one-vs-rest DE gene has to clear to enter a "
                        "cluster's gene list (default %(default)s). Without it a random "
                        "partition scores like a real one - see de_top_genes")
    p.add_argument("--de-min-lfc", type=float, default=0.0,
                   help="log fold change floor on the same genes (default %(default)s, i.e. "
                        "up-regulated is enough; Seurat's own default is 0.25)")
    p.add_argument("--no-null", action="store_true",
                   help="skip the size-matched random partition, i.e. leave the floor "
                        "unmeasured")
    p.add_argument("--presence-filter", action="store_true",
                   help="narrow the denominator to the signatures 05_9 found present in the "
                        "Harmony space, when that table is on disk")
    p.add_argument("--overwrite", action="store_true",
                   help="recompute the partitions and DE lists instead of reusing the "
                        "cache in $DATA_DIR/05_tum/pipeline_lists_<arm>_<run>.json")
    return p.parse_args()


# --------------------------------------------------------------------------- #
# The spaces
# --------------------------------------------------------------------------- #

def load_coords(name: str, harmony_run: str, drvi_run: str):
    """(cells x dimensions indexed BY CELL NAME, a one-line description) for one space.

    Realignment is by cell name everywhere in this phase and by position nowhere, which is
    the guard 05_3b's own DRVI comparison and 05_9 both use. The uncorrected PCA is the one
    array with no names of its own; it borrows the Harmony embedding's index and is DROPPED
    rather than lined up on a guess if the lengths disagree.
    """
    tum = CS.tum_dir()
    if name in ("harmony", "pca"):
        path = tum / f"embed_{harmony_run}.h5ad"
        if not path.exists():
            sys.exit(f"missing {path}: run 05_3b_harmony_run/run_harmony_tum.py first")
        e = sc.read_h5ad(path)
        cells = e.obs_names.astype(str)
        if name == "harmony":
            X = pd.DataFrame(np.asarray(e.X), index=cells)
            return X, f"05_3b: Harmony, {X.shape[1]} corrected components"
        npy = tum / f"pca_{harmony_run}.npy"
        if not npy.exists():
            print(f"[skip] {npy.name} is not on disk; the uncorrected arm is dropped",
                  flush=True)
            return None, ""
        raw = np.load(npy)
        if raw.shape[0] != len(cells):
            print(f"[skip] {npy.name} has {raw.shape[0]} rows against {len(cells)} cells in "
                  f"the Harmony embedding; the uncorrected arm is dropped", flush=True)
            return None, ""
        X = pd.DataFrame(raw, index=cells)
        return X, f"the same {X.shape[1]} components BEFORE the correction"

    path = tum / f"embed_{drvi_run}.h5ad"
    if not path.exists():
        sys.exit(f"missing {path}: run 05_3_drvi_run first")
    e = sc.read_h5ad(path)
    X = pd.DataFrame(np.asarray(e.X), index=e.obs_names.astype(str))
    return X, f"05_3: DRVI, {X.shape[1]} latent dimensions"


# --------------------------------------------------------------------------- #
# The two ways of producing a gene list
# --------------------------------------------------------------------------- #

def leiden_labels(adata, coords: pd.DataFrame, resolutions, k: int, seed: int, arm: str):
    """{resolution: labels}, one Leiden partition per resolution on one space.

    The graph is built the way 05_3 and 05_3b build theirs - k = 15, the whole space as the
    representation, no rescaling of the coordinates - so "neighbours" means here what it
    means in those two steps' UMAPs and in 05_9's Moran's I.
    """
    rep = f"X_{arm}"
    adata.obsm[rep] = coords.reindex(adata.obs_names.astype(str)).values.astype(np.float32)
    t0 = time.time()
    sc.pp.neighbors(adata, n_neighbors=k, use_rep=rep, random_state=seed)
    print(f"    graph: k = {k} on {coords.shape[1]} dimensions "
          f"({time.time() - t0:.0f}s)", flush=True)

    out = {}
    for res in resolutions:
        key = f"leiden_{arm}_{res}"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sc.tl.leiden(adata, resolution=res, key_added=key, flavor="igraph",
                         n_iterations=2, directed=False, random_state=seed)
        lab = adata.obs[key].astype(str)
        out[res] = lab
        print(f"    resolution {res:<5} -> {lab.nunique():3d} clusters "
              f"(smallest {lab.value_counts().min():,}, largest {lab.value_counts().max():,})",
              flush=True)
    return out


def de_top_genes(adata, labels: pd.Series, n_top: int, de_fdr: float, min_lfc: float,
                 key: str = "_grp"):
    """({cluster: [genes]}, per-cluster DE stats) from one-vs-rest Wilcoxon.

    THE SIGNIFICANCE FILTER IS NOT OPTIONAL, and the first run of this step is why. Handing
    the ORA the top `n_top` genes of every cluster REGARDLESS of whether the cluster has any
    marker at all gives a size-matched RANDOM partition the same recall as a real one - it
    scored 10/10 on `scie` against 10/10 for both real arms - because 200 arbitrary genes of
    a 2,000-gene background overlap a 100-gene signature about as often as 200 real ones do.
    That is not a property of clustering, it is the instrument failing to discriminate, and
    it is exactly what the null partition is in this step to catch.

    So a gene enters a cluster's list only if the one-vs-rest test calls it: BH-adjusted
    p < `de_fdr` and log fold change > `min_lfc`, the two filters any real pipeline applies
    before it opens Enrichr. `n_top` then truncates what is left, so a cluster with a
    thousand markers is still read at the same depth as a decoder direction. A cluster with
    NO marker contributes an empty list and tests nothing, which is the correct behaviour
    and is counted: `n_de_genes` is reported per cluster and the arm is not compensated for
    it.

    UP-REGULATED ONLY, and truncated by test statistic rather than by fold change. A
    cluster's ORA question is "what is this cluster's programme", so a gene it is DEPLETED
    of is not part of the answer - the same asymmetry the factor arm has, where the two
    directions of a dimension are read separately and each one's top genes are the ones it
    drives UP.

    The DE runs on the 2,000-gene HVG panel, not on the 25,133-gene object, and that is a
    decision: it is the feature set DRVI was trained on, so it is the only background both
    arms can share. Giving the cluster arm the whole transcriptome would give it genes the
    decoder could never have proposed and a different hypergeometric universe, and the
    comparison would stop being one.
    """
    adata.obs[key] = pd.Categorical(labels.reindex(adata.obs_names.astype(str)).values)
    t0 = time.time()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.tl.rank_genes_groups(adata, groupby=key, method="wilcoxon", use_raw=False,
                                n_genes=adata.n_vars, key_added="_de")
        df = sc.get.rank_genes_groups_df(adata, group=None, key="_de")
    keep = (df["logfoldchanges"] > min_lfc) & (df["pvals_adj"] < de_fdr)
    up = df[keep]
    out, stats = {}, []
    for g in adata.obs[key].cat.categories:
        g = str(g)
        grp = up[up["group"].astype(str) == g]
        out[g] = grp.nlargest(n_top, "scores")["names"].tolist()
        stats.append({"cluster": g, "n_de_genes": int(len(grp)),
                      "n_genes_tested": int(len(out[g]))})
    stats = pd.DataFrame(stats)
    empty = int((stats["n_genes_tested"] == 0).sum())
    print(f"    DE: {len(out)} clusters, FDR < {de_fdr} and lfc > {min_lfc}, "
          f"median {stats['n_de_genes'].median():.0f} markers per cluster, "
          f"{empty} with none, truncated at {n_top} ({time.time() - t0:.0f}s)", flush=True)
    return out, stats


def axis_top_genes(embed, gene_names, n_top: int):
    """{dimension-direction: [genes]} off the DRVI decoder. 05_7's own list, rebuilt here.

    Deliberately recomputed rather than read from `C.top_genes_tsv`: the point of this step
    is that both arms go through ONE code path from the gene list onward, and a list read
    from a file written by another run is a second provenance nobody can check from here.
    It is the same call 05_7 makes, so the lists are identical when the run ids match.
    """
    scores = C.interpretability_scores(embed, gene_names, key=C.SCORE_KEY)
    out = {d: scores[d].nlargest(n_top).index.tolist() for d in scores.columns}
    print(f"    decoder: {len(out)} dimension-directions x top {n_top} genes", flush=True)
    return out



# --------------------------------------------------------------------------- #
# The cache: a partition and its DE lists depend on the SPACE, never on the collection
# --------------------------------------------------------------------------- #
#
# Three collections run over the same cells, and the Leiden partitions and their one-vs-rest
# DE are identical in all three - the collection enters only at the ORA. Recomputing them per
# collection is three times the work for a byte-identical result and three chances for two of
# them to disagree about what "resolution 1.0 on the Harmony space" means. The cache is keyed
# on the parameters that CAN change it (k, list depth, the two DE filters) and is discarded
# whole when any of them moves, so a stale list cannot survive a flag change.

def cache_paths(arm: str, run: str):
    tum = CS.tum_dir()
    return (tum / f"pipeline_lists_{arm}_{run}.json",
            tum / f"pipeline_partitions_{arm}_{run}.csv.gz")


def cache_params(args, description: str = "") -> dict:
    return {"k": args.k, "n_top_genes": args.n_top_genes, "de_fdr": args.de_fdr,
            "de_min_lfc": args.de_min_lfc, "description": description}


def load_cache(arm: str, run: str, args, cells):
    """({resolution: (lists, stats)}, {resolution: labels}, description), or empties.

    The partitions are realigned on cell NAME and the cache is dropped if a single cell of
    this object is missing from it - the same guard every other cross-file read in the phase
    uses, and for the same reason.
    """
    jpath, cpath = cache_paths(arm, run)
    if args.overwrite or not jpath.exists() or not cpath.exists():
        return {}, {}, ""
    try:
        blob = json.loads(jpath.read_text())
    except Exception:                                   # noqa: BLE001 - a corrupt cache is a
        return {}, {}, ""                               # miss, never a crash
    want = cache_params(args)
    have = dict(blob.get("params", {}))
    desc = str(have.pop("description", ""))
    if {k: v for k, v in want.items() if k != "description"} != have:
        print(f"    [cache] {jpath.name}: parameters changed, recomputing", flush=True)
        return {}, {}, ""
    parts_df = pd.read_csv(cpath, index_col=0)
    parts_df.index = parts_df.index.astype(str)
    if pd.Index(cells).difference(parts_df.index).size:
        print(f"    [cache] {cpath.name}: does not cover these cells, recomputing", flush=True)
        return {}, {}, ""
    lists, parts = {}, {}
    for res_s, payload in blob.get("resolutions", {}).items():
        if res_s not in parts_df.columns:
            continue
        res = float(res_s)
        lists[res] = (payload["lists"], pd.DataFrame(payload["stats"]))
        parts[res] = parts_df[res_s].reindex(cells).astype(str)
    if lists:
        print(f"    [cache] {jpath.name}: {len(lists)} resolutions reused "
              f"({', '.join(str(r) for r in sorted(lists))})", flush=True)
    return lists, parts, desc


def save_cache(arm: str, run: str, args, description: str, lists: dict, parts: dict):
    jpath, cpath = cache_paths(arm, run)
    blob = {"params": cache_params(args, description),
            "resolutions": {str(r): {"lists": l, "stats": st.to_dict(orient="records")}
                            for r, (l, st) in sorted(lists.items())}}
    jpath.write_text(json.dumps(blob))
    pd.DataFrame({str(r): parts[r] for r in sorted(parts)}).to_csv(cpath)
    print(f"    [cache] wrote {jpath.name} and {cpath.name}", flush=True)


# --------------------------------------------------------------------------- #
# The ORA, identical for every arm
# --------------------------------------------------------------------------- #

def run_ora(lists: dict[str, list[str]], sets: dict[str, list[str]], background: list[str],
            arm: str, resolution) -> pd.DataFrame:
    """Every gene list against every signature: offline hypergeometric, declared background.

    `gp.enrich`, as 05_7 - never Enrichr's implicit all-human-genes universe. Benjamini-
    Hochberg is NOT applied here: it is applied once per arm-and-resolution by `apply_bh`,
    across every pair actually tested including the ones with no overlap, which is the only
    denominator that makes two arms of different widths comparable.
    """
    records = []
    for name, genes in lists.items():
        try:
            res = gp.enrich(gene_list=genes, gene_sets=sets, background=background,
                            outdir=None)
        except Exception as exc:                       # noqa: BLE001 - reported, not raised
            print(f"      {name}: ORA FAILED ({exc})", flush=True)
            continue
        res_df = getattr(res, "results", None)
        if res_df is None or len(res_df) == 0:
            continue                                   # no overlap at all: p = 1, added below
        df = pd.DataFrame(res_df).copy()
        df.insert(0, "gene_list", str(name))
        records.append(df)
    if not records:
        return pd.DataFrame(columns=["gene_list", "Term", "P-value", "arm", "resolution"])
    long = pd.concat(records, ignore_index=True)
    long["arm"] = arm
    long["resolution"] = resolution
    return long


def apply_bh(long: pd.DataFrame, n_lists: int, n_sets: int) -> tuple[pd.DataFrame, int]:
    """BH across all `n_lists x n_sets` pairs of ONE arm at ONE resolution.

    The pairs gseapy left out are the ones with no overlap. They are p = 1 and can never
    become significant, but they belong in the denominator, so they are added back - exactly
    as 05_7 does it. Corrected WITHIN an arm, because that is the correction each pipeline
    would actually apply to itself: an arm is not penalised for the other one's width, and
    an arm that buys its hits with more tests pays for them in its own denominator.
    """
    n_pairs = n_lists * n_sets
    if long.empty:
        return long.assign(fdr_bh=pd.Series(dtype=float),
                           significant=pd.Series(dtype=bool)), n_pairs
    pvals = np.concatenate([long["P-value"].values, np.ones(max(0, n_pairs - len(long)))])
    _, padj, _, _ = multipletests(pvals, alpha=FDR, method="fdr_bh")
    out = long.copy()
    out["fdr_bh"] = padj[:len(long)]
    out["significant"] = out["fdr_bh"] < FDR
    return out, n_pairs


# --------------------------------------------------------------------------- #
# The diagnostic that explains the gap
# --------------------------------------------------------------------------- #

def eta_squared(y: np.ndarray, codes: np.ndarray, n_groups: int) -> float:
    """Share of a per-cell score's variance that lies BETWEEN clusters. One-way eta^2.

    This is the number the whole argument rests on, and it is computed on 05_6's cached
    within-cohort z scores - the same values 05_9 measured presence on, never re-scored
    here. A programme with a high eta^2 is a programme the partition has a group for; a
    programme with a low one is a gradient ACROSS the groups, and a gradient is what a
    one-vs-rest DE has nothing to test. If the signatures the cluster arm misses are the
    low-eta^2 ones, the miss is explained rather than merely counted.
    """
    y = np.asarray(y, dtype=np.float64)
    grand = y.mean()
    ss_tot = float(((y - grand) ** 2).sum())
    if ss_tot <= 0:
        return float("nan")
    counts = np.bincount(codes, minlength=n_groups).astype(np.float64)
    sums = np.bincount(codes, weights=y, minlength=n_groups)
    ok = counts > 0
    means = np.zeros_like(counts)
    means[ok] = sums[ok] / counts[ok]
    ss_between = float((counts[ok] * (means[ok] - grand) ** 2).sum())
    return ss_between / ss_tot


def cluster_profile(labels: pd.Series, cohorts: pd.Series) -> pd.DataFrame:
    """Per cluster: size, and how much of it is one patient.

    `cohort_dominance` is here because the standard pipeline's other well-known failure is
    not about gradients at all: a cluster that is 90% one cohort makes its DE a patient
    contrast, and any programme it enriches for is that patient's. It is reported rather
    than filtered - filtering it would be deciding for the pipeline which of its clusters
    were allowed to count.
    """
    df = pd.DataFrame({"cluster": labels.values, "cohort": cohorts.reindex(labels.index).values})
    rows = []
    for cl, grp in df.groupby("cluster", observed=True):
        frac = grp["cohort"].value_counts(normalize=True)
        rows.append({"cluster": str(cl), "n_cells": len(grp),
                     "n_cohorts": int(grp["cohort"].nunique()),
                     "cohort_dominance": float(frac.iloc[0]),
                     "top_cohort": str(frac.index[0])})
    return pd.DataFrame(rows).sort_values("n_cells", ascending=False)


# --------------------------------------------------------------------------- #
# Output
# --------------------------------------------------------------------------- #

def write_table(df: pd.DataFrame, name: str, coll, run_id: str, index: bool = False):
    """C.table_path for the path, this step's own header for the first lines.

    Not `C.write_table`, for 05_9's reason: that helper stamps "<Space> run <run_id>" from
    the EMBEDDINGS registry, which is right for a table measured in one space and wrong for
    one that compares pipelines across two.
    """
    path = C.table_path(name, coll, run_id)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"# {name} | collection {coll.name} ({coll.title}) | run {run_id} "
                 f"| {HEADER_NOTE} | {STEP}\n")
        for line in C.CAVEAT.split(". "):
            if line.strip():
                fh.write(f"# CAVEAT: {line.strip().rstrip('.')}.\n")
        df.to_csv(fh, index=index)
    print(f"[table] {path}  ({df.shape[0]} x {df.shape[1]})", flush=True)
    return path


def fig_recall_curve(summary: pd.DataFrame, coll, run_id: str, n_testable: int, args):
    """Recall against the resolution, and against the list budget that bought it.

    Two panels because the budget objection is real: the left one is the scan as it was run,
    the right one puts the same points on the axis the objection is about, so a reader can
    check whether the cluster arms are merely under-tested.
    """
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))
    leiden = summary[summary["unit"] == "cluster"]
    flat = summary[summary["unit"] == "axis"]

    for ax, xcol, xlabel in ((axes[0], "resolution", "Leiden resolution"),
                             (axes[1], "n_gene_lists", "gene lists handed to the ORA")):
        for arm, grp in leiden.groupby("arm", observed=True):
            grp = grp.sort_values(xcol)
            ax.plot(grp[xcol], grp["recall_named"], "-o", ms=4, lw=1.8,
                    color=ARM_COLOUR.get(arm, INK), label=ARM_LABEL.get(arm, arm))
            ax.plot(grp[xcol], grp["recall_any"], ":", lw=1.1, alpha=0.75,
                    color=ARM_COLOUR.get(arm, INK))
        for _, row in flat.iterrows():
            ax.axhline(row["recall_named"], color=ARM_COLOUR.get(row["arm"], INK), lw=1.8,
                       ls="--", label=ARM_LABEL.get(row["arm"], row["arm"]))
            ax.axhline(row["recall_any"], color=ARM_COLOUR.get(row["arm"], INK), lw=1.0,
                       ls=":", alpha=0.75)
            if xcol == "n_gene_lists":
                ax.plot([row["n_gene_lists"]], [row["recall_named"]], "D", ms=6,
                        color=ARM_COLOUR.get(row["arm"], INK))
        null = summary[summary["arm"] == "null_partition"]
        if len(null):
            ax.axhline(float(null["recall_named"].iloc[0]), color=ARM_COLOUR["null_partition"],
                       lw=1.4, ls="-.", label=ARM_LABEL["null_partition"])
        ax.set_xlabel(xlabel)
        ax.set_ylim(-0.03, 1.03)
        ax.grid(alpha=0.25, lw=0.6)
        if xcol == "n_gene_lists":
            ax.set_xscale("log")

    axes[0].set_ylabel(f"fraction of the {n_testable} testable signatures")
    handles, labels = axes[0].get_legend_handles_labels()
    seen, h, l = set(), [], []
    for hh, ll in zip(handles, labels):
        if ll not in seen:
            seen.add(ll); h.append(hh); l.append(ll)
    axes[0].legend(h, l, fontsize=7.5, loc="lower right", frameon=False)
    fig.suptitle(f"{coll.title}: what each pipeline recovers\n"
                 f"identical ORA downstream - top {args.n_top_genes} genes per list, "
                 f"hypergeometric against the {n_testable}-set collection on the 2,000-HVG "
                 f"background, BH within arm at FDR < {FDR}.\n"
                 "SOLID: named, i.e. some unit is identified by it.   "
                 "DOTTED: any significant hit - the saturated metric, drawn to be dismissed",
                 fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    C.savefig("recall_by_resolution", STEP, coll, fig, run_id=run_id)
    plt.close(fig)


def fig_recovery_matrix(rec: pd.DataFrame, coll, run_id: str, arms: list[str], ref_res):
    """Signature x arm, coloured by -log10 FDR. The table the recall number summarises."""
    cols = [f"neglog_fdr__{a}" for a in arms if f"neglog_fdr__{a}" in rec.columns]
    if not cols:
        return
    order = coll.order(rec.index.tolist(), for_figure=True)
    order = order + [s for s in rec.index if s not in order]
    M = rec.loc[order, cols].astype(float)
    fig, ax = plt.subplots(figsize=(C.fig_span(len(cols), 1.5, 3.4),
                                    C.fig_span(len(order), 0.26, 2.2)))
    thr = -np.log10(FDR)
    im = ax.imshow(M.values, aspect="auto", cmap="Blues",
                   vmin=0, vmax=max(thr * 2, float(np.nanmax(M.values)) or thr * 2))
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M.values[i, j]
            named = bool(rec.loc[order[i], f"named__{cols[j].split('__', 1)[1]}"]) \
                if f"named__{cols[j].split('__', 1)[1]}" in rec.columns else False
            if np.isfinite(v) and v >= thr:
                ax.text(j, i, "*" if named else ".", ha="center", va="center",
                        fontsize=13 if named else 15,
                        color="#08306b" if named else "#7f7f7f")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([ARM_LABEL.get(c.split("__", 1)[1], c) for c in cols],
                       rotation=25, ha="right", fontsize=8)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=7)
    ax.set_title(f"{coll.title}: which signature each pipeline reached\n"
                 f"Leiden arms at resolution {ref_res}. "
                 f"*  some unit is IDENTIFIED BY it (its best hit)    "
                 f".  significant but never the best hit of any unit", fontsize=9)
    fig.colorbar(im, ax=ax, shrink=0.4, label="-log10 FDR (best gene list)")
    C.savefig("recovery_matrix", STEP, coll, fig, run_id=run_id)
    plt.close(fig)


def fig_eta2(rec: pd.DataFrame, coll, run_id: str, arm: str, ref_res):
    """The explanation: what the cluster arm could not name, against how clustered it is.

    x is eta^2(cluster) of the signature's own per-cell score at the reference resolution -
    how much of it lies BETWEEN clusters rather than across them. One row per signature,
    sorted, because the identities are the point: a strip plot of seventeen metaprogram names
    at two y values is unreadable, and the first version of this figure was.

    A signature low on this axis that the arm could not name is a GRADIENT the partition had
    no group for. One high on it that the arm could not name is NOT explained by this figure
    and must not be claimed to be - MP3_CELL_CYCLE_HMG_RICH is the case on this run.
    """
    xcol, ycol = f"eta2_cluster__{arm}", f"named__{arm}"
    if xcol not in rec.columns or ycol not in rec.columns:
        return
    d = rec[[xcol, ycol]].dropna(subset=[xcol]).sort_values(xcol)
    if d.empty:
        return
    named = d[ycol].astype(bool).values
    y = np.arange(len(d))
    fig, ax = plt.subplots(figsize=(8.0, C.fig_span(len(d), 0.30, 2.4)))
    ax.hlines(y, 0, d[xcol].values, color="0.85", lw=1.2, zorder=1)
    ax.scatter(d.loc[named, xcol], y[named], s=70, marker="o",
               color=ARM_COLOUR.get(arm, INK), edgecolor="white", linewidth=0.8, zorder=3,
               label=f"names a unit ({int(named.sum())})")
    ax.scatter(d.loc[~named, xcol], y[~named], s=80, marker="X", color="#c0392b",
               edgecolor="white", linewidth=0.8, zorder=3,
               label=f"never a unit's best hit ({int((~named).sum())})")
    med_hit = float(d.loc[named, xcol].median()) if named.any() else np.nan
    med_miss = float(d.loc[~named, xcol].median()) if (~named).any() else np.nan
    for v, colour, lab in ((med_hit, ARM_COLOUR.get(arm, INK), "median, named"),
                           (med_miss, "#c0392b", "median, not named")):
        if np.isfinite(v):
            ax.axvline(v, color=colour, lw=1.2, ls="--", alpha=0.7, zorder=2,
                       label=f"{lab}: {v:.3f}")
    ax.set_yticks(y)
    ax.set_yticklabels(d.index, fontsize=7.5)
    ax.set_ylim(-0.8, len(d) - 0.2)
    ax.set_xlim(left=0)
    ax.set_xlabel(r"$\eta^2$(cluster) of the per-cell score  "
                  "- share of its variance lying BETWEEN clusters")
    ax.grid(alpha=0.25, lw=0.6, axis="x")
    ax.legend(fontsize=7.5, frameon=False, loc="lower right")
    ax.set_title(f"{coll.title}: {ARM_LABEL.get(arm, arm)} at resolution {ref_res}\n"
                 "a programme the partition has no group for is a programme a one-vs-rest "
                 "DE cannot propose", fontsize=9)
    fig.tight_layout()
    C.savefig(f"eta2_vs_recovery_{arm}", STEP, coll, fig, run_id=run_id)
    plt.close(fig)


def fig_hits_per_list(per_list: pd.DataFrame, coll, run_id: str, arms: list[str], ref_res):
    """How many signatures one gene list is significant for: the entanglement read.

    A unit that hits eight signatures at once has not resolved eight programmes, it has
    found the block they share. This is the same statement 05_9 makes with `best_dim` - H1
    being the best correlate of all five immunogenicity lists - counted on the discovery
    side instead of the correlation side.
    """
    d = per_list[per_list["arm"].isin(arms)]
    if d.empty:
        return
    fig, ax = plt.subplots(figsize=(7.6, 4.2))
    order = [a for a in arms if a in set(d["arm"])]
    width = 0.8 / max(1, len(order))
    mx = int(d["n_signatures_hit"].max())
    bins = np.arange(0, mx + 2)
    for i, arm in enumerate(order):
        v = d.loc[d["arm"] == arm, "n_signatures_hit"].values
        counts, _ = np.histogram(v, bins=bins)
        ax.bar(bins[:-1] + i * width - 0.4 + width / 2, counts / max(1, counts.sum()),
               width=width, color=ARM_COLOUR.get(arm, INK),
               label=f"{ARM_LABEL.get(arm, arm)}  (n = {len(v)}, mean {v.mean():.2f})")
    ax.set_xticks(bins[:-1])
    ax.set_xlabel("signatures one gene list is significant for")
    ax.set_ylabel("fraction of that arm's gene lists")
    ax.legend(fontsize=7.5, frameon=False)
    ax.grid(alpha=0.25, lw=0.6, axis="y")
    ax.set_title(f"{coll.title}: how many programmes one unit of analysis carries\n"
                 f"Leiden arms at resolution {ref_res}", fontsize=9)
    fig.tight_layout()
    C.savefig("hits_per_gene_list", STEP, coll, fig, run_id=run_id)
    plt.close(fig)


# --------------------------------------------------------------------------- #

def main():
    args = parse_args()
    coll = SC.resolve(args)
    drvi_run = CS.run_id(method="drvi")
    harmony_run = CS.run_id(method="harmony")
    run_id = drvi_run          # names the CELLS and the gene sets, not a space; see 05_9

    CS.banner(f"{STEP} - what does a standard pipeline recover? ({coll.name})")
    print(f"collection : {coll.title}")
    print(f"question   : {coll.question}")
    print(f"arms       : {', '.join(args.arms)}")
    print(f"reference  : {drvi_run}   |   harmony: {harmony_run}")
    print(f"resolutions: {args.resolutions}   (tables at {args.ref_resolution})")
    print(f"list depth : top {args.n_top_genes} genes, both arms", flush=True)

    if args.ref_resolution not in args.resolutions:
        sys.exit(f"--ref-resolution {args.ref_resolution} is not in --resolutions "
                 f"{args.resolutions}: the per-programme tables would have no partition")

    # ---- the gene sets, and the denominator neither arm chose ---------------
    gmt = C.gmt_path(coll)
    if not gmt.exists():
        sys.exit(f"missing {gmt}: run 05_4_signatures/build_signatures_tum.py "
                 f"--collection {args.collection} first")
    sets_all = C.read_gmt(gmt)

    hvg = ad.read_h5ad(C.HVG_H5AD)
    hvg.obs_names = hvg.obs_names.astype(str)
    background = list(hvg.var_names)
    bg = set(background)
    sets_bg = {k: [g for g in v if g in bg] for k, v in sets_all.items()}
    testable = {k: v for k, v in sets_bg.items() if len(v) >= C.MIN_SIGNATURE_GENES}
    dropped = {k: len(v) for k, v in sets_bg.items() if k not in testable}

    C.banner("the denominator: every signature either arm could in principle enrich for")
    print(f"{hvg.n_obs:,} cells x {len(background):,} HVGs - the DRVI training feature set, "
          f"the DE feature set and the ORA background, all three")
    print(f"{len(sets_all)} signatures in the collection -> {len(testable)} with at least "
          f"{C.MIN_SIGNATURE_GENES} genes inside it")
    for k, n in dropped.items():
        print(f"  dropped: {k:32s} {n:3d} genes in background")

    if args.presence_filter:
        pres = C.table_path(f"signature_presence", coll, run_id)
        if not pres.exists():
            print(f"[skip] --presence-filter: {pres.name} is not on disk, the denominator "
                  f"stays the testable set", flush=True)
        else:
            t = pd.read_csv(pres, comment="#")
            keep = set(t.loc[(t["space"] == "harmony") & (t["morans_i_vs_random_z"] >= 3),
                             "readout"])
            before = len(testable)
            testable = {k: v for k, v in testable.items() if k in keep}
            print(f"--presence-filter: {before} -> {len(testable)} signatures 05_9 found "
                  f"present in the Harmony space at >= 3 sd over the matched random level")

    if not testable:
        sys.exit("no signature is testable on this background; nothing to measure")
    n_testable = len(testable)

    # ---- the cells, the cohorts, the cached per-cell scores -----------------
    cohorts = hvg.obs[C.BATCH_KEY].astype(str)
    scores_csv = C.scores_csv(coll)
    z = None
    if scores_csv.exists():
        s = pd.read_csv(scores_csv, index_col=0)
        s.index = s.index.astype(str)
        z = s[[c for c in s.columns if c.startswith("z_")]].rename(columns=lambda c: c[2:])
        z = z.reindex(hvg.obs_names)
        print(f"\n[read] {scores_csv.name}: {z.shape[1]} per-cell readouts, "
              f"05_6's within-cohort z - re-scored nowhere in this step", flush=True)
    else:
        print(f"\n[skip] {scores_csv.name} is not on disk, so eta^2(cluster) cannot be "
              f"computed; run 05_6 for the diagnostic half of this step", flush=True)

    # ---- the gene lists, arm by arm -----------------------------------------
    # {(arm, resolution): {list name: [genes]}}, plus the partitions the cluster arms used.
    gene_lists: dict[tuple, dict] = {}
    de_stats: dict[tuple, pd.DataFrame] = {}
    partitions: dict[tuple, pd.Series] = {}
    descriptions: dict[str, str] = {}

    for arm in args.arms:
        if arm not in LEIDEN_ARMS:
            continue
        C.banner(f"arm {arm}")
        cached, cparts, desc = load_cache(arm, drvi_run, args, hvg.obs_names)
        todo = [r for r in args.resolutions if r not in cached]
        if todo:
            coords, space_desc = load_coords(ARM_SPACE[arm], harmony_run, drvi_run)
            if coords is None:
                print(f"[skip] arm {arm}: its space could not be read", flush=True)
                continue
            absent = pd.Index(hvg.obs_names).difference(coords.index)
            if len(absent):
                sys.exit(f"{arm}: {len(absent)} cells of the object are not in its embedding "
                         f"(e.g. {list(absent[:3])}). Two pipelines compared on two cell "
                         f"sets is not a comparison")
            desc = space_desc
            print(f"    {desc}", flush=True)
            labs = leiden_labels(hvg, coords, todo, args.k, C.SEED, arm)
            for res, lab in labs.items():
                cparts[res] = lab.astype(str)
                cached[res] = de_top_genes(hvg, lab, args.n_top_genes, args.de_fdr,
                                           args.de_min_lfc)
            save_cache(arm, drvi_run, args, desc, cached, cparts)
        else:
            print(f"    {desc or ARM_LABEL.get(arm, arm)}", flush=True)
        descriptions[arm] = desc
        for res in args.resolutions:
            gene_lists[(arm, res)], de_stats[(arm, res)] = cached[res]
            partitions[(arm, res)] = cparts[res]

    if "drvi_axes" in args.arms:
        embed = ad.read_h5ad(CS.tum_dir() / f"embed_{drvi_run}.h5ad")
        n_van = C.n_vanished(embed)
        desc = (f"05_3: DRVI decoder, {embed.n_vars} dimensions x 2 directions "
                + (f"({n_van} vanished PRUNED)" if C.PRUNE_VANISHED
                   else f"({n_van} flagged vanished and NOT pruned)"))
        descriptions["drvi_axes"] = desc
        C.banner(f"arm drvi_axes: {desc}")
        assert embed.varm[f"{C.SCORE_KEY}_positive"].shape[1] == len(background), \
            "the embedding and the DRVI input disagree on the gene axis"
        gene_lists[("drvi_axes", np.nan)] = axis_top_genes(embed, hvg.var_names,
                                                           args.n_top_genes)

    # THE BUDGET OBJECTION, CLOSED RATHER THAN ARGUED. The decoder brings 2 x n_latent lists
    # and the widest Leiden partition brings fewer. This arm is the decoder cut down to the
    # widest cluster count in the run, taking DRVI's OWN top dimensions - `var['order']`, the
    # reconstruction-effect ranking, both directions of each - so it is the subset the model
    # itself would have offered first, not one chosen to win. If it still names what the
    # cluster arms cannot at the same number of tests, the gap is not a test budget.
    if ("drvi_axes", np.nan) in gene_lists:
        widest = max((len(v) for (a, _), v in gene_lists.items() if a in LEIDEN_ARMS),
                     default=0)
        full = gene_lists[("drvi_axes", np.nan)]
        if 0 < widest < len(full):
            keep = list(full)[:widest]                 # already in `order`, +/- interleaved
            gene_lists[("drvi_axes_matched", np.nan)] = {k: full[k] for k in keep}
            descriptions["drvi_axes_matched"] = (
                f"the first {widest} of the {len(full)} decoder directions, in DRVI's own "
                f"reconstruction-effect order - the widest Leiden partition of this run")
            print(f"\n[budget control] drvi_axes_matched: {widest} directions, "
                  f"matching the widest cluster count in the run", flush=True)

    if not args.no_null and partitions:
        ref_arm = next((a for a in args.arms if (a, args.ref_resolution) in partitions), None)
        if ref_arm is not None:
            C.banner("the floor: the reference partition shuffled among the cells")
            lab = partitions[(ref_arm, args.ref_resolution)]
            rng = np.random.default_rng(C.SEED)
            shuffled = pd.Series(rng.permutation(lab.values), index=lab.index)
            print(f"    {shuffled.nunique()} groups, every group size preserved, "
                  f"membership destroyed", flush=True)
            descriptions["null_partition"] = (
                f"the {ref_arm} partition at resolution {args.ref_resolution}, shuffled")
            gene_lists[("null_partition", args.ref_resolution)], \
                de_stats[("null_partition", args.ref_resolution)] = de_top_genes(
                    hvg, shuffled, args.n_top_genes, args.de_fdr, args.de_min_lfc)
            partitions[("null_partition", args.ref_resolution)] = shuffled

    if not gene_lists:
        sys.exit("no arm produced a gene list; nothing to test")

    # ---- one ORA, every arm -------------------------------------------------
    C.banner(f"ORA: offline hypergeometric against the {len(background):,}-HVG background, "
             f"BH within arm")
    summary_rows, per_list_rows, long_all = [], [], []

    for (arm, res), lists in gene_lists.items():
        long = run_ora(lists, testable, background, arm, res)
        long, n_pairs = apply_bh(long, len(lists), n_testable)
        long_all.append(long)

        hits = long[long["significant"]] if len(long) else long
        recovered = sorted(set(hits["Term"])) if len(hits) else []
        # A unit is IDENTIFIED BY its strongest enrichment; the distinct identities are the
        # arm's vocabulary. This is 05_8's own definition of naming a dimension, applied
        # unchanged to a cluster.
        best_per_list = (hits.loc[hits.groupby("gene_list", observed=True)["fdr_bh"].idxmin()]
                         if len(hits) else pd.DataFrame(columns=["gene_list", "Term"]))
        named = sorted(set(best_per_list["Term"])) if len(best_per_list) else []
        per_list_hits = (hits.groupby("gene_list", observed=True)["Term"].nunique()
                         if len(hits) else pd.Series(dtype=int))
        per_list_hits = per_list_hits.reindex(list(lists), fill_value=0)

        prof = None
        if (arm, res) in partitions:
            prof = cluster_profile(partitions[(arm, res)], cohorts)

        list_sizes = np.array([len(v) for v in lists.values()])
        row = {
            "arm": arm,
            "description": descriptions.get(arm, ""),
            "unit": "axis" if arm.startswith("drvi_axes") else "cluster",
            "resolution": res,
            "n_gene_lists": len(lists),
            "n_gene_lists_empty": int((list_sizes == 0).sum()),
            "median_genes_per_list": float(np.median(list_sizes)),
            "n_tests": n_pairs,
            "n_signatures_testable": n_testable,
            "n_signatures_recovered": len(recovered),
            "recall_any": len(recovered) / n_testable,
            "n_signatures_named": len(named),
            "recall_named": len(named) / n_testable,
            "n_units_with_a_name": int(len(best_per_list)),
            "n_significant_pairs": int(len(hits)),
            "mean_hits_per_list": float(per_list_hits.mean()),
            "max_hits_per_list": int(per_list_hits.max()) if len(per_list_hits) else 0,
            "frac_lists_with_no_hit": float((per_list_hits == 0).mean()),
            "median_best_fdr": (float(hits.groupby("Term")["fdr_bh"].min().median())
                                if len(hits) else float("nan")),
            "recovered": "; ".join(recovered),
            "named": "; ".join(named),
        }
        if prof is not None:
            row.update({
                "min_cluster_cells": int(prof["n_cells"].min()),
                "median_cluster_cells": float(prof["n_cells"].median()),
                "max_cohort_dominance": float(prof["cohort_dominance"].max()),
                "n_clusters_one_cohort_over_80pct":
                    int((prof["cohort_dominance"] > 0.8).sum()),
            })
        summary_rows.append(row)

        n_de = ({r["cluster"]: r["n_de_genes"]
                 for _, r in de_stats[(arm, res)].iterrows()}
                if (arm, res) in de_stats else {})
        for name, n in per_list_hits.items():
            per_list_rows.append({"arm": arm, "resolution": res, "gene_list": str(name),
                                  "n_genes_in_list": len(lists[name]),
                                  "n_de_genes": n_de.get(str(name), np.nan),
                                  "n_signatures_hit": int(n),
                                  "n_cells": (int((partitions[(arm, res)] == str(name)).sum())
                                              if (arm, res) in partitions else np.nan)})

        print(f"  {arm:16s} res {str(res):<5} {len(lists):>4} lists  "
              f"named {len(named):>3}/{n_testable} (recall_named {row['recall_named']:.3f})  "
              f"any {len(recovered):>3}/{n_testable}  "
              f"mean hits/list {row['mean_hits_per_list']:.2f}", flush=True)

    summary = pd.DataFrame(summary_rows).sort_values(["unit", "arm", "resolution"])
    # `recall` is kept as an alias of the headline so nothing downstream has to know which
    # of the two it wanted; both columns are in the table and the figures name them.
    summary["recall"] = summary["recall_named"]
    per_list = pd.DataFrame(per_list_rows)
    long_all = pd.concat(long_all, ignore_index=True) if long_all else pd.DataFrame()

    # ---- the union over the whole scan, the most generous reading -----------
    #
    # A practitioner picks ONE resolution. This row gives the cluster arms every resolution
    # at once - the union of everything they name anywhere in the scan - which is an upper
    # bound no single run of that pipeline could reach, and is labelled as one. A programme
    # absent from this union is a programme the standard pipeline does not reach on these
    # cells at any granularity it was given.
    union_rows = []
    for arm, grp in summary.groupby("arm", observed=True):
        named_union = set()
        for v in grp["named"].dropna():
            named_union |= {t for t in str(v).split("; ") if t}
        union_rows.append({
            "arm": arm,
            "n_settings": int(len(grp)),
            "n_signatures_named_union": len(named_union),
            "recall_named_union": len(named_union) / n_testable,
            "best_single_setting": float(grp["recall_named"].max()),
            "named_union": "; ".join(sorted(named_union)),
            "never_named": "; ".join(sorted(set(testable) - named_union)),
        })
    union = pd.DataFrame(union_rows).sort_values("recall_named_union", ascending=False)

    # ---- the per-programme table, at the reference resolution ---------------
    ref = {}
    extra = (["drvi_axes_matched"] if ("drvi_axes_matched", np.nan) in gene_lists else []) \
        + (["null_partition"] if not args.no_null else [])
    for arm in args.arms + extra:
        key = (arm, np.nan) if arm.startswith("drvi_axes") else (arm, args.ref_resolution)
        if key in gene_lists:
            ref[arm] = key

    rec = pd.DataFrame(index=sorted(testable))
    rec.index.name = "signature"
    rec["n_genes_in_collection"] = [len(sets_all[s]) for s in rec.index]
    rec["n_genes_in_background"] = [len(testable[s]) for s in rec.index]

    for arm, key in ref.items():
        sub = long_all[(long_all["arm"] == key[0]) &
                       ((long_all["resolution"] == key[1])
                        if not pd.isna(key[1]) else long_all["resolution"].isna())]
        best = (sub.loc[sub.groupby("Term")["fdr_bh"].idxmin()].set_index("Term")
                if len(sub) else pd.DataFrame())
        fdr = best["fdr_bh"].reindex(rec.index) if len(best) else pd.Series(index=rec.index,
                                                                           dtype=float)
        rec[f"recovered__{arm}"] = (fdr < FDR).fillna(False)
        rec[f"best_fdr__{arm}"] = fdr
        rec[f"neglog_fdr__{arm}"] = -np.log10(fdr.clip(lower=1e-300))
        rec[f"best_list__{arm}"] = (best["gene_list"].reindex(rec.index)
                                    if len(best) else pd.NA)
        sig = sub[sub["significant"]] if len(sub) else sub
        rec[f"n_lists_hitting__{arm}"] = (
            sig.groupby("Term")["gene_list"].nunique().reindex(rec.index).fillna(0)
            .astype(int) if len(sig) else 0)
        bpl = (sig.loc[sig.groupby("gene_list", observed=True)["fdr_bh"].idxmin()]
               if len(sig) else pd.DataFrame(columns=["gene_list", "Term"]))
        n_naming = (bpl.groupby("Term")["gene_list"].nunique().reindex(rec.index).fillna(0)
                    .astype(int) if len(bpl) else 0)
        rec[f"n_lists_naming__{arm}"] = n_naming
        rec[f"named__{arm}"] = (n_naming > 0) if len(bpl) else False

    # eta^2(cluster), the diagnostic
    eta_rows = []
    if z is not None:
        for (arm, res), lab in partitions.items():
            codes, uniq = pd.factorize(lab.reindex(hvg.obs_names).values)
            for sig in rec.index:
                if sig not in z.columns:
                    continue
                y = z[sig].values
                ok = np.isfinite(y)
                eta_rows.append({"arm": arm, "resolution": res, "signature": sig,
                                 "n_clusters": len(uniq),
                                 "eta2_cluster": eta_squared(y[ok], codes[ok], len(uniq))})
        eta = pd.DataFrame(eta_rows)
        for arm, key in ref.items():
            if arm.startswith("drvi_axes"):
                continue
            sub = eta[(eta["arm"] == key[0]) & (eta["resolution"] == key[1])]
            rec[f"eta2_cluster__{arm}"] = sub.set_index("signature")["eta2_cluster"].reindex(
                rec.index)
    else:
        eta = pd.DataFrame(columns=["arm", "resolution", "signature", "n_clusters",
                                    "eta2_cluster"])

    # ---- tables --------------------------------------------------------------
    C.banner("tables")
    write_table(summary, "pipeline_recall", coll, run_id, index=False)
    write_table(union, "pipeline_recall_union", coll, run_id, index=False)
    write_table(rec.reset_index(), "programme_recovery", coll, run_id, index=False)
    write_table(per_list, "gene_list_hits", coll, run_id, index=False)
    if len(eta):
        write_table(eta, "programme_cluster_eta2", coll, run_id, index=False)

    comp_rows = []
    for arm, key in ref.items():
        if key not in partitions:
            continue
        prof = cluster_profile(partitions[key], cohorts).assign(arm=arm, resolution=key[1])
        hits = long_all[(long_all["arm"] == key[0]) & (long_all["resolution"] == key[1])
                        & long_all["significant"]]
        by_cl = hits.groupby("gene_list", observed=True)["Term"].apply(
            lambda s: "; ".join(sorted(set(s)))) if len(hits) else pd.Series(dtype=str)
        if key in de_stats:
            prof = prof.merge(de_stats[key], on="cluster", how="left")
        prof["signatures_hit"] = prof["cluster"].map(by_cl).fillna("")
        prof["n_signatures_hit"] = prof["cluster"].map(
            by_cl.str.count("; ").add(1) if len(by_cl) else pd.Series(dtype=float)).fillna(0)
        comp_rows.append(prof)
    if comp_rows:
        write_table(pd.concat(comp_rows, ignore_index=True), "cluster_composition", coll,
                    run_id, index=False)

    if len(long_all):
        write_table(long_all.drop(columns=["Genes"], errors="ignore"), "ora_all_pairs", coll,
                    run_id, index=False)

    # ---- figures -------------------------------------------------------------
    C.banner("figures")
    fig_recall_curve(summary, coll, run_id, n_testable, args)
    fig_recovery_matrix(rec, coll, run_id, list(ref), args.ref_resolution)
    for arm in ref:
        if arm in LEIDEN_ARMS:
            fig_eta2(rec, coll, run_id, arm, args.ref_resolution)
    fig_hits_per_list(per_list[(per_list["resolution"] == args.ref_resolution)
                               | per_list["resolution"].isna()],
                      coll, run_id, list(ref), args.ref_resolution)

    # ---- the read ------------------------------------------------------------
    C.banner("what came out")
    print(f"  denominator: {n_testable} testable signatures of {len(sets_all)}\n")
    best = (summary[summary["unit"] == "cluster"].groupby("arm")["recall_named"].max()
            if (summary["unit"] == "cluster").any() else pd.Series(dtype=float))
    for arm, r in best.items():
        row = summary[(summary["arm"] == arm) & (summary["recall_named"] == r)].iloc[0]
        print(f"  {ARM_LABEL.get(arm, arm):34s} named {int(row['n_signatures_named']):>3}"
              f"/{n_testable} (best, at resolution {row['resolution']} with "
              f"{int(row['n_gene_lists'])} clusters)   any {int(row['n_signatures_recovered'])}"
              f"/{n_testable}   {row['mean_hits_per_list']:.2f} signatures per cluster")
    for _, row in summary[summary["unit"] == "axis"].iterrows():
        print(f"  {ARM_LABEL.get(row['arm'], row['arm']):34s} named "
              f"{int(row['n_signatures_named']):>3}/{n_testable} (with "
              f"{int(row['n_gene_lists'])} directions)   any "
              f"{int(row['n_signatures_recovered'])}/{n_testable}   "
              f"{row['mean_hits_per_list']:.2f} signatures per direction")

    print("\n  union over the whole scan - every resolution at once, an upper bound no "
          "single run\n  of the clustering pipeline could reach:")
    for _, r in union.iterrows():
        never = r["never_named"]
        print(f"    {ARM_LABEL.get(r['arm'], r['arm']):34s} "
              f"{int(r['n_signatures_named_union']):>3}/{n_testable} "
              f"(best single setting {r['best_single_setting']:.3f})")
        if never and r["arm"] in LEIDEN_ARMS:
            print(f"        never named at any resolution: {never}")

    if "drvi_axes" in ref and any(a in ref for a in LEIDEN_ARMS):
        axis_hit = set(rec.index[rec["named__drvi_axes"].astype(bool)])
        for arm in [a for a in LEIDEN_ARMS if a in ref]:
            cl_hit = set(rec.index[rec[f"named__{arm}"].astype(bool)])
            only_axis = sorted(axis_hit - cl_hit)
            only_cl = sorted(cl_hit - axis_hit)
            print(f"\n  {ARM_LABEL.get(arm, arm)} at resolution {args.ref_resolution}:")
            print(f"    named by the decoder and NOT by the clusters ({len(only_axis)}): "
                  f"{', '.join(only_axis) or '-'}")
            print(f"    named by the clusters and NOT by the decoder ({len(only_cl)}): "
                  f"{', '.join(only_cl) or '-'}")
            col = f"eta2_cluster__{arm}"
            if col in rec.columns and rec[col].notna().any():
                miss = rec.loc[[s for s in only_axis if s in rec.index], col].dropna()
                hit = rec.loc[[s for s in sorted(cl_hit) if s in rec.index], col].dropna()
                if len(miss) and len(hit):
                    print(f"    median eta^2(cluster): {miss.median():.4f} for the ones it "
                          f"could not name, {hit.median():.4f} for the ones it named")

    print("\ndone.")


if __name__ == "__main__":
    main()
