#!/usr/bin/env python3
"""05_9: is a signature PRESENT in a coordinate system? Score-first, axis-free.

THE QUESTION, AND WHY IT IS NOT 05_6 ON ANOTHER SPACE
-----------------------------------------------------
05_3b produced a Harmony embedding as a control on DRVI, and its README calls it a
deliberate dead end for one reason: 05_7/05_8 read an additive DECODER, Harmony has none,
and its dimensions are corrected principal components whose loadings belong to the PCA
*before* the correction. That argument is about Route B and it is correct.

It says nothing about the cells. A signature score is a number per CELL, computed from the
gene expression matrix by 05_6, and it has never seen a latent space - which is exactly why
`signature_common.scores_csv()` names it after the object and not after the run. So the
question "are the programmes this project cares about present in the Harmony space?" can be
asked without asking Harmony for anything it does not offer:

    take the per-cell score. Ask whether cells that are NEIGHBOURS IN THE SPACE have
    similar scores.

That is the whole instrument. It never reads a dimension's gene list, never names an axis,
and never assumes the programme lives on one coordinate. A method whose axes are not
individually interpretable can still carry a programme perfectly, spread over all of them;
this step is built so that such a method scores well rather than being punished for its
architecture.

WHAT IS MEASURED
----------------
Four numbers per (space, readout), all on the WITHIN-COHORT z score, because absolute
scores are not comparable across patients (05_6's argument, inherited):

  morans_i              spatial autocorrelation of the score on the k-NN graph of the space.
                        "Do neighbours agree?" One number, no axis, no direction.
  knn_r2_cell           out-of-sample R^2 of predicting a cell's score from the mean score
                        of its k nearest TRAINING neighbours, 5-fold. "Does the space know
                        the score at all?"
  knn_r2_cohort         the same with the folds grouped by cohort, i.e. the neighbours are
                        never from the held-out patients. "Does it transfer to a patient the
                        fold never saw?" This is the honest one, and it is the harder one.
  top_decile_fold       of the 15 neighbours of a top-decile cell, how many are top-decile
                        themselves, over the 0.10 baseline. "Is the high end a REGION?"

plus three that exist to be read against each other:

  ridge_r2_cohort       the LINEAR share: ridge on the coordinates, the same cohort folds
  nonlinear_gap         knn_r2_cohort - ridge_r2_cohort. POSITIVE means the programme sits
                        in the space in a way no linear function of the coordinates
                        reproduces; NEGATIVE means it is a linear gradient across the whole
                        space and averaging 15 neighbours is simply the worse way to read
                        it. On a corrected PCA the second is the expected sign and is not a
                        defect of anything
  space_r2_cohort       max(knn, ridge). THE COLUMN TO QUOTE for "is it present": the
                        question is the space, not which of two readers was used on it, and
                        a space that holds a programme as a clean gradient must not score
                        low for being easy

and exactly one that is the DRVI claim and is labelled as such:

  max_abs_rho, best_dim the strongest single-dimension Spearman. It is DESCRIPTIVE and it is
                        the one column that is unfair to Harmony by construction, which is
                        why it is one column and not the table. Harmony is the floor here,
                        never the thing being ranked; see 04_9's header for the full form of
                        that argument.

THE NULLS, WHICH ARE THE POINT
------------------------------
Every one of the numbers above is positive for almost any vector on almost any graph, so
none of them means anything alone. Three reference levels are computed on the same graph,
the same folds and the same code path:

  1. PERMUTATION WITHIN COHORT (`--n-perm`, default 20). The score shuffled among the cells
     of its own patient. Keeps the per-cohort mean, destroys the cell-level structure. This
     is the floor the graph itself induces, cohort composition included.
  2. SIZE- AND EXPRESSION-MATCHED RANDOM GENE SETS (`--n-random`, default 5 per signature).
     For each signature, random gene sets of the SAME size drawn from the SAME expression
     bins `sc.tl.score_genes` uses, scored by the same call and standardised the same way.
     This is the level a gene set of that size reaches on this object for no biological
     reason, and it is the number "present" has to beat. Cached; `--n-random 0` skips it.
  3. THE CONFOUNDER ARM. `n_genes_by_counts` and `pct_counts_mt`, z-scored within cohort,
     carried through as readouts. 05_3b found that Harmony CONCENTRATES depth rather than
     removing it (H2, H4), so a signature that does not clear the depth row in the Harmony
     space is a depth readout in the Harmony space, whatever its name is.

THREE SPACES, ONE SET OF CELLS
------------------------------
  harmony   05_3b's corrected 64 components - what this step was asked about
  drvi      05_3's latent space, the phase's own, as the reference level
  pca       the UNCORRECTED PCA Harmony started from (`pca_<harmony_run>.npy`), which is the
            arm that says whether the correction cost the programme anything

The cells and the scores are IDENTICAL in all three - the spaces are realigned on cell name,
never on position, and the step aborts if any of them is missing a cell. The k-NN graph is
built the same way 05_3 and 05_3b build theirs (k = 15, the whole space as the
representation, no rescaling of the coordinates), so the only thing that differs between the
arms is the space.

WHAT THIS STEP DOES NOT DO
--------------------------
It does not name a Harmony dimension, rank Harmony against DRVI as integration methods (that
is phase 02, on what they both promise), or write anything into the per-cell tables 05_6
owns. It reads 05_6's cached scores and never re-scores the real signatures.

Outputs (`<run>` is the DRVI reference run id, which names the CELLS and the SCORES; the
`space` column says which coordinate system each row was measured in):

    tables/<coll>/<run>/signature_presence_<coll>_<run>.csv      one row per (space, readout)
    tables/<coll>/<run>/space_summary_<coll>_<run>.csv           one row per space
    tables/<coll>/<run>/random_null_<coll>_<run>.csv             per random set, if computed
    figures/05_9_embedding_control/<coll>/<run>/
        morans_i_by_space_*.png          the headline, with both nulls drawn on it
        knn_r2_cohort_by_space_*.png     the cross-patient version of the same
        presence_over_null_*.png         observed minus random-set null, in null sd
        knn_vs_ridge_*.png               how much of it is linear
        umap_top_signatures_*.png        the eyeball check, harmony beside drvi

    $DATA_DIR/05_tum/random_signature_scores_<coll>_<run>.csv    the matched null, cached

Usage:
    export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
    cd 05_drvi_tumoral_epi/05_9_embedding_control

    N_LATENT=64 conda run -n benchmark-py-r python3 signature_presence_tum.py --collection gavish
    N_LATENT=64 python3 signature_presence_tum.py                      # scie, the default
    N_LATENT=64 python3 signature_presence_tum.py --collection emt
    N_LATENT=64 python3 signature_presence_tum.py --spaces harmony drvi # without the PCA arm
    N_LATENT=64 python3 signature_presence_tum.py --n-random 0          # permutation null only
    N_LATENT=64 python3 signature_presence_tum.py --overwrite           # re-draw the random sets

`CELL_SET`, `HVG_SET` and `N_LATENT` select the run exactly as in 05_4 - 05_8. `N_LATENT=64`
is not optional in practice: 05_3b only ever ran at 64, and at any other value the harmony
embedding this step needs is not on disk.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
import zlib
from typing import NamedTuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402
import pandas as pd                      # noqa: E402
import scanpy as sc                      # noqa: E402
import scipy.sparse as sp                # noqa: E402
from scipy.stats import rankdata         # noqa: E402
from sklearn.linear_model import RidgeCV                      # noqa: E402
from sklearn.metrics import r2_score                          # noqa: E402
from sklearn.model_selection import GroupKFold, KFold         # noqa: E402
from sklearn.neighbors import NearestNeighbors                # noqa: E402

UTILS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils")
sys.path.insert(0, UTILS_DIR)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "05_2_subsetting"))
import cell_set as CS                    # noqa: E402
import signature_common as C             # noqa: E402
import sig_collections as SC             # noqa: E402

STEP = "05_9_embedding_control"

# The neighbourhood. 15 is what 05_3 and 05_3b both pass to sc.pp.neighbors on their own
# space, and it is repeated here rather than re-derived so that "neighbours in the space"
# means the same thing in this table as it does in the UMAPs those two steps drew.
K_NEIGHBOURS = 15

# The folds. Five of each, which is the usual trade between bias and how long the ridge
# takes; the cohort scheme groups whole patients, so five folds means ~4 patients held out
# together on 19 cohorts.
N_FOLDS = 5

N_BINS = 25              # sc.tl.score_genes' own default, as 05_6 states it
TOP_Q = 0.90             # the decile the region check is about

# The confounder arm. Both are obs columns of every embedding this phase writes, and both are
# carried through the identical code path as readouts - a signature that does not beat them
# in a space is measuring them in that space. `n_genes_by_counts` is the one 05_3b found
# Harmony concentrates rather than removes.
REFERENCE_READOUTS = {
    "__depth": "n_genes_by_counts",
    "__mito": "pct_counts_mt",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    SC.add_argument(p)
    p.add_argument("--spaces", nargs="+", default=["harmony", "drvi", "pca"],
                   choices=["harmony", "drvi", "pca"],
                   help="which coordinate systems to measure in (default: all three). "
                        "'pca' is the uncorrected space Harmony started from")
    p.add_argument("--k", type=int, default=K_NEIGHBOURS,
                   help=f"neighbours in the graph (default {K_NEIGHBOURS}, as 05_3/05_3b)")
    p.add_argument("--n-perm", type=int, default=20,
                   help="within-cohort permutations for the graph null (default 20, 0 skips)")
    p.add_argument("--n-random", type=int, default=5,
                   help="matched random gene sets PER SIGNATURE (default 5, 0 skips). Needs "
                        "the all-genes object and is the slow part of the step")
    p.add_argument("--no-multivariate", action="store_true",
                   help="skip the ridge fits, i.e. leave the linear share unmeasured")
    p.add_argument("--overwrite", action="store_true",
                   help="redraw and re-score the random gene sets instead of reusing the csv")
    return p.parse_args()


# --------------------------------------------------------------------------- #
# The spaces
# --------------------------------------------------------------------------- #

class Space(NamedTuple):
    """One coordinate system, with everything the step needs to read it."""

    coords: pd.DataFrame          # cells x dimensions, INDEXED BY CELL NAME
    prefix: str                   # 'H', 'DR', 'PC' - only ever used to spell `best_dim`
    description: str
    obs: pd.DataFrame | None      # the metadata, when this arm carries it
    umap: np.ndarray | None       # the space's own UMAP, for the eyeball figure


def load_spaces(names: list[str], harmony_run: str, drvi_run: str) -> dict[str, Space]:
    """The requested arms, each indexed BY CELL NAME.

    The realignment onto a common cell order happens once, in `main`, against the index of
    the score table - never by position, which is the failure 05_3b's own DRVI comparison
    guards against for the same reason.
    """
    tum = CS.tum_dir()
    out: dict[str, Space] = {}
    harmony_cells = None

    if "harmony" in names or "pca" in names:
        path = tum / f"embed_{harmony_run}.h5ad"
        if not path.exists():
            sys.exit(f"missing {path}: run 05_3b_harmony_run/run_harmony_tum.py first")
        e = sc.read_h5ad(path)
        harmony_cells = e.obs_names.astype(str)
        if "harmony" in names:
            X = pd.DataFrame(np.asarray(e.X), index=harmony_cells)
            out["harmony"] = Space(X, "H",
                                   f"05_3b: Harmony, {X.shape[1]} corrected components",
                                   e.obs.copy(), np.asarray(e.obsm["X_umap"]))

    if "drvi" in names:
        path = tum / f"embed_{drvi_run}.h5ad"
        if not path.exists():
            sys.exit(f"missing {path}: run 05_3_drvi_run first")
        e = sc.read_h5ad(path)
        X = pd.DataFrame(np.asarray(e.X), index=e.obs_names.astype(str))
        out["drvi"] = Space(X, "DR", f"05_3: DRVI, {X.shape[1]} latent dimensions",
                            e.obs.copy(), np.asarray(e.obsm["X_umap"]))

    if "pca" in names:
        # The UNCORRECTED PCA, written by run_harmony_tum.py as a plain array in the row
        # order of its own embedding - it has no cell names of its own, which is why it takes
        # that embedding's index. If the two disagree in length the arm is DROPPED rather
        # than lined up on a guess: an off-by-one here would compare the right scores against
        # the wrong cells and say nothing about it.
        npy = tum / f"pca_{harmony_run}.npy"
        if not npy.exists():
            print(f"[skip] {npy.name} is not on disk, so the uncorrected arm is dropped",
                  flush=True)
        else:
            raw = np.load(npy)
            if harmony_cells is None or raw.shape[0] != len(harmony_cells):
                print(f"[skip] {npy.name} has {raw.shape[0]} rows against "
                      f"{0 if harmony_cells is None else len(harmony_cells)} cells in the "
                      f"harmony embedding; the uncorrected arm is dropped", flush=True)
            else:
                X = pd.DataFrame(raw, index=harmony_cells)
                out["pca"] = Space(X, "PC",
                                   f"the same {X.shape[1]} components BEFORE the correction",
                                   None, None)

    # Requested order, so the figures put the arms in the order the flag named them.
    return {n: out[n] for n in names if n in out}


# --------------------------------------------------------------------------- #
# The graph, the folds, and the cheap predictors built on them
# --------------------------------------------------------------------------- #

class SpaceGraph:
    """Everything about one space that does not depend on WHICH readout is measured.

    The whole reason this is an object: the k-NN searches are the expensive part and they are
    identical for every readout and every permutation. Built once per space, then a readout
    - real, permuted or random - costs a gather and a mean.
    """

    def __init__(self, X: np.ndarray, cohorts: np.ndarray, k: int, seed: int,
                 n_folds: int = N_FOLDS):
        self.X = np.ascontiguousarray(X, dtype=np.float32)
        self.k = k
        n = self.X.shape[0]

        t0 = time.time()
        nn = NearestNeighbors(n_neighbors=k + 1).fit(self.X)
        _, idx = nn.kneighbors(self.X)
        self.nbr = idx[:, 1:]                       # drop self, (n, k)

        # The binary symmetric connectivity Moran's I is computed on. Binary and not the
        # UMAP-style weighted graph on purpose: a weighting scheme is a second thing that
        # could differ between the arms, and there is nothing to gain from it here.
        rows = np.repeat(np.arange(n), k)
        W = sp.coo_matrix((np.ones(n * k, dtype=np.float32),
                           (rows, self.nbr.ravel())), shape=(n, n)).tocsr()
        self.W = ((W + W.T) > 0).astype(np.float32)

        # Per fold scheme, the neighbours of every TEST cell among the TRAIN cells only.
        self.folds = {}
        for scheme, splitter, groups in (
            ("cell", KFold(n_splits=n_folds, shuffle=True, random_state=seed), None),
            ("cohort", GroupKFold(n_splits=min(n_folds, len(np.unique(cohorts)))), cohorts),
        ):
            parts = []
            for tr, te in splitter.split(self.X, groups=groups):
                nn_tr = NearestNeighbors(n_neighbors=k).fit(self.X[tr])
                _, nb = nn_tr.kneighbors(self.X[te])
                parts.append((tr, te, tr[nb]))       # neighbours as GLOBAL indices
            self.folds[scheme] = parts
        print(f"    graph + folds in {time.time() - t0:.0f}s", flush=True)

        # Rank-transformed coordinates, so a Spearman against a readout is one matrix
        # product rather than 64 scipy calls per readout.
        R = np.apply_along_axis(rankdata, 0, self.X.astype(np.float64))
        self.Rz = (R - R.mean(0)) / np.where(R.std(0) > 0, R.std(0), 1.0)

    # ---- the four measurements, each taking a (n,) or (n_readouts, n) array

    def morans_i(self, vals: np.ndarray) -> np.ndarray:
        """scanpy's own implementation, on the binary graph above. Vectorised over rows."""
        v = np.atleast_2d(np.asarray(vals, dtype=np.float64))
        return np.atleast_1d(sc.metrics.morans_i(self.W, v))

    def knn_r2(self, y: np.ndarray, scheme: str) -> float:
        y = np.asarray(y, dtype=np.float64)
        pred = np.empty_like(y)
        pred[:] = np.nan
        for _tr, te, nb in self.folds[scheme]:
            pred[te] = y[nb].mean(axis=1)
        return float(r2_score(y, pred))

    def ridge_r2(self, y: np.ndarray, scheme: str, seed: int) -> float:
        y = np.asarray(y, dtype=np.float64)
        pred = np.empty_like(y)
        pred[:] = np.nan
        for tr, te, _nb in self.folds[scheme]:
            m = RidgeCV(alphas=np.logspace(-2, 4, 13)).fit(self.X[tr], y[tr])
            pred[te] = m.predict(self.X[te])
        return float(r2_score(y, pred))

    def top_decile_fold(self, y: np.ndarray, q: float = TOP_Q) -> tuple[float, float]:
        """(fraction of a top cell's neighbours that are top, that over the baseline)."""
        y = np.asarray(y, dtype=np.float64)
        cut = np.quantile(y, q)
        top = y >= cut
        base = float(top.mean())
        if base <= 0 or not top.any():
            return float("nan"), float("nan")
        frac = float(top[self.nbr[top]].mean())
        return frac, frac / base

    def best_dim(self, y: np.ndarray) -> tuple[float, int]:
        """(max |Spearman| over the dimensions, which one). Descriptive; see the header."""
        ry = rankdata(np.asarray(y, dtype=np.float64))
        ry = (ry - ry.mean()) / (ry.std() if ry.std() > 0 else 1.0)
        rho = (self.Rz * ry[:, None]).mean(axis=0)
        j = int(np.nanargmax(np.abs(rho)))
        return float(rho[j]), j


# --------------------------------------------------------------------------- #
# Standardisation, reused verbatim from 05_6's rule
# --------------------------------------------------------------------------- #

def z_within(df: pd.DataFrame, cohorts: pd.Series) -> pd.DataFrame:
    """Z-score each column within `cohort`, 05_6's GROUPBY and 05_6's zero-variance rule."""
    g = df.groupby(cohorts.values, observed=True)
    z = g.transform(lambda s: (s - s.mean()) / s.std(ddof=0) if s.std(ddof=0) > 0 else s * 0.0)
    return z.fillna(0.0)


def permute_within(values: np.ndarray, codes: np.ndarray, rng) -> np.ndarray:
    """Shuffle `values` among the cells of each cohort. The graph null of the header."""
    out = values.copy()
    for c in np.unique(codes):
        m = codes == c
        out[m] = rng.permutation(values[m])
    return out


# --------------------------------------------------------------------------- #
# The matched random gene sets
# --------------------------------------------------------------------------- #

def expression_bins(adata, n_bins: int = N_BINS) -> pd.Series:
    """Every gene's expression bin, built the way sc.tl.score_genes builds its control bins.

    Rank of the mean expression, cut into `n_bins` equal-frequency-ish bins. Reproducing it
    here rather than importing it is deliberate: the null has to be drawn from the SAME bins
    the real score's control set came from, and scanpy exposes that binning nowhere.
    """
    X = adata.X
    mean = np.asarray(X.mean(axis=0)).ravel()
    ranked = pd.Series(rankdata(mean, method="min"), index=adata.var_names)
    return pd.Series(pd.cut(ranked, n_bins, labels=False), index=adata.var_names)


def draw_matched_sets(genes: list[str], bins: pd.Series, n_sets: int, rng) -> list[list[str]]:
    """`n_sets` gene sets of the same size as `genes`, one gene per bin of each real gene.

    A random set of the same SIZE is not enough: a signature of highly expressed genes scores
    differently from a signature of rare ones, whatever the biology. Matching bin by bin is
    what makes the null the level of "a gene set like this one, with no meaning".
    """
    by_bin = {b: np.asarray(idx) for b, idx in bins.groupby(bins).groups.items()}
    present = [g for g in genes if g in bins.index]
    wanted = bins.loc[present].values
    banned = set(present)
    sets = []
    for _ in range(n_sets):
        pick = []
        for b in wanted:
            pool = by_bin.get(b)
            if pool is None or len(pool) == 0:
                continue
            for _try in range(20):
                g = str(rng.choice(pool))
                if g not in banned and g not in pick:
                    pick.append(g)
                    break
        if len(pick) >= C.MIN_SIGNATURE_GENES:
            sets.append(pick)
    return sets


def build_random_scores(coll, sets: dict[str, list[str]], obs_names, cohorts,
                        n_per_sig: int, cache, overwrite: bool, seed: int) -> pd.DataFrame:
    """Per-cell z scores of the matched random sets. Cached: it is the slow part of the step.

    Columns are `<SIGNATURE>__rnd<i>`, so a null column always says which real signature's
    size and expression profile it was matched to.
    """
    if cache.exists() and not overwrite:
        print(f"[have] {cache.name}, reading it back", flush=True)
        df = pd.read_csv(cache, index_col=0)
        df.index = df.index.astype(str)
        return df.reindex(obs_names)

    full = CS.path(".h5ad")
    if not full.exists():
        print(f"[skip] {full.name} is not on disk, so the matched random null cannot be "
              f"built; the permutation null still is", flush=True)
        return pd.DataFrame(index=obs_names)

    print(f"[read] {full}  (the all-genes object, for the null only)", flush=True)
    adata = sc.read_h5ad(full)
    adata = adata[obs_names].copy()
    bins = expression_bins(adata)
    rng = np.random.default_rng(seed)

    raw = {}
    for name, genes in sets.items():
        if len(genes) < C.MIN_SIGNATURE_GENES:
            continue
        for i, pick in enumerate(draw_matched_sets(genes, bins, n_per_sig, rng)):
            key = f"{name}__rnd{i}"
            sc.tl.score_genes(adata, gene_list=pick, ctrl_size=len(pick), n_bins=N_BINS,
                              random_state=seed, use_raw=False, score_name=key)
            raw[key] = adata.obs[key].astype(float).values
        print(f"  {name:28s} {len(genes):4d} genes -> {n_per_sig} matched random sets",
              flush=True)

    if not raw:
        return pd.DataFrame(index=obs_names)
    df = z_within(pd.DataFrame(raw, index=adata.obs_names.astype(str)), cohorts)
    df.to_csv(cache)
    print(f"[write] {cache}  ({df.shape[0]} x {df.shape[1]})", flush=True)
    return df.reindex(obs_names)


# --------------------------------------------------------------------------- #
# Tables
# --------------------------------------------------------------------------- #

HEADER_NOTE = ("three coordinate systems on ONE set of cells and ONE set of per-cell scores; "
               "the `space` column says which. The run id names the cells and the scores "
               "(05_6), not the space")


def write_table(df: pd.DataFrame, name: str, coll, run_id: str, index: bool = False):
    """C.table_path for the path, a header of this step's own for the first lines.

    Not `C.write_table`: that one stamps "<Space> run <run_id>" from the EMBEDDINGS registry,
    which is exactly right for a table measured in one space and exactly wrong for this one.
    """
    path = C.table_path(name, coll, run_id)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"# {name} | collection {coll.name} ({coll.title}) | run {run_id} "
                 f"| {HEADER_NOTE} | 05_9 embedding control\n")
        for line in C.CAVEAT.split(". "):
            if line.strip():
                fh.write(f"# CAVEAT: {line.strip().rstrip('.')}.\n")
        df.to_csv(fh, index=index)
    print(f"[table] {path}  ({df.shape[0]} x {df.shape[1]})", flush=True)
    return path


# --------------------------------------------------------------------------- #
# Figures
# --------------------------------------------------------------------------- #
#
# The space is the only categorical encoding in this step, it has three levels and they are
# assigned in a fixed order that does not move when an arm is dropped - `harmony` keeps its
# colour whether or not `pca` was run. Everything else is magnitude (bars) or a signed score
# on a UMAP (a diverging map with a neutral midpoint), and no figure has two y-scales.

SPACE_COLOUR = {"harmony": "#2f6f9f", "drvi": "#c25e00", "pca": "#7a7a7a"}
SPACE_LABEL = {"harmony": "Harmony (05_3b)", "drvi": "DRVI (05_3)", "pca": "uncorrected PCA"}
NULL_COLOUR = "#b0b0b0"
INK = "#222222"


def _readout_order(coll, table: pd.DataFrame) -> list[str]:
    """Collection order for the signatures, the confounder arm last and marked.

    `Collection.order` keeps only the names it has an axis for, so anything it does not know
    is appended rather than dropped - a readout missing from a figure because a registry
    entry was forgotten is the kind of silence this phase is built against.
    """
    seen = list(dict.fromkeys(table["readout"].tolist()))
    sigs = [r for r in seen if not r.startswith("__")]
    refs = [r for r in seen if r.startswith("__")]
    ordered = coll.order(sigs, for_figure=True)
    return ordered + [r for r in sigs if r not in ordered] + sorted(refs)


def _bar_panel(ax, table, order, column, spaces, null_column=None, xlabel=""):
    """Grouped horizontal bars, one group per readout and one bar per space.

    The confounder rows are drawn hollow and hatched in their space's own colour: `kind` is
    the one distinction a reader must not miss - a signature below the depth bar of its own
    space is a depth readout - and it is therefore carried by texture as well as by the row
    label, never by colour, which is spent entirely on the space.
    """
    y = np.arange(len(order))
    h = 0.8 / max(len(spaces), 1)
    is_ref = np.array([r.startswith("__") for r in order])
    for si, s in enumerate(spaces):
        sub = table[table["space"] == s].set_index("readout").reindex(order)
        pos = y + si * h - 0.4 + h / 2
        v = sub[column].values
        ax.barh(pos[~is_ref], v[~is_ref], height=h * 0.86, color=SPACE_COLOUR[s],
                label=SPACE_LABEL[s], linewidth=0)
        ax.barh(pos[is_ref], v[is_ref], height=h * 0.86, color="none",
                edgecolor=SPACE_COLOUR[s], linewidth=1.1, hatch="////")
        if null_column and null_column in sub:
            # The matched random-set level for THIS space, as a tick on the bar it qualifies.
            ax.scatter(sub[null_column].values, pos,
                       marker="|", s=90, color=NULL_COLOUR, zorder=3, linewidths=1.4)
    ax.set_yticks(y)
    ax.set_yticklabels([r.replace("__", "") + ("  (confounder)" if r.startswith("__") else "")
                        for r in order], fontsize=8)
    ax.invert_yaxis()
    ax.axvline(0, color=INK, linewidth=0.8)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.grid(axis="x", color="0.9", linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)


def fig_bars(table, coll, run_id, spaces, column, null_column, title, xlabel, name):
    order = _readout_order(coll, table)
    fig, ax = plt.subplots(figsize=(8.5, C.fig_span(len(order), 0.34, 2.2)))
    _bar_panel(ax, table, order, column, spaces, null_column, xlabel)
    handles, labels = ax.get_legend_handles_labels()
    handles.append(plt.Rectangle((0, 0), 1, 1, facecolor="none", edgecolor=INK,
                                 linewidth=1.1, hatch="////"))
    labels.append("confounder, not a signature")
    if null_column and null_column in table:
        handles.append(plt.Line2D([], [], color=NULL_COLOUR, marker="|", linestyle="none",
                                  markersize=9, markeredgewidth=1.4))
        labels.append("matched random gene sets")
    ax.legend(handles, labels, fontsize=8, frameon=False, ncol=2,
              loc="upper center", bbox_to_anchor=(0.5, -0.045 - 2.2 / len(order)))
    ax.set_title(title, fontsize=10, color=INK)
    fig.tight_layout()
    C.savefig(name, STEP, coll, fig=fig, run_id=run_id)
    plt.close(fig)


def fig_over_null(table, coll, run_id, spaces):
    """Observed Moran's I minus the matched-random mean, in units of the random sd.

    The one figure that answers the question as asked: a signature is PRESENT in a space if
    it is further from what a gene set of its size reaches there by accident than the spread
    of those accidents. Two vertical guides, at 0 (indistinguishable from a random set of the
    same size) and at 3 sd.
    """
    if "morans_i_vs_random_z" not in table:
        return
    order = _readout_order(coll, table)
    fig, ax = plt.subplots(figsize=(8.5, C.fig_span(len(order), 0.34, 2.2)))
    _bar_panel(ax, table, order, "morans_i_vs_random_z", spaces,
               xlabel="(Moran's I  -  matched random mean) / random sd")
    for x, ls in ((0.0, "-"), (3.0, ":")):
        ax.axvline(x, color=INK if x == 0 else "0.5", linewidth=0.8, linestyle=ls)
    ax.legend(fontsize=8, frameon=False, ncol=3, loc="upper center",
              bbox_to_anchor=(0.5, -0.045 - 2.2 / len(order)))
    ax.set_title("how far above a size-matched random gene set, per space", fontsize=10,
                 color=INK)
    fig.tight_layout()
    C.savefig("presence_over_null", STEP, coll, fig=fig, run_id=run_id)
    plt.close(fig)


def fig_knn_vs_ridge(table, coll, run_id, spaces):
    """The non-linearity gap. One panel per space, identical axes, the diagonal drawn."""
    cols = [c for c in ("ridge_r2_cohort", "knn_r2_cohort") if c in table]
    if len(cols) < 2 or table["ridge_r2_cohort"].isna().all():
        return
    fig, axes = plt.subplots(1, len(spaces), figsize=(3.4 * len(spaces), 3.6), sharex=True,
                             sharey=True)
    axes = np.atleast_1d(axes)
    lim = [min(0.0, float(table[cols].min().min())) - 0.03,
           float(table[cols].max().max()) + 0.05]
    for ax, s in zip(axes, spaces):
        sub = table[table["space"] == s]
        sig = sub[~sub["readout"].str.startswith("__")]
        ref = sub[sub["readout"].str.startswith("__")]
        ax.plot(lim, lim, color="0.75", linewidth=0.8, zorder=0)
        ax.scatter(sig["ridge_r2_cohort"], sig["knn_r2_cohort"], s=26,
                   color=SPACE_COLOUR[s], edgecolor="white", linewidth=0.6, label="signature")
        ax.scatter(ref["ridge_r2_cohort"], ref["knn_r2_cohort"], s=44, marker="D",
                   color="white", edgecolor=INK, linewidth=1.0, label="confounder")
        for _, r in ref.iterrows():
            ax.annotate(r["readout"].replace("__", ""),
                        (r["ridge_r2_cohort"], r["knn_r2_cohort"]),
                        textcoords="offset points", xytext=(6, -2), fontsize=7, color=INK)
        ax.set_title(SPACE_LABEL[s], fontsize=9, color=INK)
        ax.set_xlabel("ridge R$^2$ (linear in the coordinates)", fontsize=8)
        ax.grid(color="0.92", linewidth=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axes[0].set_ylabel("k-NN R$^2$ (neighbours in the space)", fontsize=8)
    axes[0].set_xlim(lim)
    axes[0].set_ylim(lim)
    axes[-1].legend(fontsize=7, frameon=False, loc="lower right")
    fig.suptitle("above the diagonal: neighbourhood structure a linear fit misses.  below: "
                 "the programme is a linear gradient, and a local average is the worse "
                 "reader of it  (cohort-held-out folds)", fontsize=8, color=INK)
    fig.tight_layout()
    C.savefig("knn_vs_ridge", STEP, coll, fig=fig, run_id=run_id)
    plt.close(fig)


def fig_umaps(table, coll, run_id, spaces, umaps, z, n_top=6):
    """The eyeball check: the best-ranked signatures on each space's own UMAP.

    Ranked by the Harmony k-NN R^2 when Harmony is in the run, because it is the space the
    step was asked about. Diverging map centred at zero: the score is a signed z.
    """
    have = [s for s in spaces if umaps.get(s) is not None]
    if not have:
        return
    rank_space = "harmony" if "harmony" in have else have[0]
    sub = table[(table["space"] == rank_space) & (~table["readout"].str.startswith("__"))]
    top = sub.sort_values("space_r2_cohort", ascending=False)["readout"].head(n_top).tolist()
    if not top:
        return
    fig, axes = plt.subplots(len(top), len(have),
                             figsize=(3.1 * len(have), 2.9 * len(top)), squeeze=False)
    vmax = float(np.nanpercentile(np.abs(z[top].values), 99))
    for r, readout in enumerate(top):
        v = z[readout].values
        o = np.argsort(np.abs(v))                 # the extremes drawn last, not hidden
        for c, s in enumerate(have):
            ax = axes[r][c]
            U = umaps[s]
            im = ax.scatter(U[o, 0], U[o, 1], c=v[o], s=1.2, cmap="RdBu_r",
                            vmin=-vmax, vmax=vmax, linewidths=0, rasterized=True)
            ax.set_xticks([]); ax.set_yticks([])
            for side in ax.spines.values():
                side.set_visible(False)
            if r == 0:
                ax.set_title(SPACE_LABEL[s], fontsize=9, color=INK)
            if c == 0:
                ax.set_ylabel(readout, fontsize=8, color=INK)
        cax = fig.colorbar(im, ax=axes[r], fraction=0.02, pad=0.01)
        cax.set_label("z within cohort", fontsize=7)
        cax.ax.tick_params(labelsize=6)
    fig.suptitle(f"the {len(top)} signatures the {SPACE_LABEL[rank_space]} space predicts "
                 f"best, on each space's own UMAP", fontsize=9, color=INK)
    C.savefig("umap_top_signatures", STEP, coll, fig=fig, run_id=run_id, dpi=200)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

def main():
    args = parse_args()
    coll = SC.resolve(args)

    drvi_run = CS.run_id(method="drvi")
    harmony_run = CS.run_id(method="harmony")

    CS.banner(f"05_9 - is a signature present in a space? ({coll.name})")
    print(f"collection : {coll.title}")
    print(f"question   : {coll.question}")
    print(f"scores     : {C.scores_csv(coll).name}  (05_6, embedding-independent)")
    print(f"reference  : {drvi_run}   |   harmony: {harmony_run}")
    print(f"k = {args.k}, {N_FOLDS} folds, {args.n_perm} permutations, "
          f"{args.n_random} matched random sets per signature\n", flush=True)

    # ---- the cells, and the scores that were computed on them ----------------
    scores_csv = C.scores_csv(coll)
    if not scores_csv.exists():
        sys.exit(f"missing {scores_csv}: run 05_6_cell_first/cell_first_tum.py "
                 f"--collection {args.collection} first")
    scores = pd.read_csv(scores_csv, index_col=0)
    scores.index = scores.index.astype(str)
    z_cols = [c for c in scores.columns if c.startswith("z_")]
    if not z_cols:
        sys.exit(f"{scores_csv.name} has no z_ columns; it is not a 05_6 score table")
    z = scores[z_cols].rename(columns=lambda c: c[2:])
    print(f"[read] {scores_csv.name}: {z.shape[0]} cells x {z.shape[1]} readouts", flush=True)

    # ---- the spaces, realigned on cell NAME ---------------------------------
    spaces_raw = load_spaces(args.spaces, harmony_run, drvi_run)
    if not spaces_raw:
        sys.exit("no coordinate system could be read; nothing to measure")

    cells = z.index
    for name, s in spaces_raw.items():
        missing = cells.difference(s.coords.index)
        if len(missing):
            sys.exit(f"{name}: {len(missing)} of the scored cells are not in its embedding "
                     f"(e.g. {list(missing[:3])}). Two spaces compared on two different cell "
                     f"sets is not a comparison")

    # The metadata comes from whichever embedding carries it; they are the same obs.
    obs = next((s.obs for s in spaces_raw.values() if s.obs is not None), None)
    if obs is None:
        sys.exit("no embedding with an obs table was read; the cohort strata are unknown")
    obs = obs.reindex(cells)
    cohorts = obs[C.BATCH_KEY].astype(str)
    print(f"{cohorts.nunique()} cohorts, {len(cells)} cells, "
          f"spaces: {', '.join(spaces_raw)}\n", flush=True)

    # ---- the confounder arm, standardised exactly like the signatures -------
    ref = {}
    for key, col in REFERENCE_READOUTS.items():
        if col in obs:
            ref[key] = obs[col].astype(float).values
    if ref:
        z = pd.concat([z, z_within(pd.DataFrame(ref, index=cells), cohorts)], axis=1)
        print(f"confounder arm: {', '.join(ref)} "
              f"({', '.join(REFERENCE_READOUTS[k] for k in ref)}), z within cohort", flush=True)

    # ---- the matched random gene sets ---------------------------------------
    rnd = pd.DataFrame(index=cells)
    if args.n_random > 0:
        gmt = C.gmt_path(coll)
        if not gmt.exists():
            print(f"[skip] {gmt.name} is not on disk, so the matched random null is skipped; "
                  f"run 05_4_signatures/build_signatures_tum.py for it", flush=True)
        else:
            sets = C.read_gmt(gmt)
            cache = CS.tum_dir() / f"random_signature_scores_{coll.name}_{drvi_run}.csv"
            print("\n>>> matched random gene sets (the null 'present' has to beat)",
                  flush=True)
            rnd = build_random_scores(coll, sets, cells, cohorts, args.n_random, cache,
                                      args.overwrite, C.SEED)

    # ---- measure, space by space -------------------------------------------
    rows, space_rows = [], []
    cohort_codes = pd.Categorical(cohorts).codes

    for name, space in spaces_raw.items():
        print(f"\n>>> {name}: {space.description}", flush=True)
        Xa = space.coords.reindex(cells).values.astype(np.float32)
        g = SpaceGraph(Xa, cohort_codes, args.k, C.SEED)

        # The null levels for this space, computed once over all readouts and permutations.
        rnd_stats = {}
        if rnd.shape[1]:
            rnd_I = g.morans_i(rnd.values.T)
            rnd_knn = {c: g.knn_r2(rnd[c].values, "cohort") for c in rnd.columns}
            for c, I in zip(rnd.columns, rnd_I):
                rnd_stats.setdefault(c.split("__rnd")[0], []).append((I, rnd_knn[c]))

        pr, mean_off = C.effective_rank(Xa)
        space_rows.append({
            "space": name, "description": space.description, "n_dims": Xa.shape[1],
            "n_cells": Xa.shape[0], "k": args.k,
            "effective_rank_of_space": round(pr, 2),
            "mean_abs_r_between_dims": round(mean_off, 4),
        })

        for readout in z.columns:
            y = z[readout].values.astype(np.float64)
            I = float(g.morans_i(y)[0])
            frac, fold = g.top_decile_fold(y)
            rho, j = g.best_dim(y)

            row = {
                "space": name, "readout": readout,
                "axis": ("confounder" if readout.startswith("__")
                         else coll.axis_of.get(readout, "")),
                "kind": "reference" if readout.startswith("__") else "signature",
                "morans_i": I,
                "knn_r2_cell": g.knn_r2(y, "cell"),
                "knn_r2_cohort": g.knn_r2(y, "cohort"),
                "top_decile_nbr_frac": frac,
                "top_decile_fold": fold,
                "max_abs_rho": abs(rho),
                "best_dim": f"{space.prefix} {j + 1}",
                "best_dim_rho": rho,
            }

            if not args.no_multivariate:
                row["ridge_r2_cohort"] = g.ridge_r2(y, "cohort", C.SEED)
                row["nonlinear_gap"] = row["knn_r2_cohort"] - row["ridge_r2_cohort"]
                # WHAT THE SPACE KNOWS, by whichever of the two readers reads it better.
                # Neither reader is the question - the question is the space - and pinning
                # the answer to the k-NN alone would understate a space that holds the
                # programme as a clean linear gradient, which is exactly what a corrected
                # PCA is built to be. This is the column to quote for "is it present".
                row["space_r2_cohort"] = max(row["knn_r2_cohort"], row["ridge_r2_cohort"])
            else:
                row["space_r2_cohort"] = row["knn_r2_cohort"]

            # 1. the graph null: the same score shuffled inside its own patient.
            # Seeded per (space, readout) rather than from one running stream, so a rerun
            # with a different --spaces or a different collection gives the same null for
            # the readouts it shares with this one.
            if args.n_perm > 0:
                rng = np.random.default_rng(
                    zlib.crc32(f"{C.SEED}|{name}|{readout}".encode()))
                perms = np.vstack([permute_within(y, cohort_codes, rng)
                                   for _ in range(args.n_perm)])
                pI = g.morans_i(perms)
                pknn = np.array([g.knn_r2(p, "cohort") for p in perms])
                pfold = np.array([g.top_decile_fold(p)[1] for p in perms])
                row.update({
                    "morans_i_perm_mean": float(pI.mean()),
                    "morans_i_perm_sd": float(pI.std(ddof=1)) if args.n_perm > 1 else np.nan,
                    "morans_i_perm_p": float((pI >= I).sum() + 1) / (args.n_perm + 1),
                    "knn_r2_cohort_perm_mean": float(pknn.mean()),
                    "top_decile_fold_perm_mean": float(np.nanmean(pfold)),
                })

            # 2. the matched random sets, which is the level that actually means something
            st = rnd_stats.get(readout)
            if st:
                a = np.array(st, dtype=float)
                m = float(a[:, 0].mean())
                s = float(a[:, 0].std(ddof=1)) if len(a) > 1 else float("nan")
                row.update({
                    "morans_i_random_mean": m,
                    "morans_i_random_sd": s,
                    "morans_i_vs_random_z": (I - m) / s if s and s > 0 else np.nan,
                    "knn_r2_cohort_random_mean": float(a[:, 1].mean()),
                    "n_random_sets": int(len(a)),
                })
            rows.append(row)
        print(f"    {len(z.columns)} readouts measured", flush=True)

    table = pd.DataFrame(rows)
    order_all = _readout_order(coll, table)
    table["__o"] = table["readout"].map({r: i for i, r in enumerate(order_all)})
    table = table.sort_values(["space", "__o"]).drop(columns="__o").reset_index(drop=True)

    write_table(table.round(5), "signature_presence", coll, drvi_run)
    write_table(pd.DataFrame(space_rows), "space_summary", coll, drvi_run)
    if rnd.shape[1]:
        write_table(pd.DataFrame({
            "column": rnd.columns,
            "matched_to": [c.split("__rnd")[0] for c in rnd.columns],
            "draw": [int(c.split("__rnd")[1]) for c in rnd.columns],
        }), "random_null", coll, drvi_run)

    # ---- the reading, in the log, before any figure -------------------------
    print("\n" + "=" * 78, flush=True)
    for name in spaces_raw:
        sub = table[(table["space"] == name) & (table["kind"] == "signature")]
        refs = table[(table["space"] == name) & (table["kind"] == "reference")]
        depth = refs[refs["readout"] == "__depth"]["space_r2_cohort"]
        above = (sub["space_r2_cohort"] > float(depth.iloc[0])).sum() if len(depth) else np.nan
        line = (f"{SPACE_LABEL[name]:22s} mean Moran's I {sub['morans_i'].mean():.3f}  "
                f"mean R2 (cohort-held-out, better reader) "
                f"{sub['space_r2_cohort'].mean():.3f}")
        if len(depth):
            line += (f"  |  depth reaches {float(depth.iloc[0]):.3f}, "
                     f"{above}/{len(sub)} signatures above it")
        print(line, flush=True)
    if "morans_i_vs_random_z" in table:
        # Only the readouts that HAVE a matched null are in the denominator: a readout that
        # is not a gene list (CytoTRACE2) has none by construction and counting it as a
        # failure would be a different claim from the one being made.
        tested = table[(table["kind"] == "signature")
                       & table["morans_i_vs_random_z"].notna()]
        cl = tested[tested["morans_i_vs_random_z"] > 3]
        untested = ((table["kind"] == "signature")
                    & table["morans_i_vs_random_z"].isna()).sum()
        print(f"\nsignatures more than 3 random-set sd above the matched null: "
              f"{len(cl)} of {len(tested)} (space, readout) pairs"
              + (f"  ({untested} pairs have no matched null: not a gene list)"
                 if untested else ""), flush=True)
    print("=" * 78 + "\n", flush=True)

    # ---- figures ------------------------------------------------------------
    spaces = list(spaces_raw)
    fig_bars(table, coll, drvi_run, spaces, "morans_i",
             "morans_i_random_mean" if "morans_i_random_mean" in table else None,
             "do neighbours in the space agree on the score?",
             "Moran's I on the k-NN graph", "morans_i_by_space")
    fig_bars(table, coll, drvi_run, spaces, "knn_r2_cohort",
             "knn_r2_cohort_random_mean" if "knn_r2_cohort_random_mean" in table else None,
             "can the space predict the score of a patient it never saw?",
             "out-of-sample k-NN R$^2$, folds grouped by cohort", "knn_r2_cohort_by_space")
    fig_over_null(table, coll, drvi_run, spaces)
    if not args.no_multivariate:
        fig_knn_vs_ridge(table, coll, drvi_run, spaces)
    fig_umaps(table, coll, drvi_run, spaces,
              {n: s.umap for n, s in spaces_raw.items()}, z)

    print(f"[ok] {coll.name}: {len(table)} (space, readout) rows over "
          f"{len(spaces)} spaces", flush=True)


if __name__ == "__main__":
    main()
