"""Shared helpers for the signature-interpretation steps of the malignant phase (05_4 - 05_8).

Lives in the phase's `utils/`, next to the scripts that use it, the same way
`02_integration_benchmark/utils/` holds `h5ad_compat.py` and `metrics_shared.py`. It is
imported by the five step scripts with the idiom those scripts use:

    UTILS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils")
    sys.path.insert(0, UTILS_DIR)
    import signature_common as C

Duplicated from `04_drvi_epithelial/utils/` rather than imported from it, the convention of
every phase in this repo: 05 reads 05_3's outputs and shares no code with 04, so each phase
still reads as a self-contained Materials & Methods section and 04 stays frozen. What is NOT
a copy is listed under "What differs from 04" below - four things, each of them a decision.

What differs from 04, and why:

  1. The object is the malignant subset (36,192 cells x 24,779 genes) and the run is
     `drvi_tum_32`, both resolved through `05_2_subsetting/cell_set.py` so that `CELL_SET`
     works here exactly as it does in 05_2 and 05_3.
  2. `CAVEAT` is new. 04's says no CNV inference has been run, which stays true of 04's
     tables; this phase ran one and has its own limits to state.
  3. `GROUP_KEY` / `LEIDEN_KEY` exist here and not in 04, because `cell_type` is the constant
     `malignant` on this object. Every step that needs to group cells by something reads them
     rather than `cell_type`. See 05_3's README for the whole argument.
  4. `EMBEDDINGS` holds DRVI only. 04 registers Harmony and PCA for its 04_9 embedding
     control; that step is not planned here, so the two entries would point at files nobody
     writes. The machinery is kept intact - adding a control means one entry plus a writer.

Nothing here computes anything on its own; it is the paths, the figure/table writers and the
interpretability-score accessor the five steps all need.

What is NOT here, deliberately: which signatures exist, which axes they sit on, and what the
target region looks like on the cell-first plane. That differs between the two readouts the
phase runs - `scie` (stemness x immunogenicity) and `emt` - and lives in `sig_collections.py`,
so the step scripts stay single files taking `--collection` instead of being duplicated. Every
writer below therefore takes a collection and scopes its output to it: nothing the EMT run
writes can land on a filename the SCIE run owns.
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

# The CELL_SET -> prefix mapping, the thresholds of the subset and the label keys live in
# 05_2_subsetting/cell_set.py and nowhere else, for the reason its own docstring gives: a
# prefix recomputed independently in each script is a silent-overwrite bug waiting to
# happen. 05_3 imports the same module the same way.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "05_2_subsetting"))
import cell_set as CS  # noqa: E402

# --------------------------------------------------------------------------- #
# The run this step reads
# --------------------------------------------------------------------------- #

N_LATENT = 32                       # the 05_3 run of this phase, see drvi_tum.ipynb
# `drvi_tum_32` under the default CELL_SET, `drvi_epicnv_32` on the control set - the same
# id 05_3 writes its model, embedding and downstream object under.
DEFAULT_RUN_ID = f"drvi_{CS.compartment()}_{N_LATENT}"

# The run every step writes under. This is the ONE mutable global of the module: 05_6 rebinds
# it through `set_embedding()` when it is asked to read a coordinate system other than DRVI's,
# so that every table and figure name follows the embedding without each writer having to be
# told about it. 05_7 and 05_8 never call `set_embedding`, so for them nothing moves: they
# read and write `drvi_tum_32` exactly as before. With only DRVI registered today, nothing
# moves for 05_6 either - the mechanism is kept, not exercised.
RUN_ID = DEFAULT_RUN_ID

SEED = 0                            # scoring / subsampling seed, as in 01_4 and 03_1
SCORE_KEY = "OOD_combined"          # 05_3's interpretability scores, as in 04_3 and 03_3

# Vanished dimensions are NOT pruned in this stage. 03_2's `plot_pruned_umap_nonimm.py`
# measured what pruning does to a DRVI space and the answer was nothing - the vanished
# dimensions carry ~1e-05 of the latent variance - so pruning buys no cleanliness, while
# dropping dimensions before the correlations and the ORA silently decides, ahead of the
# analysis, which axes are allowed to mean something. It matters more here than in 04: at
# n_latent 32 there are half as many axes to begin with, so pruning would spend a larger
# share of them. 05_6 - 05_8 therefore run on all
# `N_LATENT` dimensions and all 2 x `N_LATENT` directions; `var['vanished']` and
# `var['vanished_*_direction']` are still read and reported, so any dimension that does
# come out significant can be checked against its flag. Set this to True to get the
# pruned behaviour back, unchanged, everywhere at once.
PRUNE_VANISHED = False

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #

# utils/signature_common.py -> 05_drvi_tumoral_epi/ -> the repo root
PROJECT_DIR = Path(__file__).resolve().parents[2]
os.environ.setdefault("DATA_DIR", str(PROJECT_DIR / "datasets"))   # cell_set reads the env
DATA_DIR = CS.data_dir()

TUM_DIR = CS.tum_dir()                        # $DATA_DIR/05_tum, the heavy objects
# The signature text files are shared with 04: they are the collaborator's lists, one gene
# symbol per line, and are a property of the prior knowledge rather than of a compartment.
# Nothing writes into this directory.
SIG_DIR = Path(os.environ.get("SIGNATURE_DIR", DATA_DIR / "signatures"))

PHASE_DIR = PROJECT_DIR / "05_drvi_tumoral_epi"

# Small result tables, versioned in the repo as appendix material. Phase-level, next to
# `figures/`, and for the same reason: 05_8 reads what 05_6 and 05_7 wrote, so a per-step
# `tables/` would mean steps reaching into each other's folders. 05_1 already writes here.
TABLE_DIR = PHASE_DIR / "tables"

# Inputs, all read-only here. Through `cell_set.path()`, so `CELL_SET=epi` reads the control
# set's objects and writes under its own prefix without a second copy of this module.
FULL_H5AD = CS.path(".h5ad")                            # 05_2: all genes, log-normalised
HVG_H5AD = CS.path("_hvg_2k.h5ad")                      # 05_2: the DRVI training features
EMBED_H5AD = TUM_DIR / f"embed_{RUN_ID}.h5ad"           # 05_3: latent space + scores

# --------------------------------------------------------------------------- #
# What groups the cells, now that `cell_type` cannot
# --------------------------------------------------------------------------- #
#
# `cell_type` is the constant `malignant` on this object: it was written by 05_2 from the
# CNV call, and it is what makes phase 05 different from phase 04. Anything in this stage
# that would have grouped cells by it reads these instead, in this order. The whole
# argument - including why the borrowed label is not simply dropped - is in
# `05_3_drvi_run/README.md`; the short form:
#
#   LEIDEN_KEY   05_2's clustering, computed ON these cells. Leads. Its clusters come from
#                the expression graph; only their NUMBER was picked by maximising NMI
#                against GROUP_KEY (`clustering_tum.py`), so it is not fully independent
#                of it either.
#   GROUP_KEY    the pre-CNV CellTypist label. Obsolete as an identity - that is what this
#                phase establishes - and kept as the only non-constant CellTypist column
#                left on the subset. A landmark, never a claim.
#
# Neither is an input to anything: DRVI was trained with `batch_key='cohort'` and no
# `labels_key`, and Route A standardises within `cohort` alone (see 05_6). They enter only
# where cells have to be summarised into rows of a table.
LABEL_KEY = CS.LABEL_KEY              # 'cell_type', constant under CELL_SET=tum
GROUP_KEY = CS.PRIOR_LABEL_KEY        # 'cell_type_01_4', the pre-CNV CellTypist label
LEIDEN_KEY = "optscib_tum_leiden"     # 05_2's clustering, as named by clustering_tum.py
BATCH_KEY = CS.BATCH_KEY              # 'cohort'


def grouping_keys(obs) -> list[str]:
    """The keys of `obs` that can actually group these cells, best first.

    A constant column is not an error here, it is the shape of this subset, so it is
    dropped rather than summarised into a one-row table. Under `CELL_SET=epi` `cell_type`
    is not constant and comes back at the end of the list - one flag, one code path.
    """
    out = []
    for key in (LEIDEN_KEY, GROUP_KEY, LABEL_KEY):
        if key in obs and obs[key].astype(str).nunique() > 1:
            out.append(key)
    return out

# --------------------------------------------------------------------------- #
# The coordinate systems Route A can be read against
# --------------------------------------------------------------------------- #
#
# Route A asks where prior knowledge puts its cells along the axes of a latent space. The
# space is an ARGUMENT of that question, not part of it: A1 - A5 of 05_6 (scoring, within-
# stratum standardisation, target region, consensus, confounders) never touch an embedding,
# and only A6 does.
#
# ONLY DRVI IS REGISTERED IN THIS PHASE. 04 registers Harmony and PCA as well, written by its
# 04_9 embedding-control step; that step is not planned here, so registering them would point
# the flag at files nobody writes. The mechanism is kept rather than stripped out, because
# putting it back would be the harder change: adding a control run means one entry below plus
# a writer that produces the file in the shape 04's `run_harmony_epi.py` documents:
#
#     .X       cells x dimensions, cells in the order of shiao_tum.h5ad
#     .obs     the metadata of the subset
#     .var     `title` (the dimension label), `order` (0-based), `vanished` (bool)
#
# which is DRVI's own embedding format, so `analysis_dimensions()` and A6 read either
# without a branch. `vanished` is a DRVI concept; for every other method it is written all
# False and reported as such rather than faked.
#
# ROUTE B IS DRVI-ONLY AND STAYS THAT WAY. It reads the additive decoder off `embed.varm`,
# which a linear embedding does not have, so 05_7 and 05_8 are not parameterised and must not
# be run against one - `interpretability_scores()` raises on it by construction (the varm
# keys simply are not there).


@dataclass(frozen=True)
class Embedding:
    """One coordinate system Route A can be read against."""

    name: str            # the --embedding value
    run_id: str          # goes into every table and figure name, as `drvi_tum_32` does
    title: str           # for figure titles
    dim_prefix: str      # dimension labels: 'DR 1', 'HD 1', 'PC 1'
    n_dims: int
    description: str
    is_reference: bool = False   # True for the phase's own run: the one that owns the
                                 # embedding-INDEPENDENT tables of 05_6

    @property
    def embed_h5ad(self) -> Path:
        return TUM_DIR / f"embed_{self.run_id}.h5ad"


EMBEDDINGS = {
    "drvi": Embedding(
        name="drvi", run_id=DEFAULT_RUN_ID, title="DRVI", dim_prefix="DR",
        n_dims=N_LATENT, is_reference=True,
        description=f"05_3: DRVI, n_latent {N_LATENT}, trained on the 2,000 batch-aware "
                    "HVGs of 05_2"),
}

DEFAULT_EMBEDDING = "drvi"


def get_embedding(name: str) -> Embedding:
    if name not in EMBEDDINGS:
        raise KeyError(f"unknown embedding {name!r}; have: {', '.join(EMBEDDINGS)}")
    return EMBEDDINGS[name]


def set_embedding(emb) -> Embedding:
    """Point this module at one coordinate system, once, at the start of a step.

    Rebinds `RUN_ID` and `EMBED_H5AD`, which every writer below reads at call time, so the
    whole output naming of a step follows the flag. Deliberately a global rather than an
    argument threaded through fifteen call sites: a step interprets ONE embedding per run,
    and the run id is in every filename precisely so two runs cannot collide.

    With DRVI as the only registered space this is a no-op that keeps the call sites honest.
    """
    global RUN_ID, EMBED_H5AD
    emb = get_embedding(emb) if isinstance(emb, str) else emb
    RUN_ID = emb.run_id
    EMBED_H5AD = emb.embed_h5ad
    return emb


def add_embedding_argument(parser) -> None:
    """The `--embedding` flag, spelled the same way `--collection` is."""
    parser.add_argument(
        "--embedding", choices=sorted(EMBEDDINGS), default=DEFAULT_EMBEDDING,
        help=f"the coordinate system Route A is read against (default {DEFAULT_EMBEDDING})")

# Outputs. Everything a step writes is scoped to its collection, both in the folder and in
# the filename, so `scie` and `emt` can be re-run in any order without one touching the
# other's results - and so a table pulled out of the repo still says which readout it is from.
# The two per-cell files below are named with `DEFAULT_RUN_ID` and NOT with `RUN_ID`, and
# that is not an oversight: both are properties of the malignant OBJECT, not of any latent
# space. CytoTRACE2 is a measurement on raw counts and the signature scores come from
# `sc.tl.score_genes` on the all-genes matrix; neither has ever seen a dimension. Binding
# them to the embedding would make a control run re-score 36,192 cells for a byte-identical
# result, and leave two copies free to drift. The `drvi_tum_32` in their names is a label,
# not a dependency.
CYTOTRACE_CSV = TUM_DIR / f"cytotrace2_{DEFAULT_RUN_ID}.csv"   # per cell, 05_5: a readout, not a collection


def gmt_path(coll) -> Path:
    """The collection as actually used - mapped genes only - which is also the Appendix table."""
    return TABLE_DIR / coll.name / f"signatures_{coll.name}.gmt"


def scores_csv(coll) -> Path:
    """Per-cell raw + within-stratum z scores, 05_6. Heavy, so it lives outside the repo.

    Embedding-independent, hence `DEFAULT_RUN_ID`: see the note on CYTOTRACE_CSV above.
    """
    return TUM_DIR / f"signature_scores_{coll.name}_{DEFAULT_RUN_ID}.csv"


def enrichment_tsv(coll, n_top: int) -> Path:
    """The full ORA output of 05_7: every pair tested, before any significance filter."""
    return TUM_DIR / f"factor_first_top{n_top}_{coll.name}_{RUN_ID}.tsv"


def top_genes_tsv(n_top: int) -> Path:
    """The top-gene list per dimension-direction.

    Read off the DRVI decoder alone, so this one is deliberately NOT collection-scoped: both
    collections are tested against the same gene lists, and writing it twice would invite the
    two copies to drift apart.
    """
    return TUM_DIR / f"factor_first_top{n_top}_genes_{RUN_ID}.tsv"

# --------------------------------------------------------------------------- #
# The caveat that has to travel with every output of this step
# --------------------------------------------------------------------------- #

# 04's caveat says that no CNV inference has been run, which stays true of 04's tables and
# figures. This phase ran one, so that sentence is not the limit any more - and the honest
# replacement is not "no caveat", it is the list of what the CNV call itself does and does
# not establish. Three things, all of them consequences of how 05_1 and 05_2 were built:
#
#   * the call is a THRESHOLD on an inferCNV score, per cohort, not a validated clone
#     assignment: there is no matched normal, no bulk WGS, no per-clone reconstruction;
#   * cells no inferCNV run covered are `not_tested` and were EXCLUDED rather than assumed
#     normal, so the subset is what the call could reach, not every tumour cell present;
#   * 19 of 04's 29 cohorts clear the 200-cell floor, so the population is the patients with
#     the most callable tumour.
#
# What it buys, and this is the point of the phase: a state reported here is a state of
# aneuploid cells. It is no longer a mixture of normal and malignant epithelium described
# with a normal-breast vocabulary.
CAVEAT = (
    "Cells are the aneuploid epithelium called by inferCNV in 05_1: a per-cohort THRESHOLD "
    "on a CNV score, not a validated clone assignment - no matched normal, no bulk WGS. "
    "Cells no run covered (not_tested) were excluded rather than assumed normal, and 19 of "
    "the 29 cohorts of the epithelial compartment reach the 200-cell floor, so this is "
    "the tumour the call "
    "could reach in the patients with the most of it. The CellTypist labels carried along "
    "(cell_type_01_4) come from a NORMAL breast atlas with no malignant class and are "
    "landmarks, never identities: every state reported here is a state of malignant "
    "epithelium, named by its genes and not by that label."
)

CAVEAT_SHORT = (
    "Malignant epithelium as called by inferCNV (per-cohort threshold, 19 cohorts, not_tested\n"
    "excluded). CellTypist labels are landmarks from a normal-breast model, not identities."
)

# --------------------------------------------------------------------------- #
# Floors that apply to every collection
# --------------------------------------------------------------------------- #
#
# The signature registries themselves are in `sig_collections.py`. These two numbers are not
# collection-specific: they are what "this list is still the list it is named after" means on
# this object, and they are enforced identically for the SCIE lists and the EMT ones.

MIN_SIGNATURE_GENES = 10      # below this a signature is skipped and reported
MIN_MAPPED_FRACTION = 0.60    # below this the step stops: low coverage means NOT MEASURED

# --------------------------------------------------------------------------- #
# The target region, and the vote over its definitions
# --------------------------------------------------------------------------- #
#
# These three live here, and not in `cell_first_tum.py` where 04 wrote them, because any
# step that wants to ask a second question of the SAME consensus cell set has to be able to
# rebuild it from the same cached per-cell scores. In 04 that step is 04_9; here it is
# whatever comes after 05_8. Two copies of a quadrant rule are two chances to drift apart,
# and every such comparison rests on being asked about identical cells.


@dataclass(frozen=True)
class Cutoffs:
    """The quantiles that define 'high', 'low' and the middle band of a readout.

    Quantiles of the standardised scores rather than fixed z values: the scores are not
    normal and a fixed z would give wildly different group sizes across readouts, which
    would make the stability comparison meaningless.
    """

    high_q: float = 0.75
    low_q: float = 0.25
    mid_lo_q: float = 0.40
    mid_hi_q: float = 0.60

    @classmethod
    def from_args(cls, args) -> "Cutoffs":
        return cls(args.high_q, args.low_q, args.mid_lo_q, args.mid_hi_q)


def add_cutoff_arguments(parser) -> None:
    """The four quantile flags of 05_6.

    Kept on the shared module rather than in the step, so anything that later re-derives the
    target region spells them the same way: run 05_6 with a non-default region and a
    comparison told different numbers would be reading a different cell set than the one
    Route A reported.
    """
    d = Cutoffs()
    parser.add_argument("--high-q", type=float, default=d.high_q,
                        help=f"quantile of the within-stratum z-score above which a cell is "
                             f"'high' (default {d.high_q})")
    parser.add_argument("--low-q", type=float, default=d.low_q,
                        help=f"quantile below which a cell is 'low' (default {d.low_q})")
    parser.add_argument("--mid-lo-q", type=float, default=d.mid_lo_q,
                        help=f"lower edge of the 'mid' band, used by the EMT target region "
                             f"(default {d.mid_lo_q})")
    parser.add_argument("--mid-hi-q", type=float, default=d.mid_hi_q,
                        help=f"upper edge of the 'mid' band (default {d.mid_hi_q})")


def in_region(v: pd.Series, rule: str, cut: Cutoffs) -> pd.Series:
    """Whether each cell is at the high end, the low end, or in the middle band of a readout.

    "mid" selects the middle of an axis rather than an end of it. Nothing uses it now - the
    EMT target moved to co-expression, which is both more stable and less dependent on the
    hybrid lists - but the rule and its two flags stay, because that first definition is part
    of the record and re-running it has to remain a one-line change.
    """
    if rule == "high":
        return v >= v.quantile(cut.high_q)
    if rule == "low":
        return v <= v.quantile(cut.low_q)
    if rule == "mid":
        return (v >= v.quantile(cut.mid_lo_q)) & (v <= v.quantile(cut.mid_hi_q))
    raise ValueError(f"unknown region rule {rule!r}")


def define_target(z: pd.DataFrame, plane, cut: Cutoffs) -> pd.Series:
    """The target cell set of one plane: both readouts inside their respective regions."""
    return (in_region(z[f"z_{plane.x}"], plane.x_rule, cut)
            & in_region(z[f"z_{plane.y}"], plane.y_rule, cut))


def consensus_vote(qdf: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    """(votes per cell, called by a majority of the definitions).

    `>= ceil(n/2)`, floored at 2 so a single definition can never carry the consensus on its
    own. Changing this changes the size of the target set arithmetically before it changes
    anything biological - see the README on what CytoTRACE2 and the EMP removal did to it.
    """
    n_defs = qdf.sum(axis=1)
    return n_defs, n_defs >= max(2, int(np.ceil(qdf.shape[1] / 2)))


def effective_rank(X) -> tuple[float, float]:
    """(participation ratio of the coordinate correlation matrix, mean |r| off-diagonal).

    How many of a space's dimensions are actually independent of each other. It is here
    because the number of dimensions one readout's association is spread over cannot be read
    without it: correlated axes SHARE an association, so a space whose dimensions are
    redundant will spread any readout over more of them. A PCA is orthogonal by construction
    and scores exactly its own dimension count; a latent space need not, and DRVI does not.
    Nothing in 05_4 - 05_8 calls it yet; it is the half of 04_9's machinery that is worth
    keeping, because it is the only honest way to compare two spaces of different sizes.
    """
    r = np.corrcoef(np.asarray(X, dtype=np.float64), rowvar=False)
    r = np.nan_to_num(r, nan=0.0)            # a dead dimension has no correlation, not NaN
    np.fill_diagonal(r, 1.0)
    ev = np.linalg.eigvalsh(r)
    pr = float(ev.sum() ** 2 / np.square(ev).sum())
    off = r[~np.eye(r.shape[0], dtype=bool)]
    return pr, float(np.abs(off).mean())


# --------------------------------------------------------------------------- #
# Signature IO
# --------------------------------------------------------------------------- #


def read_signature_file(path: Path) -> list[str]:
    """One gene symbol per line, de-duplicated, order preserved.

    Several of the files are CRLF and at least one carries a UTF-8 BOM, so the
    encoding and the stripping are not optional: 'B2M\\r' matches no var_name.
    """
    seen, genes = set(), []
    with open(path, encoding="utf-8-sig") as fh:
        for line in fh:
            g = line.strip().strip('"').strip()
            if g and g not in seen:
                seen.add(g)
                genes.append(g)
    return genes


def load_signatures(coll, sig_dir: Path = SIG_DIR) -> dict[str, list[str]]:
    """One collection, as {name: [genes]}, straight from the text files, in registry order."""
    out = {}
    for name, fname in coll.files.items():
        path = sig_dir / fname
        if not path.exists():
            raise FileNotFoundError(f"signature file missing: {path}")
        genes = read_signature_file(path)
        if not genes:
            raise ValueError(f"signature file is empty: {path}")
        out[name] = genes
    return out


def write_gmt(sets: dict[str, list[str]], path: Path, descriptions: dict[str, str]) -> None:
    """GMT: name <tab> description <tab> gene <tab> gene ...

    The description field carries the provenance string, so the collection that feeds
    both routes is also the Appendix table.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        for name, genes in sets.items():
            fh.write("\t".join([name, descriptions.get(name, ""), *genes]) + "\n")


def read_gmt(path: Path) -> dict[str, list[str]]:
    """Back from GMT, dropping the description column."""
    sets = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 2:
                sets[parts[0]] = parts[2:]
    return sets


# --------------------------------------------------------------------------- #
# Table and figure writers, both stamping the caveat
# --------------------------------------------------------------------------- #


def table_dir(coll) -> Path:
    """tables/<collection>/. Created on demand, one folder per readout."""
    d = TABLE_DIR / coll.name
    d.mkdir(parents=True, exist_ok=True)
    return d


def table_path(name: str, coll, run_id: str | None = None) -> Path:
    """`run_id` defaults to the module's current one; pass it to reach another run's table.

    Nothing in this phase passes it today - it is what a step reading several runs in one
    process would need, since such a step cannot go through the module global.
    """
    return table_dir(coll) / f"{name}_{coll.name}_{run_id or RUN_ID}.csv"


def write_table(df: pd.DataFrame, name: str, coll, index: bool = True,
                run_id: str | None = None) -> Path:
    """Write a small result table to tables/<collection>/, caveat as a leading comment.

    Read it back with `pd.read_csv(path, comment='#', index_col=0)`. The comment lines
    are how the caveat travels with the table when it is pulled out of this folder.
    """
    run_id = run_id or RUN_ID
    path = table_path(name, coll, run_id)
    # The header names the space the table was computed in: `DRVI run drvi_tum_32`. A run
    # against any other registered space would say so instead, so a table pulled out of this
    # folder cannot be mistaken for the phase's own. A run id matching no registry entry
    # falls back to the bare id, which is the honest label for one.
    space = next((e.title for e in EMBEDDINGS.values() if e.run_id == run_id), None)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"# {name} | collection {coll.name} ({coll.title}) "
                 f"| {space + ' run' if space else 'run'} {run_id} "
                 f"| 05_4 - 05_8 signature interpretation\n")
        for line in CAVEAT.split(". "):
            if line.strip():
                fh.write(f"# CAVEAT: {line.strip().rstrip('.')}.\n")
        df.to_csv(fh, index=index)
    print(f"[table] {path}  ({df.shape[0]} x {df.shape[1]})")
    return path


def read_table(name: str, coll, run_id: str | None = None) -> pd.DataFrame:
    """Read back a table this stage wrote, dropping the caveat comment lines.

    05_7 and 05_8 read what earlier steps wrote; going through this rather than through a
    hand-built path is what keeps a step from silently reading the OTHER collection's table
    when it is run with the wrong flag - the file simply is not there.
    """
    return pd.read_csv(table_path(name, coll, run_id), comment="#", index_col=0)


def fig_dir(step: str, coll) -> Path:
    """figures/<step>/<collection>/, `step` being the full step-folder name, e.g. '05_6_cell_first'.

    The step folders mirror the Methods sections and stay one per step; the collection is a
    subfolder of each, so the SCIE and EMT versions of the same figure sit side by side
    without either being able to overwrite the other.
    """
    d = PHASE_DIR / "figures" / step / coll.name
    d.mkdir(parents=True, exist_ok=True)
    return d


def savefig(name: str, step: str, coll, fig=None, dpi: int = 300, caveat: bool = True,
            run_id: str | None = None):
    """Save a figure into figures/<step>/<collection>/, collection and run id appended.

    Same helper as 05_3 and 04_2 except for the footnote, which is the mandatory
    caveat: a figure showing cells or states must carry it wherever it ends up.

    `caveat=False` is for the figures that show neither. The two of 05_4 describe the
    signature collection itself - gene lists and their overlap, with no cell in them -
    so the CellTypist caveat has nothing to qualify there, and it was landing on top of
    the rotated x tick labels.
    """
    import matplotlib.pyplot as plt

    fig = plt.gcf() if fig is None else fig
    if caveat:
        fig.text(0.5, -0.02, CAVEAT_SHORT, ha="center", va="top", fontsize=6,
                 style="italic", color="0.35", linespacing=1.4)
    path = fig_dir(step, coll) / f"{name}_{coll.name}_{run_id or RUN_ID}.png"
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    print(f"[fig] {path}")
    return path


def dim_slug(dim: str) -> str:
    """'DR 7+' -> 'DR_07_pos'. Same convention as 03_3, so the two steps sort alike."""
    number, sign = dim[3:-1], dim[-1]
    return f"DR_{int(number):02d}_{'pos' if sign == '+' else 'neg'}"


def dim_sort_key(dim: str) -> tuple[int, int]:
    """Sort 'DR 7+' / 'DR 7-' numerically, positive direction first."""
    return int(dim[3:-1]), 0 if dim.endswith("+") else 1


# --------------------------------------------------------------------------- #
# DRVI interpretability scores
# --------------------------------------------------------------------------- #


def interpretability_scores(embed, gene_names, key: str = SCORE_KEY,
                            hide_vanished: bool = PRUNE_VANISHED):
    """Genes x dimension-directions scores from `embed.varm`, as DRVI's own accessor returns.

    Identical to the helper of 04_3 and 03_3: 05_3 stored the two score matrices per approach
    (one per direction), so the table `model.get_interpretability_scores(embed, adata)`
    returns is rebuilt from the embedding alone - no model, no GPU, no scvi-tools.

    One difference from 03_3, and it is deliberate: `hide_vanished` defaults to
    `PRUNE_VANISHED`, i.e. False, so every direction of every dimension comes back and
    the caller decides what to do with it. DRVI's own accessor drops a direction it
    marked vanished; here that flag is reported rather than acted on. Passing
    `hide_vanished=True` restores DRVI's behaviour: a direction is then dropped only if
    *that* direction vanished, since a dimension can carry a real program on one side and
    nothing on the other.
    """
    effect = np.concatenate([embed.varm[f"{key}_positive"], embed.varm[f"{key}_negative"]])

    info = (
        pd.concat([embed.var.assign(direction="+"), embed.var.assign(direction="-")])
        .assign(title=lambda df: df["title"] + df["direction"])
        .reset_index(drop=True)
    )
    info["keep"] = ~np.where(
        info["direction"] == "+",
        info["vanished_positive_direction"],
        info["vanished_negative_direction"],
    ) if hide_vanished else True

    return (
        pd.DataFrame(effect, columns=gene_names, index=info["title"])
        .loc[info.query("keep == True").sort_values(["order", "direction"])["title"]]
        .T
    )


def analysis_dimensions(embed) -> list[str]:
    """The dimensions this stage works on: all of them, in the embedding's own order.

    `PRUNE_VANISHED` is False, so nothing is dropped. When it is flipped back on, the
    vanished set is read programmatically from `var['vanished']` - never from a plot.
    """
    if PRUNE_VANISHED and "vanished" in embed.var:
        return embed.var.loc[~embed.var["vanished"].astype(bool), "title"].tolist()
    return embed.var["title"].tolist()


def n_vanished(embed) -> int:
    """How many dimensions DRVI flagged as vanished. Reported, not acted on.

    0 for a space that has no such notion - a PCA dimension is never 'vanished', it is only
    ever low-variance - rather than a KeyError. A hand-written embedding should carry the
    column as all-False for exactly this reason; the guard is for one that does not.
    """
    if "vanished" not in embed.var:
        return 0
    return int(embed.var["vanished"].astype(bool).sum())


def banner(title: str) -> None:
    print(f"\n{'=' * 78}\n{title}\n{'=' * 78}")
