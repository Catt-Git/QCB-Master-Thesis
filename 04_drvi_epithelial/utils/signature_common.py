"""Shared helpers for the signature-interpretation steps of the epithelial phase (04_3 - 04_7).

Lives in the phase's `utils/`, next to the scripts that use it, the same way
`02_integration_benchmark/utils/` holds `h5ad_compat.py` and `metrics_shared.py`. It is
imported by the five step scripts with the idiom those scripts use:

    UTILS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils")
    sys.path.insert(0, UTILS_DIR)
    import signature_common as C

Kept inside this phase rather than shared with 03, as everywhere else in Part 2: 04 reads
04_2's outputs but shares no code with 03_3, so the phase still reads as a self-contained
Materials & Methods section.

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
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------- #
# The run this step reads
# --------------------------------------------------------------------------- #

N_LATENT = 64                       # the 04_2 run of this phase, see drvi_epi.ipynb
DEFAULT_RUN_ID = f"drvi_epi_{N_LATENT}"

# The run every step writes under. This is the ONE mutable global of the module: 04_5 rebinds
# it through `set_embedding()` when it is asked to read a coordinate system other than DRVI's,
# so that every table and figure name follows the embedding without each writer having to be
# told about it. 04_6 and 04_7 never call `set_embedding`, so for them nothing moves: they
# read and write `drvi_epi_64` exactly as before.
RUN_ID = DEFAULT_RUN_ID

SEED = 0                            # scoring / subsampling seed, as in 01_4 and 03_1
SCORE_KEY = "OOD_combined"          # 04_2's interpretability scores, as in 03_3

# Vanished dimensions are NOT pruned in this stage. 03_2's `plot_pruned_umap_nonimm.py`
# measured what pruning does to a DRVI space and the answer was nothing - the vanished
# dimensions carry ~1e-05 of the latent variance - so pruning buys no cleanliness, while
# dropping dimensions before the correlations and the ORA silently decides, ahead of the
# analysis, which axes are allowed to mean something. 04_5 - 04_7 therefore run on all
# `N_LATENT` dimensions and all 2 x `N_LATENT` directions; `var['vanished']` and
# `var['vanished_*_direction']` are still read and reported, so any dimension that does
# come out significant can be checked against its flag. Set this to True to get the
# pruned behaviour back, unchanged, everywhere at once.
PRUNE_VANISHED = False

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #

# utils/signature_common.py -> 04_drvi_epithelial/ -> the repo root
PROJECT_DIR = Path(__file__).resolve().parents[2]
DATA_DIR = Path(os.environ.get("DATA_DIR", PROJECT_DIR / "datasets"))

EPI_DIR = DATA_DIR / "04_epi"                 # the heavy objects, outside the repo
SIG_DIR = Path(os.environ.get("SIGNATURE_DIR", DATA_DIR / "signatures"))

PHASE_DIR = PROJECT_DIR / "04_drvi_epithelial"

# Small result tables, versioned in the repo as appendix material. Phase-level, next to
# `figures/`, and for the same reason: 04_7 reads what 04_5 and 04_6 wrote, so a per-step
# `tables/` would mean steps reaching into each other's folders.
TABLE_DIR = PHASE_DIR / "tables"

# Inputs, all read-only here.
FULL_H5AD = EPI_DIR / "shiao_epi.h5ad"                  # 04_1: all genes, log-normalised
HVG_H5AD = EPI_DIR / "shiao_epi_hvg_2k.h5ad"            # 04_1: the DRVI training features
EMBED_H5AD = EPI_DIR / f"embed_{RUN_ID}.h5ad"           # 04_2: latent space + scores

# --------------------------------------------------------------------------- #
# The coordinate systems Route A can be read against
# --------------------------------------------------------------------------- #
#
# Route A asks where prior knowledge puts its cells along the axes of a latent space. The
# space is an ARGUMENT of that question, not part of it: A1 - A5 of 04_5 (scoring, within-
# stratum standardisation, target region, consensus, confounders) never touch an embedding,
# and only A6 does. Registering the alternatives here is what lets 04_9 ask the same
# question of Harmony without a second copy of the script.
#
# Adding a method means one entry below plus a writer that produces the file in the shape
# `04_9_embedding_control/run_harmony_epi.py` documents:
#
#     .X       cells x dimensions, cells in the order of shiao_epi.h5ad
#     .obs     the metadata of the compartment
#     .var     `title` (the dimension label), `order` (0-based), `vanished` (bool)
#
# which is DRVI's own embedding format, so `analysis_dimensions()` and A6 read either
# without a branch. `vanished` is a DRVI concept; for every other method it is written all
# False and reported as such rather than faked.
#
# ROUTE B IS DRVI-ONLY AND STAYS THAT WAY. It reads the additive decoder off `embed.varm`,
# which no method here has, so 04_6 and 04_7 are not parameterised and must not be run
# against these embeddings - `interpretability_scores()` raises on them by construction
# (the varm keys simply are not there).


@dataclass(frozen=True)
class Embedding:
    """One coordinate system Route A can be read against."""

    name: str            # the --embedding value
    run_id: str          # goes into every table and figure name, as `drvi_epi_64` does
    title: str           # for figure titles
    dim_prefix: str      # dimension labels: 'DR 1', 'HD 1', 'PC 1'
    n_dims: int
    description: str
    is_reference: bool = False   # True for the phase's own run: the one that owns the
                                 # embedding-INDEPENDENT tables of 04_5

    @property
    def embed_h5ad(self) -> Path:
        return EPI_DIR / f"embed_{self.run_id}.h5ad"


# 50 for both linear spaces: `scib.integration.harmony` calls `sc.tl.pca(adata)` and takes
# scanpy's default, so this is the size the phase-02 benchmark ran Harmony at. The PCA is
# nested, so it is not a handicap against DRVI's 64 on any per-readout maximum; it does have
# to be normalised away in anything that COUNTS dimensions, which is what the effective
# dimension count of `compare_embeddings_epi.py` is for.
LINEAR_N_DIMS = 50

EMBEDDINGS = {
    "drvi": Embedding(
        name="drvi", run_id=DEFAULT_RUN_ID, title="DRVI", dim_prefix="DR",
        n_dims=N_LATENT, is_reference=True,
        description="04_2: DRVI, n_latent 64, trained on the 2,000 batch-aware HVGs of 04_1"),
    "harmony": Embedding(
        name="harmony", run_id=f"harmony_epi_{LINEAR_N_DIMS}", title="Harmony",
        dim_prefix="HD", n_dims=LINEAR_N_DIMS,
        description="04_9: Harmony on the same 2,000 HVGs, batch-corrected PCs"),
    "pca": Embedding(
        name="pca", run_id=f"pca_epi_{LINEAR_N_DIMS}", title="PCA (unintegrated)",
        dim_prefix="PC", n_dims=LINEAR_N_DIMS,
        description="04_9: the PCA Harmony corrects, kept uncorrected as the null arm"),
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
# that is not an oversight: both are properties of the epithelial OBJECT, not of any latent
# space. CytoTRACE2 is a measurement on raw counts and the signature scores come from
# `sc.tl.score_genes` on the all-genes matrix; neither has ever seen a dimension. Binding
# them to the embedding would make a Harmony run re-score 74,441 cells for a byte-identical
# result, and leave two copies free to drift. The `drvi_epi_64` in their names is historical.
CYTOTRACE_CSV = EPI_DIR / f"cytotrace2_{DEFAULT_RUN_ID}.csv"   # per cell, 04_4: a readout, not a collection


def gmt_path(coll) -> Path:
    """The collection as actually used - mapped genes only - which is also the Appendix table."""
    return TABLE_DIR / coll.name / f"signatures_{coll.name}.gmt"


def scores_csv(coll) -> Path:
    """Per-cell raw + within-stratum z scores, 04_5. Heavy, so it lives outside the repo.

    Embedding-independent, hence `DEFAULT_RUN_ID`: see the note on CYTOTRACE_CSV above.
    """
    return EPI_DIR / f"signature_scores_{coll.name}_{DEFAULT_RUN_ID}.csv"


def enrichment_tsv(coll, n_top: int) -> Path:
    """The full ORA output of 04_6: every pair tested, before any significance filter."""
    return EPI_DIR / f"factor_first_top{n_top}_{coll.name}_{RUN_ID}.tsv"


def top_genes_tsv(n_top: int) -> Path:
    """The top-gene list per dimension-direction.

    Read off the DRVI decoder alone, so this one is deliberately NOT collection-scoped: both
    collections are tested against the same gene lists, and writing it twice would invite the
    two copies to drift apart.
    """
    return EPI_DIR / f"factor_first_top{n_top}_genes_{RUN_ID}.tsv"

# --------------------------------------------------------------------------- #
# The caveat that has to travel with every output of this step
# --------------------------------------------------------------------------- #

CAVEAT = (
    "Cell type labels come from CellTypist (Cells_Adult_Breast.pkl, Kumar et al. 2023), "
    "a NORMAL adult breast atlas with no malignant class: TNBC cells are assigned to the "
    "nearest normal state. This subset therefore mixes normal and malignant epithelium and "
    "no CNV inference has been run. Every state reported here is an EPITHELIAL state, "
    "not a tumour cell state."
)

CAVEAT_SHORT = (
    "Epithelial state, not a tumour cell state: CellTypist labels come from a normal breast\n"
    "atlas with no malignant class, and no CNV inference has been run."
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
# These three live here, and not in `cell_first_epi.py` where they were written, because
# 04_9 has to reconstruct the SAME consensus cell set from the same cached per-cell scores
# in order to ask a second question of it (can the whole space predict it, not just one
# dimension?). Two copies of a quadrant rule are two chances to drift apart, and the whole
# comparison rests on the two spaces being asked about identical cells.


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
    """The four quantile flags, spelled identically by 04_5 and 04_9.

    04_9 needs them because it re-derives the target region: run with a non-default region
    in 04_5 and the comparison has to be told the same numbers, or it would be comparing
    the spaces on a different cell set than the one Route A reported.
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
    because `effective_n_dims` in 04_9 - the number of dimensions one readout's association
    is spread over - cannot be read without it: correlated axes SHARE an association, so a
    space whose dimensions are redundant will spread any readout over more of them. A PCA is
    orthogonal by construction and scores exactly its own dimension count; a latent space
    need not, and DRVI does not.
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

    Only `04_9_embedding_control/compare_embeddings_epi.py` passes it, and it has to: it is
    the one step that reads several runs in the same process, so it cannot go through the
    global.
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
    # The header names the space the table was computed in. `DRVI run drvi_epi_64` is what
    # it has always said and still says; a control run of 04_5 says `Harmony run
    # harmony_epi_50` instead, so a table pulled out of this folder cannot be mistaken for
    # the phase's own. A cross-run table (04_9's comparison) matches no registry entry and
    # falls back to the bare run id, which is the honest label for one.
    space = next((e.title for e in EMBEDDINGS.values() if e.run_id == run_id), None)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"# {name} | collection {coll.name} ({coll.title}) "
                 f"| {space + ' run' if space else 'run'} {run_id} "
                 f"| 04_3 - 04_7 signature interpretation\n")
        for line in CAVEAT.split(". "):
            if line.strip():
                fh.write(f"# CAVEAT: {line.strip().rstrip('.')}.\n")
        df.to_csv(fh, index=index)
    print(f"[table] {path}  ({df.shape[0]} x {df.shape[1]})")
    return path


def read_table(name: str, coll, run_id: str | None = None) -> pd.DataFrame:
    """Read back a table this stage wrote, dropping the caveat comment lines.

    04_6 and 04_7 read what earlier steps wrote; going through this rather than through a
    hand-built path is what keeps a step from silently reading the OTHER collection's table
    when it is run with the wrong flag - the file simply is not there.
    """
    return pd.read_csv(table_path(name, coll, run_id), comment="#", index_col=0)


def fig_dir(step: str, coll) -> Path:
    """figures/<step>/<collection>/, `step` being the full step-folder name, e.g. '04_5_cell_first'.

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

    Same helper as 03_3 and 04_2 except for the footnote, which is the mandatory
    caveat: a figure showing cells or states must carry it wherever it ends up.

    `caveat=False` is for the figures that show neither. The two of 04_3 describe the
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

    Identical to the helper of 03_3: 04_2 stored the two score matrices per approach
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
    ever low-variance - rather than a KeyError. The written embeddings of 04_9 carry the
    column as all-False for exactly this reason; the guard is for anything hand-built.
    """
    if "vanished" not in embed.var:
        return 0
    return int(embed.var["vanished"].astype(bool).sum())


def banner(title: str) -> None:
    print(f"\n{'=' * 78}\n{title}\n{'=' * 78}")
