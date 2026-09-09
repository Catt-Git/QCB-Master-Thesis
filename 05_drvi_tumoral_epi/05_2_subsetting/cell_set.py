"""05_2: which cells this phase is about, and what its files are called.

Phase 05 can be run over two different cell sets, and the difference is one flag rather
than two copies of the chain:

  CELL_SET=tum     (default)  the MALIGNANT cells only - the primary line of the phase
  CELL_SET=epi                every epithelial cell, malignant ones included, under the
                              post-CNV labels - i.e. phase 04 redone with the labels fixed

`tum` is the default because of what the phase is for. 05_6/05_7 (cell-first and
factor-first) read the DRVI latent dimensions and ask which gene programmes load on them.
If DRVI is trained on all epithelium, its dominant dimensions encode malignant-versus-normal
- a contrast that is nearly CONSTANT inside the malignant subset - so those two steps would
be interpreting the leftovers. It is the same argument 04 makes for existing at all against
03: a latent space trained on a mixture describes the mixture.

`epi` is a control, not a second main line. What it answers is "how much did the wrong
labels cost phase 04?", which is a methods question worth a paragraph and not the biology.
Its results are named differently, live in different files, and never overwrite `tum`.

A second flag, `HVG_SET`, picks WHICH 2,000 genes the DRVI input is made of - unset for the
panel `reduce_data_tum.py` selected, `nomt` for the same panel with the mitochondrial genes
replaced (see `_HVG_SETS` and `hvg_no_mt.py`). It is orthogonal to CELL_SET: same cells, same
normalization, same clustering, a different gene panel and therefore a different run id.

A third, `N_LATENT`, picks the DRVI latent size. It is not a property of the cells at all -
it is 05_3's hyperparameter - and it lives here for one reason: it is the third segment of
the run id, and `run_id()` below is the only place that string is built. 05_3 takes it as
`--n-latent` and the flag wins; everything downstream of 05_3 has no flag and reads the
environment, so ONE variable moves the whole chain from `drvi_tum_32` to `drvi_tum_64`.

The three are independent, and the run id carries all three, so no two combinations can
overwrite each other:

    CELL_SET  HVG_SET  N_LATENT   run id
    tum       -        32         drvi_tum_32          the run this phase reports
    tum       -        64         drvi_tum_64
    epi       nomt     64         drvi_epicnv_64_nomt

This module exists so that the mapping CELL_SET -> file prefix lives in exactly one place.
Four scripts write into the same directory; a prefix computed independently in each of them
is a silent-overwrite bug waiting to happen.

Every path is under $DATA_DIR/05_tum/, next to what 05_1 wrote.
"""

from __future__ import annotations

import os
from pathlib import Path

# The eleven epithelial labels of the CellTypist model, verbatim as in 04_1 and in
# 01_4/fraction_reassignment.py. Used only by CELL_SET=epi; `tum` selects on cnv_status.
EPITHELIAL_CELL_TYPES = {
    "LummHR-SCGB", "LummHR-active", "LummHR-major",
    "Lumsec-HLA", "Lumsec-KIT", "Lumsec-basal", "Lumsec-lac", "Lumsec-major",
    "Lumsec-myo", "Lumsec-prol",
    "basal",
}

MALIGNANT_LABEL = "malignant"

# The immune lineage of the same CellTypist model, verbatim from
# 01_4/fraction_reassignment.py. Needed because `fraction` in shiao.h5ad is a FUNCTION of
# `cell_type`, and this phase replaces `cell_type` with the post-CNV label: the inherited
# column is stale the moment the join happens. Recomputing it is the whole of
# fraction_reassignment.py - it has no state and no model, it is a lookup - so it is inlined
# here rather than re-run as a script.
#
# It matters only for CELL_SET=epi. Under `tum` every cell was epithelial and therefore
# 'non_imm' already; under `epi` the subset can pull in cells the 01_4 annotation had called
# immune (2,821 of them changed lineage in the re-annotation), and those would arrive carrying
# 'imm'.
IMMUNE_CELL_TYPES = {
    "CD4-Tem", "CD4-Th", "CD4-Th-like", "CD4-Treg", "CD4-activated", "CD4-naive",
    "CD8-Tem", "CD8-Trm", "CD8-activated", "GD", "NK", "NK-ILCs", "NKT", "T_prol",
    "Macro-IFN", "Macro-lipo", "Macro-m1", "Macro-m1-CCL", "Macro-m2", "Macro-m2-CXCL",
    "Mast", "Mono-classical", "Mono-non-classical", "Neutrophil",
    "cDC1", "cDC2", "mDC", "pDC", "mye-prol",
    "b_naive", "bmem_switched", "bmem_unswitched", "plasma_IgA", "plasma_IgG",
}


def reassign_fraction(labels):
    """`fraction` recomputed from a post-CNV `cell_type`. 01_4's rule, plus MALIGNANT_LABEL.

    Malignant cells are non-immune by construction: only epithelial cells were eligible for
    the call in 05_1. Anything not in the immune set and not malignant is non-immune too,
    which is 01_4's own convention.
    """
    import pandas as pd
    return pd.Categorical(
        ["imm" if v in IMMUNE_CELL_TYPES else "non_imm" for v in labels],
        categories=["imm", "non_imm"],
    )

# obs columns this phase writes onto the subset. `cell_type` is the POST-CNV label, so
# every script inherited from 04 keeps reading the column it expects; the pre-CNV label is
# kept beside it under its own name because it is the only non-constant grouping the
# malignant subset has (see NMI_LABEL_KEY below).
LABEL_KEY = "cell_type"              # = cell_type_cnv; constant 'malignant' when CELL_SET=tum
PRIOR_LABEL_KEY = "cell_type_01_4"   # the contaminated CellTypist label, kept for grouping
NMI_LABEL_KEY = PRIOR_LABEL_KEY      # what 05_2 clustering maximises NMI against
STATUS_KEY = "cnv_status"
COMPARTMENT_KEY = "compartment"
BATCH_KEY = "cohort"
TREATMENT_KEY = "treatment"

# Cohort filter. 200 is inherited from 04_1, and on this subset the choice is almost free:
# the per-cohort malignant counts are bimodal (a cohort has thousands of malignant cells or
# a handful, with next to nothing between), so
#     >=  50 -> 27 cohorts, 36,901 cells
#     >= 100 -> 21 cohorts, 36,433 cells
#     >= 200 -> 19 cohorts, 36,192 cells
#     >= 300 -> 18 cohorts, 35,991 cells
# Dropping from 100 to 200 costs 241 cells out of 36,433 (0.7%) and buys back sane per-cohort
# batch parameters for DRVI, which treats `cohort` as its batch key and has to estimate them
# somewhere. 200 it is.
MIN_CELLS_PER_COHORT = 200

# 04_1 also dropped cohorts missing any of BASE/PD1/RTPD1, because a treatment phase was
# planned downstream of it. That phase is not planned any more, and on the malignant subset
# the filter is actively harmful: several cohorts clear it with two cells in a timepoint
# while Patient63 would lose 1,673 malignant cells for having none. Off here, deliberately.
DROP_INCOMPLETE_COHORTS = False
REQUIRED_TREATMENTS = ("BASE", "PD1", "RTPD1")

MIN_GENES = 100     # cell filter, as 04_1
MAX_PCT_MT = 10     # cell filter, as 04_1
MIN_CELLS = 3       # gene filter, as 04_1
N_HVGS = 2000
SEED = 0

_SETS = {
    "tum": ("shiao_tum", "tum", "the malignant cells"),
    "epi": ("shiao_epicnv", "epicnv", "every epithelial cell under the post-CNV labels"),
}

# A second, orthogonal knob: WHICH 2,000 genes the DRVI input is made of. It is not a cell
# set - the cells, the normalization and the clustering are identical - so it does not go
# through _SETS; it only appends a tag to the two files reduce_data_tum.py writes and to the
# run id 05_3 derives from them.
#
#   HVG_SET unset / ''   the panel reduce_data_tum.py selected, mitochondrial genes included
#   HVG_SET=nomt         the same panel with the 11 MT- genes dropped and the next 11 in the
#                        same batch-aware ranking put in their place, so it is still 2,000
#                        (built by hvg_no_mt.py, which is the only script that writes it)
#
# The MT- genes are transcription of the mitochondrial genome: they track how much of a
# cell's RNA came from mitochondria, which is a dissociation-stress and lysis readout rather
# than a cell state, and DRVI has no reason to spend a latent dimension on it. `pct_counts_mt`
# already filters on them upstream (MAX_PCT_MT), and the metallothioneins (MT1*, MT2A) and
# MTRNR2L12 are nuclear genes that only share the prefix - they are NOT dropped.
_HVG_SETS = {
    "": "",
    "nomt": "_nomt",
}


def cell_set() -> str:
    """The active cell set, from $CELL_SET. Defaults to 'tum'."""
    value = os.environ.get("CELL_SET", "tum").strip().lower() or "tum"
    if value not in _SETS:
        raise SystemExit(
            f"CELL_SET={value!r} is not one of {sorted(_SETS)}; "
            "unset it for the malignant subset (the default) or set CELL_SET=epi"
        )
    return value


def prefix(value: str | None = None) -> str:
    """File prefix for the active cell set: 'shiao_tum' or 'shiao_epicnv'."""
    return _SETS[value or cell_set()][0]


def compartment(value: str | None = None) -> str:
    """What goes into obs['compartment']: 'tum' or 'epicnv'."""
    return _SETS[value or cell_set()][1]


def describe(value: str | None = None) -> str:
    return _SETS[value or cell_set()][2]


def data_dir() -> Path:
    """$DATA_DIR, which the scripts never default: they abort if it is unset."""
    return Path(os.environ["DATA_DIR"]).expanduser().resolve()


def tum_dir() -> Path:
    """$DATA_DIR/05_tum, where 05_1 already wrote and where this phase keeps everything."""
    d = data_dir() / "05_tum"
    d.mkdir(parents=True, exist_ok=True)
    return d


def path(suffix: str, value: str | None = None) -> Path:
    """`path('_norm.h5ad')` -> $DATA_DIR/05_tum/shiao_tum_norm.h5ad (or shiao_epicnv_...)."""
    return tum_dir() / f"{prefix(value)}{suffix}"


def hvg_set() -> str:
    """The active HVG panel variant, from $HVG_SET. Defaults to '' (the original panel)."""
    value = os.environ.get("HVG_SET", "").strip().lower()
    if value not in _HVG_SETS:
        raise SystemExit(
            f"HVG_SET={value!r} is not one of {sorted(k for k in _HVG_SETS if k)}; "
            "unset it for the panel reduce_data_tum.py selected"
        )
    return value


def hvg_tag() -> str:
    """The filename/run-id tag of the active panel: '' or '_nomt'."""
    return _HVG_SETS[hvg_set()]


def hvg_path(suffix: str, value: str | None = None) -> Path:
    """`hvg_path('.h5ad')` -> .../shiao_epicnv_hvg_2k.h5ad, or ..._hvg_2k_nomt.h5ad."""
    return path(f"_hvg_2k{hvg_tag()}{suffix}", value)


# The size 05_3 trained the run this phase reports at. 05_3 judges it by how many dimensions
# DRVI lets vanish, which is why it is a knob and not a constant: see that step's README.
DEFAULT_N_LATENT = 32


def n_latent() -> int:
    """The DRVI latent size, from $N_LATENT. Defaults to DEFAULT_N_LATENT.

    Read by 05_4 - 05_8, which take no flag of their own: they do not train anything, they
    read what 05_3 wrote, and what they need is its NAME. 05_3 itself defaults `--n-latent`
    to this and overrides it when the flag is given.
    """
    raw = os.environ.get("N_LATENT", "").strip()
    if not raw:
        return DEFAULT_N_LATENT
    try:
        value = int(raw)
    except ValueError:
        raise SystemExit(f"N_LATENT={raw!r} is not an integer") from None
    if value <= 0:
        raise SystemExit(f"N_LATENT={value} must be positive")
    return value


def run_id(n: int | None = None, value: str | None = None) -> str:
    """`drvi_<compartment>_<n_latent><hvg_tag>` - the id 05_3 names every output after.

    THE ONLY PLACE THIS STRING IS BUILT, for the reason the whole module exists: 05_3 writes
    the model, the embedding and the downstream object under it, and 05_4 - 05_8 find them by
    reconstructing it. Two spellings of it is a step reading last week's run and saying
    nothing. `n` is for 05_3, which has the value on its command line and must not go through
    the environment to get it back.
    """
    return f"drvi_{compartment(value)}_{n or n_latent()}{hvg_tag()}"


def banner(step: str) -> None:
    """Every script prints the same three lines, so a log says which set it was run on."""
    s = cell_set()
    print(f"{step}  |  CELL_SET={s} ({describe(s)})", flush=True)
    print(f"DATA_DIR : {data_dir()}", flush=True)
    print(f"prefix   : {prefix(s)}", flush=True)
    if hvg_set():
        print(f"HVG_SET  : {hvg_set()} (panel tag {hvg_tag()!r})", flush=True)
    if os.environ.get("N_LATENT", "").strip():
        print(f"N_LATENT : {n_latent()}", flush=True)
    print(flush=True)
