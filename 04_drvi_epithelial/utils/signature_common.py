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

Nothing here computes anything on its own; it is the paths, the signature registry, the
figure/table writers and the interpretability-score accessor the five steps all need.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------- #
# The run this step reads
# --------------------------------------------------------------------------- #

N_LATENT = 64                       # the 04_2 run of this phase, see drvi_epi.ipynb
RUN_ID = f"drvi_epi_{N_LATENT}"

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

# Outputs written next to the objects they come from (per-cell, heavy).
GMT_PATH = TABLE_DIR / "lab_signatures.gmt"
SCORES_CSV = EPI_DIR / f"signature_scores_{RUN_ID}.csv"        # per cell, raw + z
CYTOTRACE_CSV = EPI_DIR / f"cytotrace2_{RUN_ID}.csv"           # per cell, 04_3 step 02

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
# The lab's signature collection
# --------------------------------------------------------------------------- #
#
# One plain-text file per signature, one gene symbol per line. `file` is the name on
# disk (two of them do not match the signature name), `axis` is which of the two
# readouts it belongs to, `provenance` goes verbatim into the .gmt description field
# and doubles as the Appendix table.
#
# NOTE. There is deliberately NO stemness consensus signature. The intersection of
# BENPORATH_ES1 / ESC_WONG / ESC_ASSOU was tested by the lab and captures
# proliferation only, which the existing S_score / G2M_score already cover. Each
# stemness signature is used on its own and none of them is primary at this stage.

SIGNATURES = [
    # name,                  file,                    axis,       provenance
    ("HALLMARK_IFNA",         "HALLMARK_IFNA.txt",      "immune",   "Interferon Alpha response, MSigDB Hallmark"),
    ("HALLMARK_IFNG",         "HALLMARK_IFNG.txt",      "immune",   "Interferon Gamma response, MSigDB Hallmark"),
    ("ISDS",                  "ISDS.txt",               "immune",   "IFN-Stem Cell-Down signature, PMC5481166"),
    ("KEGG_APM",              "KEGG_APM.txt",           "immune",   "Antigen Presentation Machinery, KEGG"),
    ("IMMUNOGENIC_CONSENSUS", "Immune consensus.txt",   "immune",   "Curated by the lab from the KEGG APM signature, retaining only immunogenic genes and immunomodulators; primary immune readout"),
    ("BENPORATH_ES1",         "BENPORATH_ES1.txt",      "stemness", "MSigDB BENPORATH_ES_1"),
    ("ESC_ASSOU",             "ESC_ASSOU.txt",          "stemness", "PMC1906587, Table S3"),
    ("ESC_WONG",              "ESC_WONG.txt",           "stemness", "MSigDB WONG_EMBRYONIC_STEM_CELL_CORE"),
    ("EMP",                   "EMP.txt",                "stemness", "Embryonic Multipotent Progenitors, PMID 29784918"),
    ("LIM_STEM",              "LIM_STEM.txt",           "stemness", "MSigDB LIM_MAMMARY_STEM_CELL_UP"),
    ("FMASC",                 "fMaSC.txt",              "stemness", "Fetal mammary stem cells, PMC3277444, Supplementary Table 2"),
]

# The immune axis is expected to be LOW in the states of interest: the project looks for
# immune *evasion*, not immunogenicity.
PRIMARY_IMMUNE = "IMMUNOGENIC_CONSENSUS"

IMMUNE_SIGS = [n for n, _, ax, _ in SIGNATURES if ax == "immune"]
STEMNESS_SIGS = [n for n, _, ax, _ in SIGNATURES if ax == "stemness"]
SIG_AXIS = {n: ax for n, _, ax, _ in SIGNATURES}
SIG_PROVENANCE = {n: p for n, _, _, p in SIGNATURES}

MIN_SIGNATURE_GENES = 10      # below this a signature is skipped and reported
MIN_MAPPED_FRACTION = 0.60    # below this the step stops: low coverage means NOT MEASURED

# --------------------------------------------------------------------------- #
# Signature IO
# --------------------------------------------------------------------------- #


def read_signature_file(path: Path) -> list[str]:
    """One gene symbol per line, de-duplicated, order preserved.

    Five of the eleven files are CRLF and at least one carries a UTF-8 BOM, so the
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


def load_signatures(sig_dir: Path = SIG_DIR) -> dict[str, list[str]]:
    """The whole collection, as {name: [genes]}, straight from the text files."""
    out = {}
    for name, fname, _, _ in SIGNATURES:
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


def write_table(df: pd.DataFrame, name: str, index: bool = True, table_dir: Path = TABLE_DIR) -> Path:
    """Write a small result table to the phase's tables/, caveat as a leading comment.

    Read it back with `pd.read_csv(path, comment='#', index_col=0)`. The comment lines
    are how the caveat travels with the table when it is pulled out of this folder.
    """
    table_dir.mkdir(parents=True, exist_ok=True)
    path = table_dir / f"{name}_{RUN_ID}.csv"
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"# {name} | DRVI run {RUN_ID} | 04_3 signature interpretation\n")
        for line in CAVEAT.split(". "):
            if line.strip():
                fh.write(f"# CAVEAT: {line.strip().rstrip('.')}.\n")
        df.to_csv(fh, index=index)
    print(f"[table] {path}  ({df.shape[0]} x {df.shape[1]})")
    return path


def fig_dir(step: str) -> Path:
    """figures/<step>/, `step` being the full step-folder name, e.g. '04_5_cell_first'."""
    d = PHASE_DIR / "figures" / step
    d.mkdir(parents=True, exist_ok=True)
    return d


def savefig(name: str, step: str, fig=None, dpi: int = 300, caveat: bool = True):
    """Save a figure into figures/<step>/, run id appended, caveat as a footnote.

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
    path = fig_dir(step) / f"{name}_{RUN_ID}.png"
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
    if PRUNE_VANISHED:
        return embed.var.loc[~embed.var["vanished"].astype(bool), "title"].tolist()
    return embed.var["title"].tolist()


def n_vanished(embed) -> int:
    """How many dimensions DRVI flagged as vanished. Reported, not acted on."""
    return int(embed.var["vanished"].astype(bool).sum())


def banner(title: str) -> None:
    print(f"\n{'=' * 78}\n{title}\n{'=' * 78}")
