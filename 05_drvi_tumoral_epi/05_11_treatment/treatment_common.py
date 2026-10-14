"""05_11: what the two treatment scripts share - the sample, the readouts, the hypotheses.

Two scripts ask the treatment question of the same latent space:

    treatment_programmes_tum.py   does a NAMED programme move with the treatment?
                                  (per-patient pseudobulk of the DRVI directions, paired test)
    milo_treatment_tum.py         does the COMPOSITION of the malignant states move?
                                  (Milo: differential abundance on the DRVI k-NN graph)

Everything they must agree on lives here, for the reason `cell_set.py` gives for itself: a
sample definition, a cell floor or a programme list written twice is two chances to drift,
and the two analyses are only readable side by side if they were asked about the same cells.

WHY THE JOINT RUN AND NOT ONE DRVI PER TIMEPOINT. `drvi_tum_64_nomt` was trained with
`batch_key='cohort'`, i.e. it corrects the PATIENT. The timepoint varies inside the patient
and was never given to the model, so whatever the treatment did to the tumour is still in the
latent space - and BASE, PD1 and RTPD1 cells sit in ONE coordinate system. Three separate runs
would give three sets of axes with no correspondence between them, trained on very different
numbers of cells, so a programme missing at RTPD1 could not be told apart from a programme
the smaller run had too few cells to learn.

THE UNIT IS THE BIOPSY, NOT THE CELL. A cell-level test on 42,000 cells calls everything
significant and is pseudoreplication: the cells of one biopsy are not independent draws. The
sample here is `cohort x treatment` - one biopsy per patient per timepoint (its CD45+ and
CD45- libraries pooled, which is right: they are the same piece of tumour) - and every
statistic is computed over those, paired within the patient.

THE PROGRAMMES ARE NOT RE-CHOSEN HERE. They are read from 05_8's
`programme_dimensions_consensus_<run>.csv` - the table `programmes_named` is drawn from - so
the directions tested are exactly the bars of that figure. The default table is the one the
thesis reports: vanished dimensions pruned, Route A bar at rho 0.30 (`OUT_TAG=_pruned_rho030`).
Set PRUNE_VANISHED / OUT_TAG in the environment to read another one; the outputs follow the
same tag, so they land beside the barplot they extend.
"""

from __future__ import annotations

import os
import re
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

# The programme table this step reads is the one `programmes_named` was drawn from. Both
# variables are resolved by signature_common at IMPORT time, so they have to be in the
# environment before the import below - the same constraint dr_umap_tum.py works around.
# `setdefault`, not assignment: an explicit export still wins.
os.environ.setdefault("PRUNE_VANISHED", "1")
os.environ.setdefault("OUT_TAG", "_pruned_rho030")

STEP_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(STEP_DIR.parent / "utils"))
import signature_common as C          # noqa: E402
import sig_collections as SC          # noqa: E402

STEP = "05_11_treatment"

# Outputs belong to the run and to no single collection, like 05_8's barplot: they land in
# figures<OUT_TAG>/05_11_treatment/consensus/<run>/ and tables<OUT_TAG>/consensus/<run>/.
CONSENSUS = SimpleNamespace(name="consensus",
                            title="all three collections, read together")

TREATMENT_KEY = C.CS.TREATMENT_KEY    # 'treatment'
RESPONSE_KEY = "response"
BATCH_KEY = C.BATCH_KEY               # 'cohort' = the patient

# Shiao et al. 2024: a biopsy before treatment, one after a single cycle of pembrolizumab,
# one after pembrolizumab plus radiotherapy. The order is the order of the trial.
TREATMENTS = ("BASE", "PD1", "RTPD1")
RESPONSES = ("NR", "R1", "R2")

# (later, earlier). RTPD1 vs PD1 is the one that isolates the radiotherapy: both biopsies are
# under anti-PD1, only one has had RT.
CONTRASTS = (("PD1", "BASE"), ("RTPD1", "BASE"), ("RTPD1", "PD1"))


def contrast_name(later: str, earlier: str) -> str:
    return f"{later}_vs_{earlier}"


CONTRAST_NAMES = tuple(contrast_name(*c) for c in CONTRASTS)

# A biopsy below this many malignant cells is not summarised: its mean is a handful of cells.
# 20 keeps most patient x timepoint pairs (the malignant counts are bimodal per biopsy, as
# they are per cohort) and is printed with what it drops. `--min-cells` moves it.
MIN_CELLS_PER_SAMPLE = 20

# A paired test with fewer pairs than this is not run. With 5 pairs the smallest two-sided
# Wilcoxon p is 0.0625 - the test cannot reject - so below it the row is reported as
# not tested rather than as not significant.
MIN_PAIRS = 6

# Categorical slots 1-3 of the reference palette, which validate all-pairs in light mode;
# fixed per entity so a response group keeps its colour in every figure of the step.
RESPONSE_COLORS = {"NR": "#2a78d6", "R1": "#eb6834", "R2": "#1baf7a"}

# Diverging blue <-> red with a grey midpoint, for paired deltas and associations.
DIVERGING = ("#184f95", "#6da7ec", "#f0efec", "#e98584", "#b02a29")

# --------------------------------------------------------------------------- #
# The hypotheses fixed BEFORE looking at the data
# --------------------------------------------------------------------------- #
#
# One (programme, contrast, expected sign) per row. They are their own Benjamini-Hochberg
# family: a handful of tests stated in advance, not 3 x 38 directions searched afterwards.
# Everything else this step tests is the exploratory family and is reported as such.
#
#   IFN-gamma / MHC-II after anti-PD1    the tumour answering the IFN-gamma of re-activated T
#                                        cells; the trial's early responders (R1) are where
#                                        it should be strongest
#   type-I IFN after RT                  radiotherapy induces type-I interferon in the tumour
#                                        cell (cytosolic DNA sensing); the patients who only
#                                        responded once RT was added (R2) are where to look
#   senescence after RT                  RT-induced senescence
#   cell cycle down after RT             RT-induced arrest
HYPOTHESES = pd.DataFrame([
    ("HALLMARK_IFNG", "PD1_vs_BASE", +1, "IFN-gamma response after anti-PD1"),
    ("MP17_INTERFERON_MHC_II_1", "PD1_vs_BASE", +1, "IFN / MHC-II after anti-PD1"),
    ("HALLMARK_IFNA", "RTPD1_vs_PD1", +1, "type-I IFN after RT"),
    ("MP19_EPITHELIAL_SENESCENCE", "RTPD1_vs_PD1", +1, "senescence after RT"),
    ("MP1_CELL_CYCLE_G2_M", "RTPD1_vs_PD1", -1, "G2/M down after RT"),
    ("MP2_CELL_CYCLE_G1_S", "RTPD1_vs_PD1", -1, "G1/S down after RT"),
], columns=["programme", "contrast", "expected_sign", "rationale"])

TREATMENT_CAVEAT = (
    "One biopsy per patient per timepoint: the timepoint cannot be separated from that "
    "biopsy's library. Responders lose tumour, so a late biopsy's malignant cells are the "
    "survivors - a change can be selection as well as reprogramming."
)

# --------------------------------------------------------------------------- #
# The embedding and the programmes
# --------------------------------------------------------------------------- #


def load_embedding():
    """05_3's embedding, and the cells x dimensions frame with columns TITLED ('DR 1' ...).

    Positional against `var['title']`, never `X[:, n]`: DRVI titles the dimensions by
    `reconstruction_effect`, so latent column n and `DR n` are different things. The same
    construction as `cell_first_tum.py` and `dr_umap_tum.py`.
    """
    import anndata as ad

    emb = C.get_embedding(C.DEFAULT_EMBEDDING)
    if not emb.embed_h5ad.exists():
        sys.exit(f"missing {emb.embed_h5ad}: run 05_3_drvi_run/run_drvi_tum.py "
                 f"with N_LATENT/HVG_SET matching {C.RUN_ID}")
    adata = ad.read_h5ad(emb.embed_h5ad)
    for key in (BATCH_KEY, TREATMENT_KEY, RESPONSE_KEY):
        if key not in adata.obs:
            sys.exit(f"{emb.embed_h5ad.name} has no obs['{key}']")
    L = pd.DataFrame(np.asarray(adata.X), index=adata.obs_names,
                     columns=adata.var["title"].values)
    print(f"[read] {emb.embed_h5ad.name}: {adata.n_obs} cells x {adata.n_vars} dimensions")
    return adata, L


def read_programmes() -> pd.DataFrame:
    """The bars of `programmes_named`: one row per programme, its directions as a list.

    Read through `signature_common.read_table`, so the folder follows $OUT_TAG exactly as
    05_8 wrote it. A programme named out of two collections keeps both in `collections`.
    """
    try:
        prog = C.read_table("programme_dimensions", CONSENSUS)
    except FileNotFoundError:
        sys.exit(f"missing programme_dimensions_consensus_{C.RUN_ID}.csv under "
                 f"{C.TABLE_DIR.name}/: run 05_8_convergence/dr_naming_tum.py with the same "
                 f"PRUNE_VANISHED / OUT_TAG first")
    prog = prog.copy()
    prog["directions"] = [[d.strip() for d in str(s).split(",") if d.strip()]
                          for s in prog["dimensions"]]
    prog["collections"] = [[c.strip() for c in str(s).split(",") if c.strip()]
                           for s in prog["collections"]]
    print(f"[read] {len(prog)} programmes on "
          f"{len({d for ds in prog['directions'] for d in ds})} directions "
          f"({C.TABLE_DIR.name}/consensus/{C.RUN_ID})")
    return prog


def direction_parts(direction: str) -> tuple[str, float]:
    """'DR 18-' -> ('DR 18', -1.0)."""
    m = re.fullmatch(r"\s*(DR \d+)\s*([+-])\s*", direction)
    if not m:
        raise ValueError(f"not a dimension-direction: {direction!r}")
    return m.group(1), (1.0 if m.group(2) == "+" else -1.0)


def oriented_directions(L: pd.DataFrame, directions) -> pd.DataFrame:
    """Cells x directions: the coordinate times the sign of the direction, z-scored over cells.

    ORIENTED because a name is a claim about one side of an axis - `DR 18-` names the
    negative side - so "higher" has to mean "more of the programme" in every column.
    Z-SCORED over all cells so a paired delta reads in one unit for every direction: standard
    deviations of the whole malignant population on that axis. The centre is irrelevant to a
    paired difference; the scale is what makes two directions comparable.
    """
    out = {}
    for d in directions:
        dim, sign = direction_parts(d)
        if dim not in L.columns:
            sys.exit(f"{d}: no dimension titled {dim!r} in the embedding")
        v = sign * L[dim].to_numpy(float)
        sd = v.std()
        out[d] = (v - v.mean()) / (sd if sd > 0 else 1.0)
    return pd.DataFrame(out, index=L.index)


def programme_composites(Z: pd.DataFrame, prog: pd.DataFrame) -> pd.DataFrame:
    """Cells x programmes: the mean of the programme's oriented, z-scored directions.

    A programme named on several directions (`HALLMARK_IFNG` on four) is ONE hypothesis, so
    it gets one readout; the directions stay available one by one in the exploratory family.
    """
    return pd.DataFrame({p: Z[r["directions"]].mean(axis=1) for p, r in prog.iterrows()},
                        index=Z.index)


# --------------------------------------------------------------------------- #
# Samples
# --------------------------------------------------------------------------- #


def sample_ids(obs: pd.DataFrame) -> pd.Series:
    """`Patient01__BASE`: one biopsy. Double underscore, safe as an R column name."""
    return obs[BATCH_KEY].astype(str) + "__" + obs[TREATMENT_KEY].astype(str)


def sample_table(obs: pd.DataFrame, min_cells: int) -> pd.DataFrame:
    """One row per biopsy: patient, timepoint, response, malignant cells, kept or not."""
    s = (obs.assign(sample=sample_ids(obs))
         .groupby("sample", observed=True)
         .agg(cohort=(BATCH_KEY, "first"), treatment=(TREATMENT_KEY, "first"),
              response=(RESPONSE_KEY, "first"), n_cells=(BATCH_KEY, "size")))
    for key in ("cohort", "treatment", "response"):
        s[key] = s[key].astype(str)
    unknown = set(s["treatment"]) - set(TREATMENTS)
    if unknown:
        sys.exit(f"unexpected treatment level(s) {sorted(unknown)}; expected {TREATMENTS}")
    s["kept"] = s["n_cells"] >= min_cells
    return s.sort_values(["response", "cohort", "treatment"])


def count_matrix(samples: pd.DataFrame) -> pd.DataFrame:
    """Patients x timepoints, malignant cells per biopsy, response as the first column."""
    m = samples.pivot_table(index="cohort", columns="treatment", values="n_cells",
                            aggfunc="sum", fill_value=0)
    m = m.reindex(columns=list(TREATMENTS), fill_value=0)
    resp = samples.groupby("cohort")["response"].first()
    m.insert(0, "response", resp.loc[m.index])
    m = m.rename_axis("cohort").reset_index()
    return m.sort_values(["response", "cohort"]).set_index("cohort")


def report_samples(samples: pd.DataFrame, min_cells: int) -> None:
    kept = samples[samples["kept"]]
    dropped = samples[~samples["kept"]]
    print(f"{len(samples)} biopsies, {len(kept)} with >= {min_cells} malignant cells "
          f"({kept['n_cells'].sum()} of {samples['n_cells'].sum()} cells)")
    if len(dropped):
        print(f"dropped ({len(dropped)}): " + ", ".join(
            f"{r.cohort} {r.treatment} ({r.n_cells})" for r in dropped.itertuples()))
    per_patient = kept.groupby("cohort")["treatment"].apply(lambda t: "+".join(
        x for x in TREATMENTS if x in set(t)))
    print("timepoints per patient after the floor: " + ", ".join(
        f"{k} {v}" for k, v in per_patient.value_counts().items()))


def bh(p) -> np.ndarray:
    """Benjamini-Hochberg over the non-missing entries; NaN stays NaN."""
    p = np.asarray(p, dtype=float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if ok.sum():
        from statsmodels.stats.multitest import multipletests
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


def diverging_cmap():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("treatment_diverging", DIVERGING)


def collection_by_slug(slug: str):
    """A collection by the slug its outputs are written under - `gavish_tnbc` included."""
    for coll in (*SC.COLLECTIONS.values(), *SC.GAVISH_VARIANTS.values()):
        if coll.name == slug:
            return coll
    raise KeyError(f"no collection writes under the slug {slug!r}")


def footnote(fig) -> None:
    """The treatment caveat, just under the phase's own two lines (`savefig` stamps those at
    y=-0.02, 6 pt). The gap is fixed in inches and converted, so it holds at any figure height.
    """
    y = -0.02 - 0.30 / fig.get_figheight()
    fig.text(0.5, y, TREATMENT_CAVEAT, ha="center", va="top", fontsize=6, style="italic",
             color="0.35", wrap=True)
