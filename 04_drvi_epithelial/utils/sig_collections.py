"""The signature collections this phase interprets, and everything that differs between them.

04_3 - 04_7 implement ONE procedure: score prior knowledge on the cells (Route A), read the
gene programme off the DRVI decoder (Route B), and call a cell state only where the two
converge (Route C). That procedure is applied twice, to two independent bodies of prior
knowledge:

  * `scie`  - stemness x immunogenicity. Is there an epithelial state that is stem-like AND
              immune-evasive? Ten lab lists on two axes, plus CytoTRACE2 as a stemness
              readout independent of any list.
  * `emt`   - the EMT axis. Which cells sit in the HYBRID (partial-EMT) state, the one the
              hysteresis project is about? Three lists per axis version, on the epithelial /
              hybrid / mesenchymal axes, plus a derived E-to-M score per version. The hybrid
              state is called by CO-EXPRESSION of the epithelial and mesenchymal programmes;
              the hybrid lists validate that call rather than making it.

Nothing about the procedure changes between them. What changes is the input lists, the axis
names, the shape of the target region on the cell-first plane, and which named risk has to be
checked before the result is believed. All of that lives here, so the five step scripts stay
single files taking `--collection {scie,emt}` rather than being duplicated per readout.

The outputs never mix: every table and figure is written to `<tables|figures>/<collection>/`
and carries the collection in its filename. The Benjamini-Hochberg correction of 04_6 is
likewise computed inside a collection, so adding the EMT lists cannot move a single SCIE
p-value.

Adding a third collection means appending a `Collection` below. No step script changes.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable


# --------------------------------------------------------------------------- #
# The pieces a collection is made of
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Signature:
    """One gene list. `file` is the name in `$DATA_DIR/signatures/`, which does not always
    match the signature name; `provenance` goes verbatim into the .gmt description field and
    doubles as the Appendix table.

    `primary=False` marks a list kept as a robustness check on a primary one rather than as an
    independent test. It is scored and reported like any other, but it is there to show that a
    result does not depend on which version of the list was used - and it is why the stability
    table across list versions is a reported result, not a diagnostic. Which KIND of check it is
    matters and is stated in `provenance`: a strict subset of the primary list is a sensitivity
    analysis, only a separately curated list is a replicate. See the EMT block below.
    """

    name: str
    file: str
    axis: str
    provenance: str
    primary: bool = True


@dataclass(frozen=True)
class Derived:
    """A readout computed from two z-scored signatures as `z(plus) - z(minus)`.

    It has no Route B counterpart by construction: it is a per-cell contrast, not a gene set,
    so 04_6 never sees it and 04_7 reports it as an independent Route A readout. The subtraction
    is done AFTER the within-stratum standardisation, never on the raw scores, so the two halves
    are on the same scale before they are differenced.
    """

    name: str
    plus: str
    minus: str
    axis: str
    description: str


@dataclass(frozen=True)
class Plane:
    """One panel of the cell-first plane: two readouts and the rule that defines the target.

    `rule` is one of:
      * "high" - at or above the high quantile
      * "low"  - at or below the low quantile
      * "mid"  - inside the middle band, i.e. NEITHER end of the axis

    "mid" describes the middle of an axis rather than an end of it. No collection defines its
    target with it today - see the EMT plane below for why the co-expression form was preferred -
    but it stays available, with `--mid-lo-q` / `--mid-hi-q`, because the partial-EMT state was
    first written that way and the comparison is in the record.
    """

    label: str
    x: str
    y: str
    x_rule: str
    y_rule: str


@dataclass(frozen=True)
class Criterion:
    """One half of the "is this axis the project's target?" test in 04_7.

    The dimension is scored by the strongest Route A correlation over `axis` (every readout on
    that axis) or over `names`, oriented by `sign`: +1 as is, -1 for a readout that has to be
    LOW in the target state, 0 for one whose magnitude matters but whose direction does not.
    """

    label: str
    sign: int
    axis: str | None = None
    names: tuple[str, ...] = ()


@dataclass(frozen=True)
class Collection:
    name: str                       # the slug in every path and filename
    title: str                      # what it is called in headings
    question: str                   # the biological question, printed by every step
    signatures: tuple[Signature, ...]
    axes: tuple[str, ...]           # column / block order in every heatmap
    planes: Callable[["Collection", list[str]], list[Plane]]
    plane_figure: str               # figure basename for the cell-first plane
    target_label: str               # what the target region is called in prose
    risks: tuple[str, ...]          # named checks 04_5 must run, see its header
    criteria: tuple[Criterion, ...]  # the 04_7 target-axis test
    derived: tuple[Derived, ...] = ()
    extra_readouts: dict[str, str] = field(default_factory=dict)  # name -> axis, joined at runtime
    depth_risk_readout: str | None = None      # the readout whose LOW group could just be shallow
    ambient_risk_axis: str | None = None       # the axis whose HIGH group could just be ambient
    extra_flags: Callable[["Collection", str, float], list[str]] = lambda c, n, r: []

    # ---------------------------------------------------------------- accessors

    @property
    def names(self) -> list[str]:
        return [s.name for s in self.signatures]

    @property
    def axis_of(self) -> dict[str, str]:
        """Every readout the collection can produce, including the derived and joined ones."""
        m = {s.name: s.axis for s in self.signatures}
        m.update({d.name: d.axis for d in self.derived})
        m.update(self.extra_readouts)
        return m

    @property
    def provenance(self) -> dict[str, str]:
        return {s.name: s.provenance for s in self.signatures}

    @property
    def files(self) -> dict[str, str]:
        return {s.name: s.file for s in self.signatures}

    def by_axis(self, axis: str) -> list[str]:
        """Every readout on one axis, registry order, derived and joined ones included."""
        return [n for n, a in self.axis_of.items() if a == axis]

    def primary_names(self) -> list[str]:
        return [s.name for s in self.signatures if s.primary]

    def order(self, present: list[str]) -> list[str]:
        """Reorder readouts by axis, then by registry order. The column order of every heatmap."""
        rank = {n: i for i, n in enumerate(self.names + [d.name for d in self.derived]
                                          + list(self.extra_readouts))}
        axis_rank = {a: i for i, a in enumerate(self.axes)}
        return sorted([n for n in present if n in self.axis_of],
                      key=lambda n: (axis_rank.get(self.axis_of[n], len(self.axes)),
                                     rank.get(n, len(rank))))

    def block_edges(self, ordered: list[str]) -> list[int]:
        """Positions where the axis changes, i.e. where a heatmap needs a separator line.

        The blocks are what stop a reader treating the columns as independent tests: the
        lists inside one block overlap heavily and are not ten (or nine) separate
        hypotheses, which is exactly what the Jaccard matrix of 04_3 is there to show.
        """
        edges, prev = [], None
        for i, n in enumerate(ordered):
            a = self.axis_of.get(n)
            if prev is not None and a != prev:
                edges.append(i)
            prev = a
        return edges


# --------------------------------------------------------------------------- #
# scie - stemness x immunogenicity
# --------------------------------------------------------------------------- #
#
# NOTE. There is deliberately NO stemness consensus signature. The intersection of
# BENPORATH_ES1 / ESC_WONG / ESC_ASSOU was tested by the lab and captures proliferation only,
# which the existing S_score / G2M_score already cover. Each stemness signature is used on its
# own and none of them is primary at this stage.

PRIMARY_IMMUNE = "IMMUNOGENIC_CONSENSUS"

_SCIE_SIGNATURES = (
    Signature("HALLMARK_IFNA",         "HALLMARK_IFNA.txt",      "immune",   "Interferon Alpha response, MSigDB Hallmark"),
    Signature("HALLMARK_IFNG",         "HALLMARK_IFNG.txt",      "immune",   "Interferon Gamma response, MSigDB Hallmark"),
    Signature("ISDS",                  "ISDS.txt",               "immune",   "IFN-Stem Cell-Down signature, PMC5481166"),
    Signature("KEGG_APM",              "KEGG_APM.txt",           "immune",   "Antigen Presentation Machinery, KEGG"),
    Signature("IMMUNOGENIC_CONSENSUS", "Immune consensus.txt",   "immune",   "Curated by the lab from the KEGG APM signature, retaining only immunogenic genes and immunomodulators; primary immune readout"),
    Signature("BENPORATH_ES1",         "BENPORATH_ES1.txt",      "stemness", "MSigDB BENPORATH_ES_1"),
    Signature("ESC_ASSOU",             "ESC_ASSOU.txt",          "stemness", "PMC1906587, Table S3"),
    Signature("ESC_WONG",              "ESC_WONG.txt",           "stemness", "MSigDB WONG_EMBRYONIC_STEM_CELL_CORE"),
    # EMP (Embryonic Multipotent Progenitors, PMID 29784918, `EMP.txt`) WAS the fourth
    # embryonic readout and was dropped. Not for coverage - 13 of its 15 symbols map, well
    # above the 60% floor - but because of what the 13 measure on THIS data. Three genes carry
    # 92% of the score: RPSA 64.8%, FN1 16.7%, MFAP2 10.6%. RPSA is a 40S ribosomal protein
    # and the score follows it (Spearman 0.70) rather than the progenitor markers the list is
    # named for; FN1 and MFAP2 are ECM, which is why the score tracked the mesenchymal axis of
    # the emt collection (rho 0.27-0.31) as closely as the other stemness lists (0.21-0.27).
    # The genes that make it an EMP list are below the droplet noise floor: NDNF detected in
    # 0.2% of cells, IGF2BP1 0.3%, FRAS1 1.1%, EPHA7 1.2%, SOX11 3.3% - and 20.9% of cells
    # have zero counts across all 13. The consequence is visible in the results it produced:
    # in 128 dimension-directions EMP was never significant on either route (max |rho| 0.104
    # against 0.30-0.35 for every other stemness readout; best-of-row 3x on A and 3x on B, all
    # non-significant), so it contributed no verdict, while costing stability - dropping it
    # raises the median pairwise Jaccard of the stemness quadrants from 0.184 to 0.255.
    # Restoring it is one line. Its own tables have been overwritten by the re-run, so these
    # numbers and the README section "Why EMP was dropped" are the record of what it did.
    Signature("LIM_STEM",              "LIM_STEM.txt",           "stemness", "MSigDB LIM_MAMMARY_STEM_CELL_UP"),
    Signature("FMASC",                 "fMaSC.txt",              "stemness", "Fetal mammary stem cells, PMC3277444, Supplementary Table 2"),
)


def _scie_planes(coll: "Collection", available: list[str]) -> list[Plane]:
    """One panel per stemness readout, all against the one primary immune readout.

    No stemness signature is primary, so the quadrant is defined once per stemness readout and
    the stability of the resulting cell set ACROSS those definitions is itself the result. The
    immune axis is fixed: the lab curated one primary immunogenic list and the other three are
    largely nested inside it.
    """
    stem = [n for n in coll.by_axis("stemness") if n in available]
    return [Plane(label=s, x=s, y=PRIMARY_IMMUNE, x_rule="high", y_rule="low") for s in stem]


def _scie_flags(coll: "Collection", claimed: str, a_rho: float) -> list[str]:
    # Immune evasion is defined by the ABSENCE of a signal, which shallow sequencing mimics
    # perfectly. A dimension claimed on a negative correlation with an immune list is exactly
    # that situation and is flagged wherever it is reported.
    if coll.axis_of.get(claimed) == "immune" and a_rho < 0:
        return ["immune_low_is_absence_of_signal"]
    return []


SCIE = Collection(
    name="scie",
    title="stemness x immunogenicity (SCIE)",
    question="is there an epithelial state that is stem-like AND immune-evasive?",
    signatures=_SCIE_SIGNATURES,
    axes=("immune", "stemness"),
    planes=_scie_planes,
    plane_figure="stemness_immunogenicity_plane",
    target_label="stem-high / immunogenic-low",
    risks=("cell_cycle", "depth"),
    criteria=(
        Criterion("stem_rho", sign=+1, axis="stemness"),
        # evasion = LOW immunogenicity, hence the sign
        Criterion("immunogenic_low_rho", sign=-1, names=(PRIMARY_IMMUNE,)),
    ),
    # CytoTRACE2 is joined by 04_5 from the csv 04_4 writes, when that csv exists. It is the
    # only stemness evidence in the whole stage that does not come from a gene list, so the
    # quadrant is materially weaker without it and 04_5 says so rather than failing quietly.
    extra_readouts={"CytoTRACE2": "stemness"},
    depth_risk_readout=PRIMARY_IMMUNE,
    extra_flags=_scie_flags,
)


# --------------------------------------------------------------------------- #
# emt - the epithelial-to-mesenchymal axis and its hybrid state
# --------------------------------------------------------------------------- #
#
# Three versions of the same three lists, from the collaborator (see
# `$DATA_DIR/signatures/EMT_LISTS_NOTES.md` for the symbol normalisation and for what was
# substituted). B is the primary triad. A (experimentally validated only) and C (Tomas' list,
# translated from mouse) are robustness checks on it, not independent tests: B vs C on the
# mesenchymal axis is Jaccard 0.76.
#
# The three are NESTED, and the wording above has to respect it. `A \ B` is empty on all three
# axes - A is a strict subset of B, and 13/15, 11/12 and 19/19 of it are also inside C. B in turn
# carries only 2 / 2 / 1 genes that no other list has, so B minus those genes IS the >=2-of-3
# majority-vote set (26 / 18 / 32 genes). That, and not its size, is why B is primary.
#
# It follows that A is NOT a replicate of B: a subset cannot replicate its superset. A vs B is a
# SENSITIVITY ANALYSIS on the non-validated genes - one variable changes, the validation status -
# while B vs C changes curator, species and size at once. Only C is a replicate.
#
# Nesting at the gene level is not redundancy at the readout level, which is why A is scored and
# not dropped. Scores are means over the set, so the extra genes move the score rather than
# extending it: A and B agree on only 32% of the cells they call (Jaccard 0.3215), LESS than B
# agrees with C (0.3851), and `EMT_A_MESENCHYMAL` is the strongest Route A correlate in 33 of the
# 128 dimension-directions against 8 for the `EMT_B_MESENCHYMAL` that contains it. Read that last
# fact next to `signature_concentration`: A_MESENCHYMAL has an effective n of 3.75 and VIM carries
# 38% of its variance (B: 6.92 and 26%, C: 9.50 and 21%), so the shorter list wins by being the
# purer VIM readout - which is the gene `mesenchymal_may_be_ambient_or_doublet` exists for, not a
# reason to promote it.
#
# The state of interest is the HYBRID one, and it is not an end of the axis. It is defined here
# by CO-EXPRESSION: a cell that is high on the epithelial programme AND high on the mesenchymal
# one at the same time. That is the textbook operationalisation of hybrid E/M, and on this
# dataset it is also the only one that survives its own robustness check.
#
# It was first written the other way - the E-to-M score inside a middle band, crossed with the
# hybrid list high - and that definition failed: the cell set it called had a Jaccard of only
# 0.08 - 0.13 across the three list versions. The diagnosis was not biological. The continuous
# scores agree well across versions (EMT_SCORE Spearman 0.72 - 0.87, epithelial 0.78 - 0.90);
# what does not survive is intersecting two NARROW quantile selections, because a 20%-wide band
# around the median of a continuous score is where the density is highest and where two
# correlated scores disagree most about rank. The instability was an artefact of the cutoffs.
#
# Co-expression fixes it on three counts, and the third is the one that matters:
#   * stability across list versions rises to 0.24 - 0.39 (median 0.32), 2.5x the old definition;
#   * it drops the hybrid list from the DEFINITION, which is the right way round - the hybrid
#     lists are the least reproducible of the three axes (Spearman 0.49 - 0.68 between versions,
#     against 0.78 - 0.90 for the epithelial ones) and the least experimentally settled;
#   * the middle band is then RECOVERED rather than assumed. The E-to-M score of the co-expression
#     set falls at percentiles 38 / 50 / 63 of the whole compartment without that ever having been
#     imposed, and the set is enriched 1.3 - 2.1x for hybrid-high cells. The hybrid lists became
#     the validation instead of the definition, which is what they are actually good for.

_EMT_SIGNATURES = (
    Signature("EMT_B_EPITHELIAL",  "EMT_B_EPITHELIAL.txt",  "epithelial",  "Collaborator list B (validated + non-validated), epithelial markers; primary triad, the >=2-of-3 consensus of A/B/C", primary=True),
    Signature("EMT_B_HYBRID",      "EMT_B_HYBRID.txt",      "hybrid",      "Collaborator list B, hybrid/partial-EMT markers; primary triad, the >=2-of-3 consensus of A/B/C", primary=True),
    Signature("EMT_B_MESENCHYMAL", "EMT_B_MESENCHYMAL.txt", "mesenchymal", "Collaborator list B, mesenchymal markers; primary triad, the >=2-of-3 consensus of A/B/C", primary=True),
    Signature("EMT_A_EPITHELIAL",  "EMT_A_EPITHELIAL.txt",  "epithelial",  "Collaborator list A (experimentally validated only), epithelial markers; strict subset of B, its validated core - a sensitivity analysis on the non-validated genes, not a replicate", primary=False),
    Signature("EMT_A_HYBRID",      "EMT_A_HYBRID.txt",      "hybrid",      "Collaborator list A, hybrid/partial-EMT markers; strict subset of B, its validated core - a sensitivity analysis, not a replicate", primary=False),
    Signature("EMT_A_MESENCHYMAL", "EMT_A_MESENCHYMAL.txt", "mesenchymal", "Collaborator list A, mesenchymal markers; strict subset of B, its validated core - a sensitivity analysis, not a replicate", primary=False),
    Signature("EMT_C_EPITHELIAL",  "EMT_C_EPITHELIAL.txt",  "epithelial",  "Collaborator list C (Tomas), epithelial markers, mouse symbols mapped to human; the one separately curated list, hence the only true robustness replicate of B", primary=False),
    Signature("EMT_C_HYBRID",      "EMT_C_HYBRID.txt",      "hybrid",      "Collaborator list C (Tomas), hybrid/partial-EMT markers, mouse symbols mapped to human; separately curated, the only true robustness replicate of B", primary=False),
    Signature("EMT_C_MESENCHYMAL", "EMT_C_MESENCHYMAL.txt", "mesenchymal", "Collaborator list C (Tomas), mesenchymal markers, mouse symbols mapped to human; separately curated, the only true robustness replicate of B", primary=False),
)

_EMT_DERIVED = tuple(
    Derived(
        name=f"EMT_SCORE_{v}",
        plus=f"EMT_{v}_MESENCHYMAL",
        minus=f"EMT_{v}_EPITHELIAL",
        axis="emt_score",
        description=f"z(EMT_{v}_MESENCHYMAL) - z(EMT_{v}_EPITHELIAL), the E-to-M position of a cell on list version {v}",
    )
    for v in ("B", "A", "C")
)


def _emt_planes(coll: "Collection", available: list[str]) -> list[Plane]:
    """One panel per list version: epithelial against mesenchymal, of the SAME version.

    The target is the top-right corner, i.e. both programmes high at once. This is the classic
    E/M plane, and reading it is the point: a clean, complete transition would leave that corner
    empty.

    Mixing versions across the two axes would confound "does the result depend on the list?"
    with "does it depend on the axis?", which is the whole point of carrying A and C.
    """
    out = []
    for v in ("B", "A", "C"):
        x, y = f"EMT_{v}_EPITHELIAL", f"EMT_{v}_MESENCHYMAL"
        if x in available and y in available:
            out.append(Plane(label=f"list {v}", x=x, y=y, x_rule="high", y_rule="high"))
    return out


def _emt_flags(coll: "Collection", claimed: str, a_rho: float) -> list[str]:
    # In an epithelial compartment, high VIM / FN1 / SPARC / ACTA2 is as easily fibroblast
    # ambient RNA or an epithelial-fibroblast doublet as it is a transition. Anything claimed
    # on the mesenchymal axis carries the flag wherever it is reported; 04_5 runs the actual
    # doublet check behind it.
    if coll.axis_of.get(claimed) == "mesenchymal":
        return ["mesenchymal_may_be_ambient_or_doublet"]
    return []


EMT = Collection(
    name="emt",
    title="epithelial-mesenchymal transition (EMT)",
    question="which cells sit in the hybrid, partial-EMT state?",
    signatures=_EMT_SIGNATURES,
    axes=("epithelial", "hybrid", "mesenchymal", "emt_score"),
    planes=_emt_planes,
    plane_figure="emt_coexpression_plane",
    target_label="hybrid state: epithelial AND mesenchymal programmes both high (co-expression)",
    risks=("cell_cycle", "ambient"),
    criteria=(
        Criterion("hybrid_rho", sign=+1, axis="hybrid"),
        # An axis is an E-to-M axis whichever way round it is oriented, so magnitude only.
        Criterion("emt_axis_rho", sign=0, axis="emt_score"),
    ),
    derived=_EMT_DERIVED,
    ambient_risk_axis="mesenchymal",
    extra_flags=_emt_flags,
)


# --------------------------------------------------------------------------- #

COLLECTIONS = {c.name: c for c in (SCIE, EMT)}
DEFAULT_COLLECTION = "scie"


def get(name: str) -> Collection:
    if name not in COLLECTIONS:
        raise KeyError(f"unknown collection {name!r}; have: {', '.join(COLLECTIONS)}")
    return COLLECTIONS[name]


def add_argument(parser) -> None:
    """The `--collection` flag, spelled identically by all five step scripts."""
    parser.add_argument("--collection", choices=sorted(COLLECTIONS), default=DEFAULT_COLLECTION,
                        help=f"which signature collection to interpret (default {DEFAULT_COLLECTION})")
