"""The signature collections this phase interprets, and everything that differs between them.

04_3 - 04_7 implement ONE procedure: score prior knowledge on the cells (Route A), read the
gene programme off the DRVI decoder (Route B), and call a cell state only where the two
converge (Route C). That procedure is applied to three independent bodies of prior knowledge,
the third of which uses only the part of it that does not need a target region:

  * `scie`  - stemness x immunogenicity. Is there an epithelial state that is stem-like AND
              immune-evasive? Ten lab lists on two axes, plus CytoTRACE2 as a stemness
              readout independent of any list.
  * `emt`   - the EMT axis. Which cells sit in the HYBRID (partial-EMT) state, the one the
              hysteresis project is about? Three lists per axis version, on the epithelial /
              hybrid / mesenchymal axes, plus a derived E-to-M score per version. The hybrid
              state is called by CO-EXPRESSION of the epithelial and mesenchymal programmes;
              the hybrid lists validate that call rather than making it.
  * `gavish` - the recurrent pan-cancer metaprograms of Gavish et al. 2023, used as a
              VOCABULARY rather than as a hypothesis: it defines no target region and calls no
              cells, it names the latent dimensions from outside this dataset and gives the
              other two collections an independent check. It comes in two widths: the 22
              relevant to a triple-negative breast carcinoma, which is what it scores by
              default, and all 41 with `--all-metaprograms`. The wider one also carries the
              lineages that cannot be here, as negative controls; the default no longer does.

Nothing about the procedure changes between them. What changes is the input lists, the axis
names, the shape of the target region on the cell-first plane - or whether there is one at all,
see `Collection.has_target` - and which named risk has to be checked before the result is
believed. All of that lives here, so the step scripts stay single files taking
`--collection {scie,emt,gavish}` rather than being duplicated per readout.

The outputs never mix: every table and figure is written to `<tables|figures>/<collection>/`
and carries the collection in its filename. The Benjamini-Hochberg correction of 04_6 is
likewise computed inside a collection, so adding the EMT lists cannot move a single SCIE
p-value, and neither can adding the forty-one Gavish metaprograms.

Adding a collection means appending a `Collection` below and nothing else - with one exception,
now on the record: a collection WITHOUT a target region needed the steps taught to skip what
depends on one. That is `has_target`, it was added for `gavish`, and it is now part of the
contract rather than a special case. The two steps that exist only to interrogate a target
region, 04_8 and 04_9, stop with a message on such a collection instead of running.
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

    A NEGATIVE CONTROL is the third kind, and it is the one the flag does not name: a list that
    is scored so that the OTHER results have a floor to be read against, and whose expected
    value is zero. `primary=False` marks those too, for the same reason - it is not an
    independent test of anything - and `provenance` again says which kind it is. See the
    lineage metaprograms of the GAVISH block.
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
    # THE TARGET REGION, AND WHAT A COLLECTION IS WITHOUT ONE.
    #
    # A collection that asks whether a NAMED STATE EXISTS defines all five: the planes that
    # define the region, what to call it, the named risks 04_5 must check, and the
    # criteria 04_7 tests the axes against. `scie` and `emt` do.
    #
    # A collection can instead be a VOCABULARY: a body of prior knowledge used to NAME the
    # dimensions rather than to call cells. It leaves all five at their defaults, `has_target`
    # is then False, and each step drops exactly what depends on a region - the cell-first
    # plane, the consensus quadrant and everything computed on it, the target-axis test - and
    # keeps the two heatmaps, the per-cell scores and the convergence table, which need no
    # region to mean anything. `gavish` is one; see its block at the bottom of this file.
    planes: Callable[["Collection", list[str]], list[Plane]] | None = None
    plane_figure: str | None = None  # figure basename for the cell-first plane
    target_label: str | None = None  # what the target region is called in prose
    risks: tuple[str, ...] = ()      # named checks 04_5 must run, see its header
    criteria: tuple[Criterion, ...] = ()  # the 04_7 target-axis test
    derived: tuple[Derived, ...] = ()
    extra_readouts: dict[str, str] = field(default_factory=dict)  # name -> axis, joined at runtime
    depth_risk_readout: str | None = None      # the readout whose LOW group could just be shallow
    ambient_risk_axis: str | None = None       # the axis whose HIGH group could just be ambient
    extra_flags: Callable[["Collection", str, float], list[str]] = lambda c, n, r: []

    # ---------------------------------------------------------------- accessors

    @property
    def has_target(self) -> bool:
        """Does this collection define a region of the cell-first plane to call cells in?

        The single question every step asks before running anything downstream of a quadrant.
        It is derived from `planes` rather than stored, so a collection cannot end up claiming
        a target it has no way to draw.
        """
        return self.planes is not None

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
# gavish - the recurrent pan-cancer metaprograms, as a vocabulary
# --------------------------------------------------------------------------- #
#
# Gavish et al. 2023 (Nature 618:598-606, "Hallmarks of transcriptional intratumour
# heterogeneity across a thousand tumours") ran NMF on the malignant cells of ~1,000 tumours
# across 24 cancer types and clustered the resulting programmes into 41 METAPROGRAMS: the
# expression programmes that recur across patients and across cancer types rather than in a
# single tumour. Each is published as a list of 50 genes.
#
# THIS COLLECTION IS A VOCABULARY, NOT A HYPOTHESIS, and that is the whole difference between
# it and the two above. `scie` and `emt` each ask whether a NAMED state exists, and each
# defines a region of the cell-first plane to call cells in. This one asks the question in the
# opposite direction - given a latent dimension, what is it? - so it defines no plane, no
# target region, no named risks and no criteria: `has_target` is False and every step skips
# what depends on them. What it does produce is the half that a latent space actually needs:
#
#   * Route A  every metaprogram scored per cell, standardised within cohort, correlated
#              against every dimension -> the dimensions x metaprograms heatmap;
#   * Route B  ORA of each dimension's top decoder genes against the same sets - the
#              standard way NMF and latent programmes are named in this literature;
#   * Route C  the per-dimension convergence table, which is the point of running it: a
#              dimension on which the two routes independently land on the SAME metaprogram
#              is a dimension that has a name coming from outside this dataset.
#
# It is also the outside check on the other two collections. `scie` and `emt` are lab and
# collaborator lists, chosen because the project is about those states; MP12 - MP16 and
# MP17 - MP18 were derived with no knowledge of this project, from other tumours and other
# cancer types. An EMT axis that appears on the collaborator's lists AND on Gavish's EMT
# metaprograms is an axis that does not depend on whose EMT list was used.
#
# WHAT IS MISSING, AND IT IS NOT MISSING HERE. Eight of the forty columns of
# `datasets/GAVISH.csv` carry 48 or 49 genes rather than 50, which is why the provenance below
# says "as exported" and why the coverage table is the place to read the per-list count.
#
# MP1 IS THE ONE LIST THAT DOES NOT COME FROM THAT CSV. The export has 40 columns and not 41:
# MP1 (Cell Cycle - G2/M) is absent from it, and from the MSigDB release the CSV was taken
# from - `GAVISH_3CA_MALIGNANT_METAPROGRAM_1_CELL_CYCLE_G2_M` does not exist there while every
# other number does. That left a hole in the vocabulary rather than in the export: a dimension
# that is specifically G2/M had nothing to match and read as MP2 or as nothing at all, on a
# compartment where proliferation is one of the states most likely to take a dimension of its
# own. It is now on disk, taken from the authors' own object - `MP_list.RDS` of
# github.com/tiroshlab/3ca, `ITH_hallmarks/MPs_distribution/`, entry `$Cancer[[1]]`, 50 genes
# in the paper's order - and written by `utils/gavish_extraction.py` alongside the forty from
# the CSV. The two sources are the same list: the MP2 of that object is identical
# gene-for-gene to `MP2_CELL_CYCLE_G1_S.txt` apart from two symbols MSigDB updated
# (HIST1H4C -> H4C3, KIAA0101 -> PCLAF).
#
# MP1 KEEPS THE PAPER'S SPELLING, AND THAT IS WHY ALL 50 OF ITS GENES MAP. Applying the same
# update would rename HIST1H4C to H4C3, and this dataset is on an older reference: it has
# HIST1H4C and not H4C3, so the updated spelling would cost MP1 a gene - as it already costs
# MP2 and MP3 one each. The coverage table is where that is visible for every list.
#
# THE NAMES ARE THE PAPER'S, THE AXES ARE OURS. Gavish numbers and names the metaprograms;
# grouping them into the fifteen families below is a decision taken here, and it is not
# cosmetic. `axis` is what orders and blocks every heatmap AND what Route C reads as
# `same_family`, so two metaprograms placed on one axis are two metaprograms whose agreement
# across the routes will be counted as convergence. The families are therefore drawn
# conservatively - MP4 (chromatin) is not folded into the cycle, MP38 / MP39 are not folded
# into metabolism - and a metaprogram with no relative sits on an axis of its own rather than
# in a bin of leftovers.
#
# THE LINEAGE METAPROGRAMS ARE KEPT HERE, AS NEGATIVE CONTROLS. Eleven of the 41 describe lineages
# that cannot be in a breast epithelial or malignant-epithelial compartment: the neural five
# (MP25 - MP29, glioma and oligodendrocyte), skin pigmentation (MP32), and the haematopoietic
# five (MP33 - MP37, erythrocytes, platelets, immunoglobulin). They are marked
# `primary=False` and they are not there to be found. They are there because forty correlated
# scores need a floor: if a dimension lands on MP25 (astrocytes) as strongly as on MP12 (EMT),
# then the MP12 reading is worth nothing either, and without the controls in the same table
# there is nothing to say so. Their expected result is a near-zero row; a non-zero one is a
# finding about the scoring rather than about the tumour - MP36 (IG) in particular is this
# dataset's ambient-immunoglobulin readout, and it is the one to look at first.
#
# All of that is true of THIS collection, the full one. The TNBC subset below carries none of
# them any more; the block that defines it says what that costs and how to get the floor back.

_GAVISH_PAPER = "Gavish et al. 2023, Nature 618:598-606, pan-cancer malignant metaprogram"

# The two families that cannot be present in this compartment. Membership of one is what
# makes a metaprogram a control rather than a test, in one place.
GAVISH_CONTROL_AXES = ("neural", "other_lineage")


def _mp(stem: str, axis: str, label: str, note: str = "", name: str | None = None,
        source: str = "the published gene list as exported to datasets/GAVISH.csv") -> Signature:
    """One metaprogram. `stem` is the filename `utils/gavish_extraction.py` wrote.

    The files live in a subdirectory of `$DATA_DIR/signatures/`, which `Signature.file`
    carries verbatim - `load_signatures` joins it to the signature directory and never
    assumed a flat layout.

    `name` overrides the default only for MP3, whose file is `MP3_CELL_CYLCE_HMG_RICH.txt`:
    the typo comes from the column header of GAVISH.csv and renaming the file would break a
    re-export rather than fix anything, so the readout carries the corrected spelling and the
    path keeps what is on disk.
    """
    # "as exported" and not "the 50 published genes": the metaprograms are 50 genes each in
    # the paper, and eight of the columns of GAVISH.csv carry 48 or 49. The per-list count is
    # in the coverage table of 04_3 / 05_4, which is where it belongs.
    provenance = f"{_GAVISH_PAPER} MP{stem.split('_')[0][2:]} ({label}), {source}"
    if note:
        provenance += f"; {note}"
    return Signature(name=name or stem, file=f"GAVISH_metaprograms/{stem}.txt", axis=axis,
                     provenance=provenance, primary=axis not in GAVISH_CONTROL_AXES)


_CONTROL_NOTE = ("lineage absent from a breast epithelial compartment by construction - kept "
                 "as an internal negative control, its expected correlation is zero")

_GAVISH_SIGNATURES = (
    _mp("MP1_CELL_CYCLE_G2_M",           "cell_cycle",         "Cell Cycle - G2/M",
        source="the published gene list, taken from the authors' MP_list.RDS rather than from "
               "GAVISH.csv, which does not carry it - see the note above"),
    _mp("MP2_CELL_CYCLE_G1_S",           "cell_cycle",         "Cell Cycle - G1/S"),
    _mp("MP3_CELL_CYLCE_HMG_RICH",       "cell_cycle",         "Cell Cycle - HMG-rich",
        name="MP3_CELL_CYCLE_HMG_RICH"),
    _mp("MP4_CHROMATIN",                 "chromatin",          "Chromatin"),
    _mp("MP5_STRESS",                    "stress",             "Stress"),
    _mp("MP6_HYPOXIA",                   "stress",             "Hypoxia"),
    _mp("MP7_STRESS_IN_VITRO",           "stress",             "Stress (in vitro)",
        note="derived from cultured cells, so a high score on fresh tissue is a dissociation "
             "readout before it is a biological one"),
    _mp("MP8_PROTEASOMAL_DEGRADATION",   "proteostasis",       "Proteasomal degradation"),
    _mp("MP9_UNFOLDED_PROTEIN_RESPONSE", "proteostasis",       "Unfolded protein response"),
    _mp("MP10_PROTEIN_MATURATION",       "proteostasis",       "Protein maturation"),
    _mp("MP11_TRANSLATION_INITIATION",   "proteostasis",       "Translation initiation",
        note="ribosomal-protein heavy, so it tracks library complexity as readily as a "
             "programme - read it next to rho_n_genes_by_counts in the confounder table"),
    _mp("MP12_EMT_1",                    "emt",                "EMT-I"),
    _mp("MP13_EMT_2",                    "emt",                "EMT-II"),
    _mp("MP14_EMT_3",                    "emt",                "EMT-III"),
    _mp("MP15_EMT_4",                    "emt",                "EMT-IV"),
    _mp("MP16_MES_GLIOMA",               "emt",                "MES (glioma)",
        note="the glioma mesenchymal programme, on the EMT axis because it is the "
             "mesenchymal one, not because gliomas are expected here"),
    _mp("MP17_INTERFERON_MHC_II_1",      "immune",             "Interferon/MHC-II (I)"),
    _mp("MP18_INTERFERON_MHC_II_2",      "immune",             "Interferon/MHC-II (II)"),
    _mp("MP19_EPITHELIAL_SENESCENCE",    "senescence",         "Epithelial Senescence"),
    _mp("MP20_MYC",                      "metabolic",          "MYC"),
    _mp("MP21_RESPIRATION",              "metabolic",          "Respiration"),
    _mp("MP38_GLUTATHIONE",              "detoxification",     "Glutathione"),
    _mp("MP39_METAL_RESPONSE",           "detoxification",     "Metal-response"),
    _mp("MP22_SECRETED_1",               "secreted",           "Secreted-I"),
    _mp("MP23_SECRETED_2",               "secreted",           "Secreted-II"),
    _mp("MP24_CILIA",                    "cilia",              "Cilia"),
    _mp("MP25_ASTROCYTES",               "neural",             "Astrocytes",           _CONTROL_NOTE),
    _mp("MP26_NPC_GLIOMA",               "neural",             "NPC (glioma)",         _CONTROL_NOTE),
    _mp("MP27_OLIGO_PROGENITOR",         "neural",             "Oligo Progenitor",     _CONTROL_NOTE),
    _mp("MP28_OLIGO_NORMAL",             "neural",             "Oligo normal",         _CONTROL_NOTE),
    _mp("MP29_NPC_OPC",                  "neural",             "NPC/OPC",              _CONTROL_NOTE),
    _mp("MP30_PDAC_CLASSICAL",           "epithelial_lineage", "PDAC-classical",
        note="a carcinoma epithelial-identity programme from another organ: not a control, "
             "an out-of-tissue epithelial reference"),
    _mp("MP31_ALVEOLAR",                 "epithelial_lineage", "Alveolar",
        note="lung epithelial identity; an out-of-tissue epithelial reference, as MP30"),
    _mp("MP40_PDAC_RELATED",             "epithelial_lineage", "PDAC-related",
        note="an out-of-tissue epithelial reference, as MP30"),
    _mp("MP32_SKIN_PIGMENTATION",        "other_lineage",      "Skin-pigmentation",    _CONTROL_NOTE),
    _mp("MP33_RBCS",                     "other_lineage",      "RBCs",                 _CONTROL_NOTE),
    _mp("MP34_PLATELET_ACTIVATION",      "other_lineage",      "Platelet activation",  _CONTROL_NOTE),
    _mp("MP35_HEMATO_RELATED_1",         "other_lineage",      "Hemato-related-I",     _CONTROL_NOTE),
    _mp("MP36_IG",                       "other_lineage",      "IG",
        note="immunoglobulin; the ambient-RNA readout of this dataset as much as a control - "
             "plasma-cell transcripts are the classic soup, so read it first"),
    _mp("MP37_HEMATO_RELATED_2",         "other_lineage",      "Hemato-related-II",    _CONTROL_NOTE),
    _mp("MP41_UNASSIGNED",               "unassigned",         "Unassigned",
        note="the paper's own residual metaprogram: it has no interpretation to lend, and it "
             "is here so that a dimension matching nothing else can match it instead"),
)

# The column order of every gavish heatmap, kept in one place because the two variants below
# share it: each takes the families it actually has, in this order.
_GAVISH_AXES = ("cell_cycle", "chromatin", "stress", "proteostasis", "emt", "immune",
                "senescence", "metabolic", "detoxification", "secreted", "cilia", "neural",
                "epithelial_lineage", "other_lineage", "unassigned")


# --------------------------------------------------------------------------- #
# The TNBC-relevant subset, which is what `--collection gavish` scores by default
# --------------------------------------------------------------------------- #
#
# WHY A SUBSET AT ALL. Forty-one scores against thirty-two or sixty-four dimension-directions
# is not a free vocabulary: every metaprogram is a column of both heatmaps, a row of the
# coverage and Jaccard tables, and - this is the part that costs something - one more test
# inside the Benjamini-Hochberg correction of 04_6 / 05_7. Nineteen of the forty-one name a
# lineage or a tissue that a triple-negative breast carcinoma cannot express, or were measured
# to carry nothing here, and they are spending that budget to confirm what is already known.
# Scoring the states a TNBC malignant epithelial cell can actually be in leaves the same result
# better resolved and the figures readable.
#
# WHAT IS IN, AND ON WHAT GROUND. Twenty-one metaprograms, every one of them a state that has
# been reported in breast - and in most cases specifically in basal-like / triple-negative -
# malignant cells: the whole cycle (MP1 G2/M, MP2 G1/S, MP3 HMG-rich), chromatin (MP4), stress
# and hypoxia (MP5, MP6), the proteostasis block (MP8 - MP10), all four EMT programmes
# (MP12 - MP15, the axis this project is about), interferon / MHC-II (MP17, MP18, the
# immune-visibility axis `scie` asks about from the other side), epithelial senescence (MP19),
# MYC (MP20, the classic basal-like amplification), respiration (MP21), the two secreted
# programmes (MP22, MP23), and metal-response (MP39, a chemoresistance programme in this
# disease).
#
# PLUS ONE THAT IS ONLY THERE TO BOUND THEM. MP7 (Stress in vitro) is the specificity control
# on the stress axis. MP5 and MP6 are scored, and the honest question about any dimension that
# matches them is whether it is in-vivo stress or dissociation. MP7 was derived in culture; a
# dimension that matches MP5 and MP7 equally is answering that question the wrong way. It is
# the one non-state list the subset carries.
#
# WHAT IS OUT. Cilia (MP24), which breast epithelium does not have; the glioma mesenchymal
# programme (MP16), whose EMT-like genes are already covered four times over by MP12 - MP15;
# the out-of-tissue epithelial identities (MP30, MP31, MP40); the neural five (MP25 - MP29);
# skin pigmentation (MP32); the haematopoietic five (MP33 - MP37); and the paper's own
# residual (MP41). That is every lineage control, and it is a decision taken by reading the
# full 41 and naming the ones that do not belong on this compartment - not the default of an
# earlier version of this file, which kept five of them.
#
# AND TWO THAT WERE MEASURED OUT RATHER THAN ARGUED OUT: MP38 and MP11. Both were in the set
# on the reasoning above - MP38 as a chemoresistance programme, MP11 as proteostasis - and
# both fail against a null. The null is 200 random 50-gene lists, matched bin-for-bin to the
# expression profile of the real metaprograms, scored with the Route A settings of 04_5 / 05_6
# and correlated against the same 64 dimension-directions of `drvi_tum_64_nomt`; its max |rho|
# has p50 0.225, p95 0.325, p99 0.372 and never exceeded 0.432 in 200 draws.
#
#   * MP38 (Glutathione) reaches |rho| 0.165, the 25th percentile of that null: a random list
#     beats it three times in four. Its genes say why - ACSM2A/B, AGXT2, CUBN, AMN, CLTRN,
#     FOLR1, AQP1 are proximal-tubule renal, so this is an out-of-tissue lineage wearing the
#     name of a metabolic programme, and it belongs with MP30 / MP31 / MP40.
#   * MP11 (Translation initiation) reaches 0.216, the 42nd percentile - nothing on Route A -
#     while firing 18 significant pairs on Route B. That split is the finding: EIF2/EIF3
#     genes land in the top decoder genes of several dimensions without the per-cell score
#     tracking any of them, which is what a library-complexity artefact looks like from both
#     sides. It was already flagged for this in its note below; the null settles it.
#
# The null is NOT in this file and NOT a threshold applied anywhere in the pipeline. It is a
# Methods number, and the two lists above are the only decision taken with it. Note also what
# it does NOT justify: the number of genes a list has inside the 2,000-HVG panel predicts
# nothing about whether it can be found. MP21 has 2 of them and fires 23 times on Route B;
# MP2 has 28 and fires 6. `MIN_SIGNATURE_GENES` is a floor on the MAPPED count, never on the
# HVG one, and the HVG count is reported as a warning for exactly this reason.
#
# WHAT DROPPING THE CONTROLS COSTS, WRITTEN DOWN SO THAT IT IS NOT REDISCOVERED. The block
# above says why the full collection keeps them: forty correlated scores need a floor, or the
# ones that look large have nothing to be large against. This variant no longer has that
# floor, and three things follow from it:
#
#   * the row that says what a correlation of nothing looks like on this data, at this depth,
#     is gone. `--all-metaprograms` still has it, and that is the run to do once per latent
#     space rather than never - it is the same tables under the other slug, so it costs a
#     command and overwrites nothing.
#   * the two ambient-RNA readouts go with it, MP36 (IG) and MP33 (RBCs). This is the least
#     costly of the three on this dataset: the 05 diagnosis put the ambient contribution about
#     two orders of magnitude below the dimensions it could have explained, and SoupX in phase
#     06 is a check on the same question that does not go through a metaprogram at all.
#   * MP41 (Unassigned) was the sink that let a dimension matching no real programme match it
#     rather than be pushed onto the nearest one. Without it, "this dimension matches nothing"
#     has to be read off the effect size and the significance columns instead of off a
#     competing row.
#
# `primary=False` and the `lineage_control_claim_bounds_the_rest` flag therefore never fire on
# this variant. They are not dead code - they fire on `gavish`, which is the same machinery on
# the same tables.
#
# CHANGING THE SET IS THE ONE LINE BELOW. The first candidates to add back are MP30 and MP40:
# they are pancreatic in name only, their genes are the generic secretory-epithelial ones
# (TFF, AGR2, mucins, CEACAM), and those are luminal-breast genes - so a dimension carrying
# epithelial identity rather than EMT has somewhere to land if they are in, and reads as
# "matches nothing" if they are out. They were left out because this compartment is
# basal-like, not because the reading would be wrong. After them, MP16 and MP24: they were
# out before the review that removed the controls and were not named by it either way.
#
# THE TWO VARIANTS NEVER OVERWRITE EACH OTHER. They are two collections with two slugs -
# `gavish_tnbc` and `gavish` - so tables and figures land in two folders and carry the slug in
# their filenames, exactly as `scie` and `emt` do. Nothing produced with all forty has to be
# deleted or re-run to use the subset, and a table cannot be read as the other set's.

_GAVISH_TNBC_MPS = (
    # the twenty-one states a TNBC malignant epithelial cell can be in
    1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 39,
    # and the one that is there to bound them: the in-vitro stress programme, which is what
    # separates a stress dimension from a dissociation one
    7,
)


def _mp_number(sig: Signature) -> int:
    """The MP number of a metaprogram signature, read off the filename stem.

    The stem is the one identifier that is stable: `name` carries a corrected spelling for
    MP3 and the labels are prose. `MP3_CELL_CYLCE_HMG_RICH.txt` -> 3.
    """
    return int(sig.file.rsplit("/", 1)[-1].split("_")[0][2:])


_GAVISH_BY_NUMBER = {_mp_number(s): s for s in _GAVISH_SIGNATURES}

# A number in the list above with no metaprogram behind it is a typo, and it would otherwise
# surface as a quietly shorter collection rather than as an error.
_unknown = sorted(set(_GAVISH_TNBC_MPS) - set(_GAVISH_BY_NUMBER))
if _unknown:
    raise ValueError(f"_GAVISH_TNBC_MPS names metaprograms that are not in the registry: "
                     f"{', '.join(f'MP{n}' for n in _unknown)}")

# Registry order, not the order of `_GAVISH_TNBC_MPS`: the tuple above is grouped by intent so
# that it can be read, and the collection is ordered by family like every other one.
_GAVISH_TNBC_SIGNATURES = tuple(s for s in _GAVISH_SIGNATURES
                                if _mp_number(s) in set(_GAVISH_TNBC_MPS))


def _gavish_flags(coll: "Collection", claimed: str, a_rho: float) -> list[str]:
    """The two things that would make a metaprogram match mean something other than it says.

    Neither is a verdict. Both are printed next to the claim so that the reader does not have
    to remember which of the lists is which.
    """
    axis = coll.axis_of.get(claimed)
    flags = []
    if axis == "emt":
        # The same risk the `emt` collection carries, for the same reason: a mesenchymal
        # programme read on an epithelial compartment is ambient fibroblast RNA or a doublet
        # until the doublet check says otherwise.
        flags.append("mesenchymal_may_be_ambient_or_doublet")
    if axis in GAVISH_CONTROL_AXES:
        # A dimension whose best match is a lineage that cannot be here is not a discovery.
        # It bounds what any OTHER match in the same table is worth.
        flags.append("lineage_control_claim_bounds_the_rest")
    return flags


def _gavish_collection(name: str, scope: str, signatures: tuple[Signature, ...]) -> Collection:
    """One of the two metaprogram collections. They differ ONLY in which lists they carry.

    Written as a factory rather than as two literals so that they cannot drift: the question,
    the flags, the family order and the absence of a target region are the same object twice.
    """
    return Collection(
        name=name,
        title=f"pan-cancer malignant metaprograms (Gavish, {scope})",
        question=("which recurrent pan-cancer metaprogram, if any, does each latent dimension "
                  "carry?"),
        signatures=signatures,
        # Only the families this variant actually has, in the shared order: an axis with no
        # member would otherwise show up as an empty block edge in every heatmap.
        axes=tuple(a for a in _GAVISH_AXES if any(s.axis == a for s in signatures)),
        # No planes, no target label, no risks, no criteria: this collection names dimensions,
        # it does not call cells. `has_target` is False and the steps say so rather than
        # inventing a region nobody asked for.
        extra_flags=_gavish_flags,
    )


GAVISH = _gavish_collection("gavish", "all 41", _GAVISH_SIGNATURES)
GAVISH_TNBC = _gavish_collection("gavish_tnbc", "TNBC-relevant", _GAVISH_TNBC_SIGNATURES)

# --------------------------------------------------------------------------- #

COLLECTIONS = {c.name: c for c in (SCIE, EMT, GAVISH)}
DEFAULT_COLLECTION = "scie"

# `gavish_tnbc` is deliberately NOT a `--collection` value. There is one name for the
# metaprogram vocabulary, `gavish`, and one flag deciding how much of it is scored - so a
# command cannot ask for the subset and the full set at once, and `--all-metaprograms` reads
# as what it is: a widening of the default, not a different body of prior knowledge. The two
# still have two slugs on disk, which is what `resolve` returns and what every path is built
# from; see the block above.
GAVISH_VARIANTS = {False: GAVISH_TNBC, True: GAVISH}


def get(name: str) -> Collection:
    if name not in COLLECTIONS:
        raise KeyError(f"unknown collection {name!r}; have: {', '.join(COLLECTIONS)}")
    return COLLECTIONS[name]


def resolve(args) -> Collection:
    """The collection a step actually runs on: `--collection`, narrowed by `--all-metaprograms`.

    Every step goes through this rather than through `get(args.collection)`, because the slug
    the outputs are named for is decided here and nowhere else.
    """
    coll = get(args.collection)
    if coll is GAVISH:
        return GAVISH_VARIANTS[bool(getattr(args, "all_metaprograms", False))]
    return coll


def add_argument(parser) -> None:
    """The `--collection` flag, spelled identically by all five step scripts."""
    parser.add_argument("--collection", choices=sorted(COLLECTIONS), default=DEFAULT_COLLECTION,
                        help=f"which signature collection to interpret (default {DEFAULT_COLLECTION})")
    parser.add_argument("--all-metaprograms", action="store_true",
                        help="with --collection gavish, score all 41 metaprograms instead of "
                             f"the {len(_GAVISH_TNBC_SIGNATURES)} TNBC-relevant ones; outputs "
                             "are written under the 'gavish' slug rather than 'gavish_tnbc'. "
                             "No effect on the other collections")
