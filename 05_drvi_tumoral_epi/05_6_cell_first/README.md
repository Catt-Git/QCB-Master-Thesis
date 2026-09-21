# 05_6_cell_first — Route A

Prior knowledge defines the state, DRVI is the coordinate system. Counterpart of
[04_5_cell_first](../../04_drvi_epithelial/04_5_cell_first/), on the malignant subset.

**The unit of analysis is the cell.** The question is whether the cells prior knowledge calls
stem-like, immune-evasive or hybrid-EMT occupy a distinct position along any latent dimension.
This is the only route that assigns cells to states — Route B can tell you a dimension is
enriched for stemness genes and still not tell you which cells, how many, or in which patients —
so the deliverable of the project comes from here.

```bash
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
python3 cell_first_tum.py                            # scie, the default
python3 cell_first_tum.py --collection emt           # the EMT lists
python3 cell_first_tum.py --collection gavish        # the metaprograms: scores only, no target
python3 cell_first_tum.py --collection gavish --all-metaprograms   # all 41, not just the TNBC 22
python3 cell_first_tum.py --high-q 0.80 --low-q 0.20 # a stricter target region
python3 cell_first_tum.py --overwrite                # re-score instead of reusing the csv

# the two runs this phase reports, each on all three collections:
N_LATENT=64 python3 cell_first_tum.py [--collection emt|gavish]                # drvi_tum_64_nomt
CELL_SET=epi N_LATENT=64 python3 cell_first_tum.py [--collection emt|gavish]   # drvi_epicnv_64_nomt
```

`CELL_SET`, `N_LATENT` and `HVG_SET` select **which 05_3 run** this reads — `drvi_tum_32_nomt`
by default, `N_LATENT=64` for `drvi_tum_64_nomt`, `CELL_SET=epi N_LATENT=64` for
`drvi_epicnv_64_nomt`. `HVG_SET` no longer has to be spelled out: `nomt` is the default panel
(`HVG_SET=withmt` is the explicit opt-out). See [the phase README](../README.md#which-run-05_4---05_8-read); the
run id is in the name of everything written.

Needs `05_4` (the `.gmt`) and `05_3` (the embedding). `05_5` is optional and its absence is
reported rather than fatal.

Scoring is on the **unintegrated, all-genes** object (`shiao_tum.h5ad`, 42,096 × 25,133): a
150-gene signature reduced to whatever survived HVG selection is no longer that signature.

## The one deliberate departure from 04: the standardisation strata

```python
GROUPBY = ["cohort"]        # 04_5 uses ["cohort", "cell_type"]
```

The reason is **not** that `cell_type` is constant here. That is true, and it would make the
second key a no-op rather than a mistake. The reason is what the second key would be if it were
not constant.

In 04, `cell_type` was a **lineage**: luminal against basal are different cells, and
standardising within them asks "inside a cell type, does this dimension track stemness?" — a
sensible question when the groups are different kinds of cell. Inside one malignant compartment
the available groupings (`cell_type_01_4`, the leiden partition) are **states of the same
tumour**, and state is the quantity this phase measures. Standardising within a state removes
the contrast being looked for, by construction.

The concrete danger, spelled out because it is not hypothetical: **`Lumsec-prol` means
proliferating and is 53% of this subset**, and the cell cycle is a named risk of the `scie`
collection. Standardising within that label would subtract the proliferation axis before
anything is measured — which would not remove the confounder, it would *hide* it, and measuring
it is exactly what the G1 recomputation below is for.

The pre-CNV label and the leiden partition are therefore reported as **covariates** of the
target set (`quadrant_per_*`), never regressed out of it.

## How it fails, and what the script does about it

| failure | the check |
|---|---|
| the stemness lists are embryonic and proliferation-heavy, and half this subset is `Lumsec-prol`, so "stem-high" can just mean "cycling" | the target region is **recomputed inside G1 alone**, cutoffs re-derived there, and the two cell sets compared by Jaccard. Runs for every collection |
| immune evasion is defined by **absence** of signal, which shallow sequencing mimics perfectly | depth of the immunogenic-low group against the rest: Mann-Whitney + AUROC of "shallower predicts evasive" |
| a high mesenchymal score can be fibroblast ambient RNA or a doublet | `doublet_score` (Scrublet) of the mesenchymal-high group, AUROC and Spearman |
| absolute scores are not comparable across patients | standardisation within `cohort` |
| it is confirmatory by construction | no fix inside Route A — it is why Route B exists |

**What the malignant subset changes for the ambient risk.** A cell had to be called aneuploid to
be in this object, so an actual fibroblast is not among them. Ambient RNA is contamination and
not identity, so it survives — the check stays, but half the risk is gone by construction, which
is the one thing 04 could not say.

## The target region

No single definition is primary. The region is defined once per plane — one per stemness readout
for `scie`, one per list version for `emt` — and the **stability of the resulting cell set across
those definitions is itself a reported result** (`quadrant_stability`, Jaccard). A consensus cell
is one called by a majority of the definitions, floored at 2 so no single definition can carry it
alone.

The EMT target is **co-expression** (epithelial-high AND mesenchymal-high), not a band in the
middle of an E–M axis; the argument, and the failed first version, are in
`utils/sig_collections.py`. Its stability numbers there were measured on **04's** cells and have
to be re-read from this phase's own table before being quoted of the tumour.

## What comes out

Embedding-independent (they describe the cells, not the space):

```
../tables/<collection>/<run_id>/confounders_*.csv              readouts vs depth, mito, S, G2M
../tables/<collection>/<run_id>/quadrant_stability_*.csv       Jaccard across definitions
../tables/<collection>/<run_id>/quadrant_vote_distribution_*.csv
../tables/<collection>/<run_id>/quadrant_per_patient_*.csv     a state in one patient is a patient effect
../tables/<collection>/<run_id>/quadrant_per_<grouping>_*.csv  leiden, cell_type_01_4, phase
../tables/<collection>/<run_id>/confounder_checks_*.csv        every named risk, as numbers
$DATA_DIR/05_tum/signature_scores_<collection>_<run_id>.csv   per cell, raw + z
```

Embedding-dependent:

```
../tables/<collection>/<run_id>/dim_signature_spearman_*.csv   dimensions x readouts (the input of 05_8)
../tables/<collection>/<run_id>/dim_target_effect_size_*.csv   AUROC and SMD of the target per dimension
../tables/<collection>/<run_id>/dimension_row_order_*.csv      the row order every later heatmap uses
```

**On a collection with no target region** (`gavish`, where `Collection.has_target` is False)
only `confounders`, `signature_scores`, `dim_signature_spearman` and `dimension_row_order` are
written, plus the confounder and dimensions × readouts heatmaps. Everything with `quadrant` or
`target` in its name is skipped rather than emptied, and the step prints the list of what it
skipped and what it produced before it starts — an all-False consensus written to disk would
read like a result of zero cells instead of like a question nobody asked.

Figures in `../figures/05_6_cell_first/<collection>/<run_id>/`.

One row per **dimension**, not per direction: Route A correlates the latent coordinate, which has
no direction of its own. The direction is the **sign of rho** — that is what makes the two routes
joinable in 05_8, and it is written on the colorbar so the heatmap cannot be misread.

### The derived axis has its own figure

`dim_signature_heatmap` draws the signature **programmes** only — for `emt` that is nine columns
in three blocks, not twelve in four. The derived readouts
(`EMT_SCORE_v = z(EMT_v_MESENCHYMAL) − z(EMT_v_EPITHELIAL)`) were a fourth block there and are
now in `emt_score_axis_*.png`, because a contrast of two columns standing beside those two
columns reads as a fourth *programme* when it is a *coordinate along two of them*. Collections
with no derived readouts (`scie`, `gavish`) are unaffected and get no extra figure.

The derived figure has three panels: **A** every dimension against the coordinate, one column
per list version, so a real axis is a row that agrees across all three; **B** the 14 strongest
dimensions ranked on |rho| — `sign=0` is 05_8's criterion for this axis, an axis being an E-to-M
axis whichever way DRVI oriented it — with one marker per version, so a short bar means the three
curations agree; **C** where the target cells sit on the coordinate, per cell.

**Panel C is the check the co-expression definition is owed.** The target is *co-expression* —
both programmes high at once — and such a cell is near **zero** on a contrast of the two by
construction. Its median lands at percentile **48 / 49 / 50** (B / A / C) of all malignant cells,
which was never imposed: it is the middle band *recovered*, and it is why `EMT_SCORE` is a
readout of the axis and not the definition of the state.

**Colour.** Red–blue (`vlag`) means exactly one thing across this phase — the sign of a
correlation, i.e. which **dimension-direction** is being talked about. The poles of the E-to-M
axis are a property of the *cell* and have nothing to do with the sign DRVI happened to give a
dimension, so reusing red and blue for them would put two unrelated meanings on one pair of hues
in figures that sit side by side. The derived figure uses **purple (epithelial pole) and green
(mesenchymal pole)** instead, in all three panels, and contains no red or blue on purpose.

Vanished dimensions are correlated like every other (`PRUNE_VANISHED = False`): dropping
dimensions first would decide, ahead of the analysis, which axes are allowed to mean something.
