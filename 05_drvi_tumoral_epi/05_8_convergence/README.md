# 05_8_convergence — Route C

Where the two routes agree. **The main result of this stage.** Counterpart of
[04_7_convergence](../../04_drvi_epithelial/04_7_convergence/).

```bash
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
python3 convergence_tum.py                     # scie, the default
python3 convergence_tum.py --collection emt    # the EMT lists
python3 convergence_tum.py --collection gavish # the metaprograms: no target axis, see below
N_LATENT=64 python3 convergence_tum.py         # join the drvi_tum_64 run's tables instead
python3 convergence_tum.py --rho-min 0.30      # a stricter cell-level bar
```

`CELL_SET`, `N_LATENT` and `HVG_SET` select **which 05_3 run** this reads — `drvi_tum_32` by
default, `N_LATENT=64` for `drvi_tum_64`, `CELL_SET=epi HVG_SET=nomt N_LATENT=64` for
`drvi_epicnv_64_nomt`. See [the phase README](../README.md#which-run-05_4---05_8-read); the
run id is in the name of everything written.

Reads five tables written by 05_6 and 05_7, all from `../tables/<collection>/<run_id>/`, so a run can only
ever join a collection with itself.

## Why agreement is the criterion

Route A and Route B traverse the same mapping in opposite directions, and neither set they map
between is ground truth. What makes the pair worth running is that **their failure modes do not
overlap**:

- a dimension can pass Route A by coincidence among heavily correlated per-cell scores;
- a dimension can pass Route B by gene-set overlap with no cellular counterpart at all;
- it is unlikely to pass **both** for the wrong reason.

So agreement is the criterion for calling a dimension a genuine cell state, and disagreement is
informative rather than a failure:

| `verdict` | | reading |
|---|---|---|
| `factor_only` | **B but not A** | the axis carries the gene programme and Route A did not clear the bar. Which of the two things that means is in the `A_vs_null` column, not in the name: **`below`** the collection's noise floor is the classic case — no coherent group of cells sits on the axis, a candidate patient-specific or technical effect — while **`above`** it is a real association the bar dropped for being weak |
| `cell_only` | **A but not B** | the model separates the cells but does not encode the programme cleanly on a single axis: the state is real, the axis is not its description |
| `convergent` | **A and B**, same family | a name |
| `both_different_family` | **A and B**, different family | both routes strong and pointing at different programmes: a signal, not a name. Only `gavish_tnbc` produces these in quantity — it is the only collection with ten families for the two routes to disagree across |
| `neither` | | |

**The names say which routes answered, and nothing else.** That is a correction, not a style:
the B-only category used to be called `factor_only_candidate_patient_or_technical`, which put
the *reading* in the label while the reading depends on where `ROUTE_A_RHO_MIN` sits. At the
0.20 bar the two coincided; at 0.30 the category also catches rows like `DR 41-`, where both
routes land on `MP19_EPITHELIAL_SENESCENCE` (ρ 0.281, FDR 4e-22) and the old name called it a
technical artefact. `A_vs_null` now carries that distinction as a fact of its own.

All three are reported separately and **nothing is promoted on a single route**.

That rule does more work here than in 04. This phase has no non-constant biological annotation to
fall back on — `cell_type` is `malignant` everywhere — so a dimension that passes only Route A
cannot be sanity-checked against a lineage label. Convergence with a gene programme is what
stands in for it.

## How the two routes are joined

Route A is computed **per dimension**, so its direction is the **sign of the Spearman
correlation**: a signature correlating positively with DR 7 is a statement about `DR 7+`, and
negatively about `DR 7-`. Route B is already per direction. That is what makes them joinable.

One row per dimension **and** direction, for every dimension of the run: nothing is pruned
anywhere in this stage, so the table is 2 × 32 = **64 rows** and a dimension DRVI wrote off can
still be read. This step never decides that itself — it takes the dimension list from the
`dimension_row_order` table 05_6 writes, so under `PRUNE_VANISHED=1` it is 2 × 56 = 112 rows
without a line of this script knowing why. See
[the phase README](../README.md#the-vanished-dimension-control) for what the control found.

A readout is flagged as confounded when its raw score correlates with a technical or cycle
covariate above 0.30 (`DEPTH_FLAG`, `CYCLE_FLAG`). Both are conventions, not derived: chosen so
the flag fires on the couplings the confounder table actually shows and stays quiet on the rest.
A convergent dimension carrying a flag is not a result, it is a lead.

## What comes out

```
../tables/<collection>/<run_id>/convergence_*.csv    one row per dimension-direction: both routes, the flags
../tables/<collection>/<run_id>/target_axes_*.csv    the axes that pass the collection's own criteria
```

Figures in `../figures/05_8_convergence/<collection>/<run_id>/`: the two routes side by side on the same
row order, and the convergence scatter.

`target_axes` is the answer to the question each collection was built to ask — "is there a
malignant state that is stem-like and immune-evasive?" for `scie`, "which cells are in the hybrid
E/M state?" for `emt` — restricted to the axes on which both routes agree.

**`gavish` writes no `target_axes`, and that is the whole point of it.** It states no target: it
is a vocabulary, so it has no state it is looking for, and `Collection.has_target` is False. The
`convergence` table *is* the result there — a dimension-direction where Route A and Route B
independently land on the same metaprogram is a dimension with a name that came from outside this
dataset. Two columns are absent from that table for the same reason (`A_auroc_target_this_side`,
`A_standardised_mean_difference`): there is no cell set for a dimension to separate.

## Naming: one bar, not a ladder — `dr_naming_tum.py`

A name is a yes or a no. A direction where **Route A ρ ≥ 0.30 and Route B FDR < 0.05 on the same
signature family** carries the programme; one that does not, does not. Whether ρ came out at 0.31
or at 0.84 is a statement about how tightly the axis tracks the score, not about whether the name
holds, so the magnitude is **reported** (`A_rho`, `max_rho`) and never used to grade the result.
This replaced a four-tier null ladder that was being read as a second filter, which it never was.

**Why 0.30.** `A_rho` is a dimension's *best* gene set, a maximum over the collection's K sets,
and that is exactly what [05_9's per-direction null](../05_9_embedding_control/) measures: 0.30
sits at its p97.5 on `scie`, ~p98 on `gavish_tnbc` and above p99 on `emt` — 1–3 directions in 112
reach it by chance on the two collections that carry the result. Not the per-list null of
`sig_collections` (p95 0.325, p99 0.372), which answers the other question.

**The bar is applied in `dr_naming`, to convergence tables written at a lower one, and that is
exact rather than a shortcut.** `convergence_*.csv` stores `A_rho`, `B_significant` and
`same_family`, and `convergence_tum.py` builds its verdict out of those three and nothing else;
raising a threshold can only turn a `convergent` row into a non-convergent one, never the
reverse. `dr_naming_tum.py --rho-min` below the stored bar is refused, because `same_family` was
never computed for the rows that would come back.

On `drvi_tum_64_nomt`: **35 bars, 28 directions, 27 dimensions, 19 programmes** — three more
directions than the ladder gave (`DR 26-` MP13 EMT 2, `DR 19-` LIM STEM, `DR 40-` MP19 epithelial
senescence, all with ρ between 0.30 and 0.325).

**The `_pruned_rho030` sweep confirms it.** That variant re-runs 05_6–05_8 with
`ROUTE_A_RHO_MIN=0.30`, i.e. the bar baked into the convergence tables instead of applied here,
and its consensus tables and all 36 consensus figures come out **byte-identical** to the ones
under `tables_pruned/` and `figures_pruned/`. Both are kept: they are two independent routes to
the same answer, and the sweep is what makes the shortcut checkable rather than asserted.

The sweep is *not* redundant upstream of this step. Exactly four figures move between the two
tags — the three `convergence_scatter_*` (coloured by a verdict the bar defines) and
`05_6/emt/dim_signature_heatmap` (which draws the bar as a dashed line). Everything else in
`figures_pruned_rho030/` is identical to `figures_pruned/`, because nothing else in 05_4–05_7
reads the number.

## The picture of a name — `dr_umap_tum.py`

`dr_naming_tum.py`'s barplot says **that** a direction carries a programme. This says **where**:
one figure per bar segment, two panels on the same integrated UMAP.

```bash
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
N_LATENT=64 HVG_SET=nomt python3 dr_umap_tum.py                       # all 35 pairs
N_LATENT=64 HVG_SET=nomt python3 dr_umap_tum.py --only "DR 50+,DR 4+" # a few
N_LATENT=64 HVG_SET=nomt python3 dr_umap_tum.py --score raw           # not the within-cohort z
```

| panel | what is coloured |
|---|---|
| **left** | every cell by its **within-cohort z** for that programme — the quantity Route A correlates, `z_<programme>` from 05_6's per-cell score table |
| **right** | every cell by its **coordinate on that dimension** |

**One side only, on both panels.** At or below zero is one flat grey; above it, a single-hue ramp
clipped at the 99th percentile of the positive values. Left is purple. Right is **red for a `+`
direction and blue for a `-` one**, so the colour names the side of the axis — for a `-`
direction the coordinate is negated before the clip and the colourbar counts down from 0, so the
blue cells are the ones with a genuinely negative coordinate.

A diverging map was tried first and was the wrong instrument: it spends half its range on the
cells that do not have the programme, which are exactly the cells neither route is making a claim
about. Grey is the honest colour for *not in it*, and it leaves the eye nothing to do but compare
the two coloured sets. The UMAP is 05_3's, computed on the latent itself (`use_rep="X"`), not on
the HVG space.

**Which pairs are drawn is not decided here.** The list is read out of
`programme_dimensions_consensus_<run>.csv`, the table the barplot is built from, so the figures
are its bars by construction. Same default as `dr_naming_tum.py`: `tables_pruned/`, with
`--keep-vanished` moving both. Nothing is measured; the ρ and the FDR in each subtitle are read
from `convergence_*` and `dim_signature_spearman_*`, and the one number computed on the spot is
the Spearman of the two vectors actually drawn.

**What the run prints that the barplot cannot show.** A bar's `max_rho` is Route A's best
correlation *over the whole collection*; the programme the bar is *labelled* with is Route B's
best enrichment, and the bar asks for the same **family**, not the same signature. On this run
**10 of the 35 bars** are labelled with a signature that is not Route A's best there, and on
three the labelled programme's own ρ is **below the bar** — `MP13_EMT_2` on `DR 26-` (+0.272, the
bar is cleared by `MP12_EMT_1` at +0.323), `ESC_WONG` on `DR 3+` (+0.274, `FMASC` +0.371),
`LIM_STEM` on `DR 32-` (+0.289, `FMASC` +0.334). Those three directions are over the bar on a
sibling of the same family, not on the list they are named after.
