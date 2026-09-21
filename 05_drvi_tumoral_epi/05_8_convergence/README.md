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

| | reading |
|---|---|
| **B but not A** | the axis carries the gene programme but no coherent group of cells sits on it — a candidate patient-specific or technical effect |
| **A but not B** | the model separates the cells but does not encode the programme cleanly on a single axis: the state is real, the axis is not its description |
| **A and B** | convergent |

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
