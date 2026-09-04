# 05_8_convergence — Route C

Where the two routes agree. **The main result of this stage.** Counterpart of
[04_7_convergence](../../04_drvi_epithelial/04_7_convergence/).

```bash
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
python3 convergence_tum.py                     # scie, the default
python3 convergence_tum.py --collection emt    # the EMT lists
python3 convergence_tum.py --rho-min 0.30      # a stricter cell-level bar
```

Reads five tables written by 05_6 and 05_7, all from `../tables/<collection>/`, so a run can only
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
still be read.

A readout is flagged as confounded when its raw score correlates with a technical or cycle
covariate above 0.30 (`DEPTH_FLAG`, `CYCLE_FLAG`). Both are conventions, not derived: chosen so
the flag fires on the couplings the confounder table actually shows and stays quiet on the rest.
A convergent dimension carrying a flag is not a result, it is a lead.

## What comes out

```
../tables/<collection>/convergence_*.csv    one row per dimension-direction: both routes, the flags
../tables/<collection>/target_axes_*.csv    the axes that pass the collection's own criteria
```

Figures in `../figures/05_8_convergence/<collection>/`: the two routes side by side on the same
row order, and the convergence scatter.

`target_axes` is the answer to the question each collection was built to ask — "is there a
malignant state that is stem-like and immune-evasive?" for `scie`, "which cells are in the hybrid
E/M state?" for `emt` — restricted to the axes on which both routes agree.
