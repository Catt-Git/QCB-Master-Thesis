# 06_4 - the decision applied, and the labels recomputed

Last step of the phase.

```bash
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
python3 recelltypist_soupx.py
```

```
06_amb/shiao_soupx_all_cells.h5ad  +  06_amb/ambient_keep.csv
                    |
                    v
06_amb/shiao.h5ad          <- the drop-in replacement for datasets/shiao.h5ad
06_amb/cell_annotation_soupx.csv   old label, new label, kept/dropped, one row per cell
```

## Why the labels are recomputed and not inherited

`cell_type` in `shiao.h5ad` was written by `01_4/celltypist_annotation.py` on the
**uncorrected** counts, and ambient RNA is exactly the kind of signal that moves a CellTypist
call. Phase 05 then spends those labels on its two most consequential decisions -
`05_1/prepare_infercnv_input.py` uses them to choose which cells are tested for aneuploidy
*and* to build inferCNV's diploid baseline - so a label distorted by the soup does not just
mislabel a cell, it moves the CNV call of every cell in that cohort. Correcting the counts
and keeping the old labels would have fixed the smaller of the two problems.

It is also nearly free. `01_4` builds its own temporary CP10K+log1p matrix from raw counts,
because that is CellTypist's training normalisation and neither phase's `.X` matches it, so
the re-run needs the corrected counts and nothing else - **not** the 470 GB scran job of
`01_pre_processing/submit_preprocessing_all.slurm`.

This is `01_4`'s procedure verbatim: the same `Cells_Adult_Breast.pkl`, the same
`majority_voting=True`. The phase-01 PCA, neighbour graph and UMAP are dropped from the
temporary copy first, so the over-clustering behind the majority vote is built on *this*
population and *this* matrix - the same reasoning as
`05_1/recelltypist_nonmalignant.py`.

## What is kept beside the new labels

| column | what |
|---|---|
| `cell_type` | the new CellTypist label, on corrected counts |
| `cell_type_01_4` | the phase-01 label. Same name, same status phase 05 gives it |
| `celltypist_predicted` | the new per-cell prediction before the majority vote |
| `fraction` | recomputed from the new labels by 01_4's rule - it is a *function* of `cell_type` and would otherwise be stale the moment the labels change |
| `fraction_01_4` | the phase-01 value, so the two can be crosstabbed |

The script refuses to write over `$DATA_DIR/shiao.h5ad`, prints how many cells changed label,
and asserts the whole vocabulary falls in the lineage sets `05_1` expects - a label outside
them would make `prepare_infercnv_input.py` stop two steps later.

## Then

```bash
cd ../utils && ./link_shadow_data_dir.sh
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets/06_amb
```

and phase 05 runs verbatim on the corrected object, into `06_amb/05_tum/`, with no edit to
any file under `05_drvi_tumoral_epi/`. Read the caveat about figures and tables in the phase
README before doing it.
