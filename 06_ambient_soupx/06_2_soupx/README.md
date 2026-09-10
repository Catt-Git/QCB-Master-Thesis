# 06_2 - SoupX, channel by channel

Three stages, chained by `soupx_all.sh`, two conda environments:

| stage | script | env | what |
|---|---|---|---|
| prepare | `prepare_soupx_input.py` | `benchmark-py-r` | `shiao.h5ad` -> per-channel matrix, clusters, soup |
| soupx | `run_soupx.R` | `soupx-r` | rho per channel, and the counts it removes |
| assemble | `assemble_soupx.py` | `benchmark-py-r` | the removals put back into one object |

```bash
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
./soupx_all.sh --clean          # everything, resuming, reclaiming the handoff at the end
./soupx_all.sh --dry-run        # print what would run
./soupx_all.sh --samples P01_A_P P02_B_P
./soupx_all.sh --plot           # autoEstCont's diagnostic, one png per channel
```

## Per channel, because the soup is

SoupX works on one 10x run at a time: the ambient RNA in `P01_A_P` came from the cells lysed
in that tube. `shiao.h5ad` holds 148 of them, so stage 1 takes it apart along `obs['sample']`
and stage 2 is a loop. SoupX is fast - seconds to a minute per channel, against inferCNV's
tens of minutes per cohort - so the loop is sequential and the whole of 06_2 is a coffee
break.

## Three decisions worth knowing

**1. The clusters `autoEstCont` conditions on are `optscib_unintegrated_leiden`, not
`cell_type`.** The estimator needs a grouping at roughly cell-type resolution to ask "which
cells express a gene they should not". `cell_type` would be the obvious choice and is the
wrong one: the CellTypist labels are among the things ambient RNA is suspected of having
distorted, and `06_4` re-runs them on the corrected counts. Estimating the correction from
the labels it is meant to fix closes a loop. The 32-cluster leiden partition of `01_5` is used
the way `05_2/clustering_tum.py` uses `cell_type_01_4` for its NMI target - one grouping to
condition on, never a biological claim - and `--cluster-key` overrides it for anyone checking
the estimate does not depend on that choice.

**2. `run_soupx.R` writes the removals, not the adjusted matrix.** The counts layer has
915,882,723 non-zero entries; 148 adjusted matrices in Matrix Market would be ~16 GB of text
on a laptop with ~48 GB free. SoupX only ever *subtracts*, so `toc - adjusted` is
non-negative, sparse and lossless - the script asserts the non-negativity rather than
assuming it - and `assemble_soupx.py` applies it back. The handoff in the other direction is
gzipped for the same reason (~4 GB instead of ~16), and `--clean` deletes it once the
assembly exists.

**3. rho is estimated with a threshold walk, and a failure is recorded rather than hidden.**
`autoEstCont`'s default `tfidfMin = 1.0` can leave a small or homogeneous channel with no
usable marker gene, and the function then errors. Rather than losing the channel, the
threshold is walked down through `1.0, 0.8, 0.5, 0.3` and the value that worked goes into
`rho.csv`. If none work, the channel is corrected at SoupX's own prior of 0.05 and marked
`method = 'fallback_prior'` - corrected at a number that was **not measured**, which is why
`06_3` draws those channels apart and has a switch to drop them.

`forceAccept = TRUE` is passed on purpose: without it an estimate landing on the edge of
SoupX's `contaminationRange` (0.01-0.8) is refused. A channel that really is 80% soup is a
result this phase exists to find, and the place to look at it is the notebook.

## What comes out

```
$DATA_DIR/06_amb/input/genes.tsv                    shared row space, all 30,869 genes
$DATA_DIR/06_amb/input/<sample>/...                 the handoff (removed by --clean)
$DATA_DIR/06_amb/adjusted/<sample>/removed.mtx.gz   what SoupX took out
$DATA_DIR/06_amb/adjusted/<sample>/rho.csv          the estimate and how it went
$DATA_DIR/06_amb/soupx_rho.csv                      those merged, one row per channel
$DATA_DIR/06_amb/soupx_cells.csv                    one row per cell - what 06_3 reads
$DATA_DIR/06_amb/shiao_soupx_all_cells.h5ad         619,693 cells, nothing dropped yet
```

Nothing is dropped in this step. Which cells are too far gone is a threshold read off a
distribution, so it belongs in `06_3`'s notebook - the same split `05_1` makes between
`run_infercnv.R` and `call_malignant.ipynb`.

See `assemble_soupx.py`'s docstring for why the assembled object holds corrected counts in
`.X` instead of a normalisation, and which phase-01 columns survive as landmarks.
