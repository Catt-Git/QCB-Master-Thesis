# 06_ambient_soupx

Seventh phase of the thesis. It has one job: build the version of `shiao.h5ad` that phase 05
should have started from, with the ambient RNA taken out.

## Why it exists

Reading the DRVI dimensions of the malignant epithelial subset (`05_3`) left the impression
that a large share of the cells that survived into it are carrying ambient RNA rather than
their own transcriptome. That is a specific, checkable claim about droplet data: in a
dissociated solid tumour, cell-free mRNA from lysed cells is distributed into every droplet
of the emulsion, so a droplet with little RNA of its own ends up expressing the *average of
its channel* - which in this dataset means immune transcripts in epithelial cells, and
epithelial transcripts in immune cells.

Phase 05 is more exposed to it than any other phase in this repo, for two reasons that
compound:

1. **The CellTypist labels decide the whole run.** `05_1/prepare_infercnv_input.py` uses
   `cell_type` twice - to choose which cells are tested for aneuploidy (the epithelium) and
   to build inferCNV's diploid baseline (T/NK and myeloid). A label pushed around by ambient
   RNA does not just mislabel one cell; it moves the reference every CNV call in that cohort
   is measured against.
2. **The malignant call is a threshold on a continuous score.** A cell whose profile is
   partly its channel's average has a residual profile that is partly the channel's average
   too, which is exactly the failure mode `cnv_corr` was added to catch and cannot fully
   catch on its own.

Phase 04's README already names this as a live risk in a different context: it lists
"fibroblast ambient RNA / doublets" as a way a cell can score mesenchymal without being in a
hybrid state, and 05's README repeats that ambient RNA *survives* the malignant subsetting
because "it is contamination, not identity". This phase is what turns that acknowledged
limitation into a correction.

## What it does not do

**Phase 05 is not modified.** No file under `05_drvi_tumoral_epi/` is edited by this phase,
and nothing here writes into `$DATA_DIR/05_tum/`. Phase 06 is a fourth branch off
`shiao.h5ad` - after 03, 04 and 05 - and owns `$DATA_DIR/06_amb/` alone.

The way 05 is re-run on the corrected data instead is a **shadow `DATA_DIR`**, and it needs
no code change at all. Every script of phase 05 reads `$DATA_DIR/shiao.h5ad` and writes under
`$DATA_DIR/05_tum/`; `06_4` therefore names its output `shiao.h5ad` and puts it in
`06_amb/`, and `utils/link_shadow_data_dir.sh` symlinks the three auxiliary inputs 05 also
expects (`Cells_Adult_Breast.pkl`, `regev_lab_cell_cycle_genes.txt`, `signatures/`) beside
it. Pointing `DATA_DIR` at `06_amb/` then runs phase 05, verbatim, on the corrected object,
into `06_amb/05_tum/`.

The one thing that mechanism cannot redirect: 05 writes its **figures and tables into the
repo**, at `05_drvi_tumoral_epi/figures/` and `.../tables/`, and those paths come from where
the scripts live rather than from `DATA_DIR`. A re-run overwrites them. Copy them aside first
if the uncorrected versions matter - `utils/link_shadow_data_dir.sh` says how, and
`git status` shows anything that moved either way.

## Why SoupX

| | needs | verdict |
|---|---|---|
| **SoupX** | the empty-droplet profile of each channel | **chosen**: the background is *measured*, and only a ~300 KB summary per channel has to come off the cluster |
| CellBender | the full unfiltered matrices + a GPU | the cluster has no GPU (see `05_3/submit_drvi_tum.slurm`), the local card is a 4 GB T1000, and 148 channels would be days; it also redoes cell calling, which would fight with `01_2`'s scrublet and QC decisions |
| decontX | nothing beyond `shiao.h5ad` | cheapest, but it *infers* the background from the cells themselves instead of measuring it in the droplets that have no cells |

SoupX's cost is that it needs something `datasets/` does not contain. `shiao.h5ad` was built
from `filtered_feature_bc_matrix.h5` (`00_6_build_h5ad/build_combined_h5ad.py`), so the empty
droplets - the only place the ambient profile can be read - were dropped at the very first
step. They are still in the Cell Ranger output on the cluster, next to the FASTQ. Hence the
shape of this phase: **one cluster step that reduces millions of droplets per channel to one
table of 30k numbers, and everything else local.**

## The chain

| step | where | what |
|---|---|---|
| `06_1_soup_profile/` | **cluster** | `raw_feature_bc_matrix.h5` -> `06_amb/soup/<sample>.csv.gz`, the ambient profile of each channel |
| `06_2_soupx/` | local | `shiao.h5ad` + those profiles -> SoupX per channel -> `06_amb/shiao_soupx_all_cells.h5ad` |
| `06_3_ambient_qc/` | local, **notebook** | how much soup, and which cells are mostly it -> `06_amb/ambient_keep.csv` |
| `06_4_recelltypist/` | local | applies that call, re-runs CellTypist on the corrected counts -> `06_amb/shiao.h5ad` |

Same division of labour as `05_1`: the scripts produce a continuous quantity, the notebook
chooses where to cut it. `06_2/soupx_all.sh` chains the three stages of 06_2 the way
`05_1/infercnv_all.sh` chains its two.

### Why CellTypist is re-run (`06_4`)

Because the labels are the thing most likely to have been wrong, and phase 05 spends them on
its two most consequential decisions (above). Re-annotating is also nearly free here:
`01_4/celltypist_annotation.py` builds its own CP10K+log1p matrix from raw counts, so the
re-run needs the corrected counts and nothing else - **not** the 470 GB scran job of
`01_pre_processing/submit_preprocessing_all.slurm`. The phase-01 label is kept beside the new
one as `cell_type_01_4`, the same name and the same status phase 05 already gives a
CellTypist label that is no longer current.

### What the corrected object is, exactly

`06_amb/shiao.h5ad` is a drop-in for `datasets/shiao.h5ad` with three differences, all of
them deliberate and all of them recorded in `uns['ambient_soupx']`:

- **`layers['counts']` are the corrected counts**, integer (`roundToInt = TRUE`, because 05_2
  asserts integer counts and scran/DRVI/inferCNV all model counts).
- **`.X` is the same matrix as `layers['counts']`, not a normalisation.** Correcting the
  counts invalidates the phase-01 scran size factors that were estimated from them, and
  re-running scran on 619,693 cells is that 470 GB job. Carrying a stale normalisation that
  still looked valid would be the worse option, so `size_factors` is dropped and `.X` holds
  counts. This costs nothing downstream: `05_2/subset_and_qc.ipynb` copies
  `layers['counts']` into `.X` and drops the size factors as its first act, then re-runs
  scran on the subset it builds.
- **Fewer cells**, by however much `06_3` decided to drop.

The QC columns (`total_counts`, `n_genes_by_counts`, `pct_counts_mt`, ...) are recomputed on
the corrected counts, with the phase-01 values kept under `<name>_precorrection`. The
phase-01 PCA, UMAP, neighbour graph, leiden sweep, cell-cycle scores and scrublet columns are
kept unchanged as **landmarks**: they describe the uncorrected data, they are named in the
object for what they are, and `05_2` drops all of them anyway.

## Running it

```bash
mamba env create -f ../environments/soupx-r.yml       # once, for 06_2's R half

# --- on the cluster (once; ~hours for 148 channels, then a few tens of MB come back) ---
cd 06_ambient_soupx/06_1_soup_profile && mkdir -p logs
export DATA_DIR=/users/genomics/albertoc/Tesi/hopes_and_dreams/datasets
sbatch --export=ALL,DATA_DIR=$DATA_DIR submit_soup_profile.slurm
rsync -av <cluster>:$DATA_DIR/06_amb/soup/           ~/Desktop/QCB-Master-Thesis/datasets/06_amb/soup/
rsync -av <cluster>:$DATA_DIR/06_amb/soup_census.csv ~/Desktop/QCB-Master-Thesis/datasets/06_amb/

# --- local ---
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
cd ../06_2_soupx && ./soupx_all.sh --clean      # prepare + 148 channels + assemble
                                                # then 06_3_ambient_qc/ambient_qc.ipynb
cd ../06_4_recelltypist && python3 recelltypist_soupx.py

# --- phase 05 on the corrected object, without touching phase 05 ---
cd ../utils && ./link_shadow_data_dir.sh
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets/06_amb
cd ../../05_drvi_tumoral_epi/05_1_infercnv && ./infercnv_all.sh --threads 12
# ...and the rest of 05 exactly as its own README describes it
```

### Disk

The counts layer of `shiao.h5ad` has **915,882,723 non-zero entries**, and `/home` has ~48 GB
free. Two consequences shape 06_2 and are worth knowing before starting:

- the Python -> R handoff is **gzipped** Matrix Market (~4 GB for all 148 channels instead of
  ~16 GB), and `soupx_all.sh --clean` removes `06_amb/input/` once the assembly succeeds;
- `run_soupx.R` writes back only the **removals** (`toc - adjusted`), not the adjusted matrix.
  SoupX only ever subtracts, so the delta is non-negative and lossless - the script asserts
  it rather than assuming it - and it is a fraction of the size.

Budget roughly: ~4 GB for the handoff (reclaimable), ~4 GB for
`shiao_soupx_all_cells.h5ad` (reclaimable once `06_4` has run), ~4 GB for
`06_amb/shiao.h5ad`, and ~8.5 GB more if phase 05 is then re-run into `06_amb/05_tum/`.

## Limits to state

- **rho is per channel, not per cell.** SoupX estimates one contamination fraction per
  emulsion, which is what the method is; the per-cell quantity `06_3` thresholds on
  (`soupx_frac_removed`) is how much of that channel's soup the correction actually took out
  of a given droplet, not an independent per-cell estimate.
- **`autoEstCont` can fail on a channel**, and on a dataset of 148 mostly small channels it
  sometimes does. `run_soupx.R` walks `tfidfMin` down before giving up, records which value
  worked, and falls back to SoupX's own prior of 0.05 with `method = 'fallback_prior'` for
  the channels where nothing worked. Those channels are *corrected at a number that was not
  measured*; they are flagged in `soupx_rho.csv`, drawn separately in the notebook, and
  `06_3` has a switch to drop them outright.
- **The empty-droplet definition is inherited, not derived.** `06_1` sums droplets with
  1-100 UMIs because that is SoupX's own `soupRange`; it is an assumption about where cells
  stop, and it is the one parameter of this phase most worth varying if a channel looks odd.
- **The correction is not a doublet filter and not a QC filter.** Phase 01 did those. A cell
  dropped in `06_3` is one whose remaining signal is too thin to carry an identity, which is
  a different statement.
