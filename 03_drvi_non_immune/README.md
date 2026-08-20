# 03_drvi_non_immune

Fourth phase of the thesis and the start of **Part 2 (biological interpretation)**. Phase 02
answered a technical question - which batch-correction method best preserves biology on this
dataset - on all 619,693 cells at once. This phase changes question: it isolates the **non-immune
compartment** and asks what DRVI's latent dimensions say about it.

The compartment is defined by `fraction == 'non_imm'`, the biological immune/non-immune class
recoded in 01_4 from the CellTypist annotation (**not** the technical CD45 sort, which stays in
`dataset_origin` and is not used here).

Input from phase 01: `shiao.h5ad` (619,693 cells x 30,869 genes).
Integration `batch_key = 'cohort'` (34 patients), biological `label_key = 'cell_type'`.

## Repository layout

```
03_drvi_non_immune/
├── README.md
├── 03_1_subsetting/               # from shiao.h5ad to the non-immune object
├── 03_2_drvi_run/                 # DRVI on that object, n_latent = 64 (notebook + headless run)
│                                 #   + the vanished-dimension pruning check
├── 03_3_enrichment/               # what the latent dimensions mean: GSEApy / Enrichr
└── figures/                       # one folder per step: 03_1_*, 03_2_<run_id>, 03_3_*
```

## 03_1_subsetting

A full re-run of the phase-01 pre-processing on the subset. Nothing derived is inherited: scran
size factors, HVGs, PCA, neighbours, UMAP and leiden were all computed on a population that was
71% immune, and none of them describe these cells. Only the raw counts and the metadata carry over.

The split between notebooks and scripts follows 01: notebooks where decisions are taken by looking
at figures, scripts where the computation is deterministic and long.

### Execution order

| # | File | What it does | Where |
|---|------|--------------|-------|
| 1 | `subset_and_qc.ipynb` | Subset `fraction == 'non_imm'`; restore `.X` = raw counts; drop every phase-01 derived slot (PCA/UMAP/neighbours/HVG/size factors/leiden/cell cycle) and the emptied categories; recompute QC metrics; re-apply the cell thresholds; gene filter `min_cells=3`; batch census. | local (notebook) |
| 2 | `subsetting_all.sh` | Driver for steps 3-6: runs them in sequence, resumes from the last completed one, logs to `logs/`. | local |
| 3 | `scran_norm_nonimm.py` | scran normalization re-estimated on the subset (scib: `quickCluster → computeSumFactors → logNormCounts`); adds `size_factors`, clears `.raw`, casts `.X` to float32. | local |
| 4 | `cell_cycle_score_nonimm.py` | Tirosh/Regev cell-cycle re-scoring → `S_score`, `G2M_score`, `phase`. | local |
| 5 | `reduce_data_nonimm.py` | Batch-aware HVG (`cohort`, 2000) + PCA(50) + neighbours + UMAP; writes the HVG list **and** the 2,000-gene DRVI input. | local |
| 6 | `clustering_nonimm.py` | Leiden optimal-resolution sweep (0.1-2.0) vs `cell_type` NMI → `shiao_nonimm.h5ad`. | local |
| 7 | `visualization_nonimm.ipynb` | Figures on the final object: composition, cell cycle, normalization, HVG overlap with phase 01, UMAPs, leiden sweep. Read-only. | local (notebook) |

> **Note on step 2.** The scripts are duplicated from 01_3/01_4/01_5 rather than called with
> different paths, so this phase reads as a self-contained Materials & Methods section and the
> phase-01 scripts stay frozen. The parameters are identical except where stated below.

> **Note on where.** 01_3 (scran) needed the cluster; here it does not. 176,610 cells instead of
> 619,693, so the whole chain runs locally in sequence and there is no SLURM wrapper.

### QC figures (`figures/03_1_subset_and_qc/`, from step 1)

Same panels as 01_2 (`n_genes_by_counts`, `total_counts`, `pct_counts_mt`) plus
`pct_counts_ribo`, at **three rounds** - `before_filter`, `after_cell_and_gene_filter`,
`after_cohort_drop` - and in **two versions** each:

- `violin_nonimm_qc_<round>_unintegrated.png` - clean, the distribution alone;
- `violin_nonimm_qc_<round>_unintegrated_jitter.png` - the same with the cells drawn on top.

The jittered version strips a random 20,000-cell subsample (fixed seed), not the whole
compartment: at 176k points the strip becomes a solid black band that hides the violin it is
supposed to annotate. Plus `scatter_nonimm_counts_vs_genes_unintegrated.png`, the complexity
scatter coloured by mitochondrial content.

### The low-complexity tail (step 1, `low_gene_table`)

The cell filter is a single cut at `n_genes_by_counts > 100`, and a single cut says nothing about
how close the surviving cells sit to it. The notebook prints the distribution in ranges just
before applying the filter, so the threshold is defended by the shape of the tail rather than
inherited:

| n_genes_by_counts | cells | % | cum % | median pct_mt |
|---|---:|---:|---:|---:|
| ≤ 100 | **0** | 0.00 | 0.00 | - |
| 101-150 | 548 | 0.31 | 0.31 | 0.22 |
| 151-200 | 394 | 0.22 | 0.53 | 0.40 |
| 201-250 | 653 | 0.37 | 0.90 | 1.57 |
| 251-300 | 906 | 0.51 | 1.42 | 2.83 |
| 301-400 | 8,464 | 4.79 | 6.21 | 3.34 |
| 401-500 | 13,547 | 7.67 | 13.88 | 3.09 |
| 501-750 | 26,727 | 15.13 | 29.01 | 3.32 |
| 751-1,000 | 19,260 | 10.91 | 39.92 | 3.35 |
| 1,001-2,000 | 52,662 | 29.82 | 69.74 | 3.36 |
| 2,001-4,000 | 43,551 | 24.66 | 94.40 | 3.31 |
| > 4,000 | 9,898 | 5.60 | 100.00 | 4.34 |

176,610 cells; min 101, median 1,289, max 9,355 detected genes.

Two things to read out of it:

- **The top row is empty.** Not one non-immune cell carries 100 genes or fewer, so re-applying the
  01_2 cell filter removes nothing on this compartment - the subset inherits a population 01_2 had
  already cleaned. The filter stays in the notebook as a re-applied check, not as a step that does
  work here.
- **The gene-poorest ranges are also the mitochondria-poorest** (0.22% and 0.40% median, against
  ~3.3% for the bulk). A gene-poor cell is usually a dying one and dying cells are mt-rich, so the
  inversion says these ~950 cells are shallow droplets - few reads of everything, mitochondrial
  transcripts included - rather than damaged cells. That is also why `pct_counts_mt < 10` does not
  catch them: there is nothing mitochondrial there to flag. They are 0.53% of the compartment and
  are kept.

### Usage

```bash
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets

# 1. the notebook, interactively (kernel: benchmark-py-r)
#    -> $DATA_DIR/03_nonimm/shiao_nonimm_raw.h5ad

# 2. the headless chain
cd 03_drvi_non_immune/03_1_subsetting
./subsetting_all.sh                 # resumes: a step whose output exists is skipped
./subsetting_all.sh --force         # re-run everything
./subsetting_all.sh --dry-run       # print what would run
./subsetting_all.sh reduce cluster  # only the named steps (norm, cc, reduce, cluster)

# 3. visualization_nonimm.ipynb, interactively
# 4. 03_2_drvi_run/drvi_nonimm.ipynb, interactively (GPU)
#    or the same run headless: 03_2_drvi_run/run_drvi_nonimm.py (local or SLURM)
# 5. 03_3_enrichment/enrichment_nonimm.ipynb, interactively (no GPU, needs the network)
```

`CC_GENES` defaults to `$DATA_DIR/regev_lab_cell_cycle_genes.txt`; export it to override.
The active conda environment must be `benchmark-py-r` (scib + rpy2 + R scran): the driver does
not activate it.

## 03_2_drvi_run

DRVI on the object 03_1 produced. Phase 02 ran DRVI on all 619,693 cells to answer a *technical*
question - how well it integrates `cohort` against the other nine methods. Here the question is
biological: what the latent dimensions say about the non-immune compartment on its own.

A notebook, `drvi_nonimm.ipynb`, for the same reason DRVI is a notebook in 02_2: the latent
size is chosen by eye, from how many dimensions vanish. **`N_LATENT = 64` is the run of this
phase.** The size chosen in 02_2 does not carry over - 176k cells and 18 labels are not the whole
dataset - so it was re-chosen here. `N_LATENT` in the configuration cell is the only thing to edit:
the run id (`drvi_nonimm_<N>`), the model, the embedding and the figure folder all follow from it,
so a re-run at another size never overwrites this one.

Three differences from the 02_2 notebook:

- the input is `shiao_nonimm_hvg_2k.h5ad` (176,610 cells x 2,000 HVGs selected **on this
  compartment**), not the whole-dataset `shiao_hvg_2k.h5ad`;
- nothing is written for the scib benchmark - 02_3 and 02_4 do not run on this compartment - so
  the outputs are the latent space and the full-gene object carrying it, for 03_3;
- `fraction` is constant (`non_imm`) after the subset and is no longer plotted, and
  `dataset_origin` (the CD45 sort) is left out as everywhere in this phase.

**Model.** `DRVI` with `n_latent=64`, `encoder_dims=[256, 128]`, `decoder_dims=[128, 256]` (against
the `[128, 128]` defaults) and `dispersion='gene-batch'`, trained on `layers['counts']` with
`batch_key='cohort'`, `SEED=123`, up to 400 epochs with early stopping after 50 without
improvement (scvi-tools' default 0.9/0.1 split). Faster than the 02_2 runs: same architecture on
176,610 cells instead of 619,693.

**Outputs.** `model_<run_id>.pt` and `embed_<run_id>.h5ad` (the latent space with the dimension
stats and the OOD/IND scores, what 03_3 reads), plus `shiao_nonimm_<run_id>.h5ad` - the definitive
03_1 object, all genes, with the latent space added as `obsm['X_drvi']`. The embedding carries no
genes by construction, so that third file is the compartment itself in the DRVI space, for any
downstream step that needs genes and latent coordinates in the same object.

`OVERWRITE = False` in the configuration cell reuses the model and embedding already on disk, so
re-running the notebook to redraw a figure does not retrain.

### The same run headless (`run_drvi_nonimm.py`, `submit_drvi_nonimm.slurm`)

Choosing the latent size needs a pair of eyes; training at a chosen size does not, and it is the
long half of the phase. `run_drvi_nonimm.py` is the notebook without the kernel - same input, same
architecture, same `SEED`, same early stopping, the same three outputs and the same 25 figures - so
a size can be trained wherever it is convenient and the notebook re-opened afterwards with
`OVERWRITE = False`, which reads model and embedding from disk instead of recomputing them. As in
`subsetting_all.sh`, an output already there is reported `[have]` and reused, so a crash after
training never costs the training.

```bash
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
cd 03_drvi_non_immune/03_2_drvi_run

python run_drvi_nonimm.py                 # n_latent 64, the run of this phase
python run_drvi_nonimm.py --n-latent 32   # another size, side by side with it
python run_drvi_nonimm.py --overwrite     # retrain and rewrite everything

# on the cluster (CPU there, hence partition `long`; 64G is enough for 176k cells)
mkdir -p logs
sbatch --export=ALL,DATA_DIR=$DATA_DIR submit_drvi_nonimm.slurm --n-latent 64
```

`submit_drvi_nonimm.slurm` passes everything after the script path to the script unchanged,
activates `$DRVI_ENV` (default `catalano_env`, which pins the same `drvi-py==0.2.7` as the local
`benchmark-py-r`, so the model it writes loads back here) and logs to `03_2_drvi_run/logs/`. The
inputs to copy up first are `$DATA_DIR/03_nonimm/shiao_nonimm_hvg_2k.h5ad` and - only if the
downstream object is wanted there - `shiao_nonimm.h5ad`, by hand:
`02_integration_benchmark/utils/sync_to_cluster.sh` only walks the phase-02 grid. A missing
`shiao_nonimm.h5ad` is announced at the start and costs nothing else: the run trains and writes
the model and the embedding either way. Only `embed_<run_id>.h5ad` has to come back, since that is
all 03_3 reads.

The number to read in the log before deciding whether to re-run at another size is the vanished
count, printed as it is in the notebook (`12 vanished / 64 latent dimensions`).

### Figures (`figures/03_2_drvi_nonimm_64/`, notebook or `run_drvi_nonimm.py`)

Same set as the DRVI figures of 02_2, each name suffixed with the run id, next to the `03_1_*`
folders:

- `umap_<key>_<run_id>.png` - one UMAP of the DRVI space per metadata/QC key, plus
  `umap_combined_<run_id>.png` and `umap_per_cell_type_<run_id>.png` (one panel per label, since
  the three largest - Fibro-matrix, Lumsec-prol, LummHR-major - dominate the single panel).
  Each has its direct counterpart in `figures/03_1_visualization/umap_*_nonimm.png`: same cells in
  the *unintegrated* PCA space, same keys and same palettes, so only the space changes.
- `latent_dimension_stats[_rmVanished]_<run_id>.png` - per-dimension reconstruction effect, with
  and without the vanished dimensions: the plot behind the latent-size choice.
- `latent_dims_in_umap_<run_id>.png` and `latent_dims_in_heatmap_<key>_<run_id>.png` - each
  non-vanished dimension on the UMAP, and how the dimensions respond to `cell_type` (also sorted
  by label), `cohort`, `treatment`, `response`, `phase`.
- `ood_*_<run_id>.png` / `ind_linear_weighted_mean_<run_id>.png` - interpretability scores. OOD
  comes from the decoder reconstructions (fast, favours the genes *specific* to a dimension,
  `OOD_min/max` being its two halves); IND averages each factor's effect over all cells (broader,
  a gene shared by several dimensions keeps a high score in all of them).

### Pruning the vanished dimensions (`plot_pruned_umap_nonimm.py`)

The DRVI paper defines a latent dimension as **vanished** when its maximum absolute value is
below 1, and Supplemental Note 7 assumes the vanished ones are pruned *before* anything else is
evaluated. `run_drvi_nonimm.py` does not prune: it builds its neighbour graph with `use_rep='X'`
on all 64 dimensions, so every UMAP in the folder above is an unpruned one. This script is the
check on that, the phase-03 counterpart of `02_3_plot_method_umap/plot_drvi_pruned_umap.py`
(kept self-contained rather than imported across phases, as everything else here is).

```bash
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
cd 03_drvi_non_immune/03_2_drvi_run
python plot_pruned_umap_nonimm.py --report-only   # the numbers, in seconds
python plot_pruned_umap_nonimm.py                 # + the figure
python plot_pruned_umap_nonimm.py --no-prune      # the control, see below
```

**The result on `drvi_nonimm_64`: pruning changes nothing measurable**, exactly as on the full
dataset at 128.

- The threshold is not a judgement call. 12 dimensions vanish and 52 are kept, separated by two
  orders of magnitude - vanished at max |z| ≤ 0.034, kept at max |z| ≥ 3.80 - so 0.1, the 0.5
  `set_latent_dimension_stats` was called with, and the paper's 1 select the identical set. The
  stored `var['vanished']` needed no correction.
- The 12 vanished dimensions carry **1.4e-05 of the total latent variance**.
- Pairwise euclidean distances move by at most **2.8e-04 relative** (mean 8.0e-06), and the
  **exact 15-NN overlap is 0.99997** on 20,000 sampled cells.

**The figures cannot settle this; the numbers do.** UMAP's approximate neighbour search and its
spectral initialisation are sensitive enough that the same data laid out again lands in a globally
different arrangement, at the same seed and on the same cells in the same order. `--no-prune`
re-runs the identical code path on all 64 dimensions as a control, and the Procrustes disparity
between 50,000-cell layouts shows how little any of it means:

| pair | disparity, this run (nonimm 64) | disparity, 02 (drvi 128) |
|---|---|---|
| pruned ↔ its own control | 0.648 | 0.176 |
| pruned ↔ the pre-existing UMAP | 0.203 | 0.426 |
| control ↔ the pre-existing UMAP | 0.663 | 0.385 |

The ordering **reverses between the two phases**: here the pruned/control pair is the furthest
apart of the three, in 02 it is the closest. Run-to-run variation swamps whatever the pruning
does, so no pair of these layouts can be read as a before/after - including the control pair. The
control's job is to demonstrate that instability, not to supply a matching picture.

The evidence that pruning is inert is the 15-NN overlap and the distance deviation above, both
computed on the embedding itself. The figures only show that the biology survives: cohorts mixed,
cell types separated, in every version.

The script **writes nothing over anything**: it reads `03_nonimm/embed_drvi_nonimm_64.h5ad`,
caches each layout as `03_nonimm/umap_<run_id>_{rmVanished,allDims_control}.npy` rather than
inside an `.h5ad`, and only adds files to `figures/03_2_<run_id>/`.

## 03_3_enrichment

First downstream analysis of the 03_2 run, and the step that turns the latent dimensions from gene
lists into biology. A single notebook, `enrichment_nonimm.ipynb`: each dimension becomes a ranked
gene list, its top `N_TOP_GENES` genes go to Enrichr through `gseapy`, and the significant terms
are drawn as one barplot per dimension.

**Directions.** A dimension is read separately in its two directions - `DR 12+` is what goes up
when the dimension increases, `DR 12-` what goes up when it decreases - because the two are
different programs. At `N_LATENT = 64` that gives 128 candidate lists, of which **104** survive
(DRVI marks the rest as vanished in that direction and `get_interpretability_scores` drops them).

**Scores.** `OOD_combined`, the out-of-distribution scores of 03_2, which favour the genes
*specific* to a dimension. `IND_linear_weighted_mean` is the alternative (`SCORE_KEY` in the
configuration cell): broader, but it keeps a shared gene high in every dimension it belongs to and
would push the same generic terms into every barplot.

**No model is loaded.** The scores are already in `embed_drvi_nonimm_64.h5ad`
(`.varm['OOD_combined_positive' / '_negative']`), so the notebook rebuilds the genes x directions
table straight from the embedding - the same table `model.get_interpretability_scores(embed, adata)`
returns, checked value by value - and needs neither a GPU nor `scvi-tools`. It reads
`shiao_nonimm_hvg_2k.h5ad` `backed='r'`, for its `var_names` alone.

**Background.** The universe is the **2,000 HVGs**, not the whole genome (`HVG_BACKGROUND = True`,
which makes `gseapy` run the hypergeometric test locally instead of calling the web endpoint). The
genes were only ever rankable among the HVGs, so testing against all human genes would call any
term rich in variable genes - ECM, cycle, interferon - enriched before any dimension is considered.

**Libraries.** `MSigDB_Hallmark_2020` (50 non-redundant programs: EMT, hypoxia, interferon, E2F/G2M
- the first read of a dimension), `Reactome_Pathways_2024` (the level below: collagen formation,
ECM organization, integrin signalling, the resolution a fibroblast-dominated compartment needs),
`GO_Biological_Process_2025` (widest coverage, the only one that reliably annotates the vascular
and perivascular populations) and `KEGG_2021_Human`. GO CC/MF are left out (where a protein sits,
not what a program does), and so are the cell-identity libraries (CellMarker, PanglaoDB, Azimuth):
which cell type a dimension belongs to is the question the `cell_type` heatmaps of 03_2 answer.

`N_TOP_GENES` is the only parameter to edit. 200 is the default - deep enough for a program-level
read, since Hallmark and Reactome terms carry 50-200 genes, and still the top 10% of the 2,000
HVGs. `OVERWRITE = False` re-reads the enrichment table from disk, so redrawing the barplots does
not re-query Enrichr.

### Figures (`figures/03_3_enrichment_<N_TOP_GENES>/`)

One folder per depth, so runs at 50 / 200 / 500 sit side by side, and the run id in every file
name as in 03_2:

- `DR_<nn>_<pos|neg>_enrichr_barplot_<run_id>.png` - the terms below FDR 0.05, up to 6 per library,
  bars coloured by library and scaled by −log10 of the adjusted p-value. Long term names are
  truncated at 60 characters; the full ones are in the results table.

At `N_TOP_GENES = 200`: **102 of the 104** dimension-directions have at least one term at FDR
&lt; 0.05 and get a barplot; `DR 20+` and `DR 41+` have none and are skipped. A direction with no
enriched term is a result, not a failure - it is DRVI having found an axis no annotated gene set
describes.

### Result tables

Both written to `$DATA_DIR/03_nonimm/`, not to the repo, and both carrying `N_TOP_GENES` and the
run id in the name:

- `top<N>_genes_<run_id>.tsv` - the ranked gene lists themselves, one column per direction.
- `enrichr_top<N>_<run_id>.tsv` - every tested term, with the dimension as the first column
  (253,391 rows at 200, 5,261 of them significant). This is the cache `OVERWRITE = False` reads.

## Data location (`DATA_DIR`)

As in every phase, the heavy objects live **outside** the repo and every script resolves `DATA_DIR`
from the environment. This phase writes into its own subfolder, so the non-immune chain is never
confused with the whole-dataset one:

    $DATA_DIR/
    ├── shiao.h5ad                                # input, from 01_5
    └── 03_nonimm/
        ├── shiao_nonimm_raw.h5ad                 # 1  (subset + refilter, raw counts)
        ├── shiao_nonimm_norm.h5ad                # 3  (scran log-norm)
        ├── shiao_nonimm_norm_cc.h5ad             # 4  (+ cell cycle)
        ├── shiao_nonimm_reduced.h5ad             # 5  (+ HVG/PCA/neighbours/UMAP)
        ├── shiao_nonimm_hvg_2k_list.csv          # 5  (selected HVG symbols)
        ├── shiao_nonimm_hvg_2k.h5ad              # 5  (DRVI input for 03_2)
        ├── shiao_nonimm.h5ad                     # 6  (+ leiden) - definitive non-immune object
        ├── shiao_nonimm_leiden_resolution_profile.csv   # 6  (resolution vs NMI)
        │
        ├── model_drvi_nonimm_64.pt               # 03_2  the trained model, one flat file per run
        ├── embed_drvi_nonimm_64.h5ad             # 03_2  latent space + interpretability scores
        ├── shiao_nonimm_drvi_nonimm_64.h5ad      # 03_2  the 03_1 object + obsm['X_drvi'] (03_3
        │                                         #       itself works from the embedding alone)
        │
        ├── top200_genes_drvi_nonimm_64.tsv       # 03_3  ranked gene list per dimension-direction
        └── enrichr_top200_drvi_nonimm_64.tsv     # 03_3  every tested term (the barplot cache)

## Object conventions (carried through every step)

- `.X` = **raw integer counts** in `shiao_nonimm_raw.h5ad`, **scran log1p-normalized** from step 3
  onward (already in log space - do **not** re-`log1p`). `.layers['counts']` = raw integer counts,
  kept unchanged throughout, and the layer DRVI trains on.
- `var_names` = gene symbols (Ensembl in `.var['gene_ids']`).
- `batch_key = 'cohort'` (patient), `label_key = 'cell_type'` (CellTypist, 18 labels observed here).
- `fraction` is **constant** (`non_imm`) by construction and is no longer usable as a covariate.
  `dataset_origin` is carried along untouched but plays no role in this phase.

## Key parameters (verbatim, ready for Materials & Methods)

**Subsetting and re-filtering (step 1)**
- Subset: `fraction == 'non_imm'` (biological lineage from CellTypist, recoded in 01_4).
- `.X` reset to `.layers['counts']`; `obsm`, `varm`, `obsp`, the derived `uns` entries (including
  the `*_colors` palettes) and the derived `obs`/`var` columns removed; unused categories dropped.
- Cell filter (keep iff): `pct_counts_mt < 10` **and** `n_genes_by_counts > 100` - the same
  thresholds as 01_2, re-applied rather than assumed.
- Gene filter: `sc.pp.filter_genes(min_cells=3)`. Stricter than 01_2's `min_cells=1`: on 176k cells
  a gene seen in one or two of them only inflates the HVG ranking's denominator. Zero-count genes
  would in any case break scran.
- Scrublet is **not** re-run: doublet detection was done per biopsy on the full object, where the
  mixture of lineages is what makes a doublet detectable at all.
- Batch census reported; `DROP_SMALL_COHORTS = False` (no patient dropped by default).

**Normalization (step 3)** - identical to 01_3
- `scib.preprocessing.normalize`: `min_mean=0.1`, `log=True`, `precluster=True` (leiden
  `quickCluster`), `sparsify=False`; `SEED=0`.
- Re-estimated rather than inherited: scran deconvolves size factors from pools of cells, so
  removing 443,083 immune cells changes every estimate.

**Cell cycle (step 4)** - identical to 01_4
- `sc.tl.score_genes_cell_cycle`, Tirosh/Regev 97-gene signature (43 S + 54 G2M), `random_state=0`.
- Re-scored rather than inherited: the reference gene sets are sampled from expression bins built
  on the current object, so the same signature gives different scores on a different population.

**Feature selection + reduction (step 5)** - identical to 01_5
- `scib.preprocessing.reduce_data`: `batch_key='cohort'`, `flavor='cell_ranger'`,
  `n_top_genes=2000`, `n_bins=20`, PCA 50 comps (`svd_solver='arpack'`), neighbours + UMAP on
  `X_pca`. HVG selection is per-patient then merged (`overwrite_hvg=True`).
- The overlap with the phase-01 HVG list is printed by the script and plotted in step 7: a low
  overlap is the expected outcome and the justification for the whole phase.

**Clustering (step 6)** - one deliberate difference from 01_5
- `scib.clustering.cluster_optimal_resolution`, `label_key='cell_type'`, leiden
  (`flavor='igraph'`, `n_iterations=2`), selected by max NMI; every per-resolution column
  (`optscib_nonimm_leiden_<res>`) is kept.
- Resolution grid **0.1-2.0**, the full scib range, not capped at 1.0 as in 01_5. That cap was
  tuned on the 48 labels of the full object; here only 18 are observed but they carry finer
  sub-structure (fibroblast, luminal and vascular subtypes), so the optimum can sit above 1.0.
  The script warns if the optimum lands on the last grid point.

## Established numbers

Measured on `shiao.h5ad` before this phase runs:

- Non-immune cells: **176,610** of 619,693 (28.5%).
- All **34** cohorts are represented; the smallest is Patient30 with **314** cells, then Patient66
  (467) and Patient20 (479).
- **18** of the 48 CellTypist labels observed in the full object have non-immune cells here. The
  five largest: Fibro-matrix (48,181), Lumsec-prol (24,647), LummHR-major (18,401), Lumsec-basal
  (17,729), Vas-venous (15,869).

- After the filters: the cell filter removes no cell, `min_cells=3` leaves **27,853** genes of
  30,869, and dropping Patient01 and Patient30 (incomplete timepoints) brings the object to
  **174,487 cells x 27,826 genes** over **32** cohorts - the definitive `shiao_nonimm.h5ad`.

To be filled in from the run: HVG overlap with 01_5, selected leiden resolution and cluster
count.

## Critical methodological notes

- The subset is taken from the **unintegrated** object, not from a phase-02 integration output.
  DRVI performs its own batch correction from raw counts in 03_2, so an already-corrected matrix
  would be the wrong input; the benchmark of phase 02 informs how to read the batch structure, it
  does not feed this chain.
- CellTypist is **not** re-run: `cell_type` is inherited from 01_4, and sub-annotation of the
  compartment belongs to 03_2. `fraction_reassignment` is likewise meaningless here, `fraction`
  being constant after the subset.
- Step 1 writes through `02_integration_benchmark/utils/h5ad_compat.py`. `shiao.h5ad` was written
  by a pandas-3 stack, so its `obs`/`var` indices come back as pandas `StringArray`s, which anndata
  refuses to write and no reader older than 0.11 understands - and DRVI's environment is one of
  those. The helper downcasts them once, here, and every later step inherits the fix.
- Removing the immune cells removes the largest source of variance in the object. Everything
  downstream of that - size factors, HVGs, PCA, neighbours, clusters - has to be recomputed for
  this reason and no other; inheriting any of it would silently describe the wrong population.
