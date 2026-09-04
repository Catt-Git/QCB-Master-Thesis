# 04_drvi_epithelial

Fifth phase of the thesis, and the second one inside **Part 2 (biological interpretation)**.
Phase 03 isolated the **non-immune** compartment and read DRVI's latent dimensions on it. That
compartment turned out to be 57% stromal and vascular - Fibro-matrix alone is its largest label,
with 48,181 cells - so its HVGs, its PCA and its latent dimensions describe fibroblasts and
endothelium at least as much as epithelium. This phase goes one lineage deeper and keeps the
**epithelial** cells only.

The compartment is defined by the `cell_type` label: the **11 epithelial labels** of the
CellTypist model annotated in 01_4 (`Cells_Adult_Breast.pkl`, Kumar et al. 2023) - `LummHR-*`
(luminal hormone-responsive), `Lumsec-*` (luminal secretory) and `basal` (basal/myoepithelial).
It is the same group of labels that `01_4/fraction_reassignment.py` counts as non-immune, so the
epithelial set is by construction *inside* the compartment of 03; the notebook asserts exactly
that against `fraction` before subsetting.

Input from phase 01: `shiao.h5ad` (619,693 cells x 30,869 genes) - **not** the 03 object, so this
phase does not depend on that chain.
Integration `batch_key = 'cohort'`, biological `label_key = 'cell_type'` (10 labels observed).

## ⚠ What this phase can and cannot say

`cell_type` comes from CellTypist with `Cells_Adult_Breast.pkl`, derived from the Kumar et al.
2023 **normal** adult breast atlas. That model **has no malignant class**: TNBC cells are assigned
to the nearest normal state. The epithelial subset therefore mixes normal and malignant epithelium
in a way the labels cannot distinguish, and **no CNV inference has been run**.

Every state reported here is an **epithelial state, not a tumour cell state**. The caveat is not
left to this README: `signature_common.CAVEAT` is written as `#` comment lines into the header of
every `.csv` in `tables/`, and `CAVEAT_SHORT` is drawn as a footnote on every figure that shows
cells or states, so it travels with the output when it is pulled out of the repo. Read the tables
back with `pd.read_csv(path, comment='#', index_col=0)`. The two 04_3 figures are the exception
(`caveat=False`): they describe the gene-set collection itself, with no cell in them, so there is
nothing there for the caveat to qualify.

Any step that would require malignant status to be meaningful is not in this phase, and the
question "is this a tumour cell state?" cannot be answered from these data.

> **This caveat is a statement about this phase, and it stays true of it.** The question it
> refuses is answered in [05_drvi_tumoral_epi](../05_drvi_tumoral_epi/), which runs inferCNV
> against a per-patient immune reference, calls `malignant` / `non_malignant` per cell, re-runs
> the 01_4 CellTypist annotation on the non-malignant cells only, and then redoes this
> phase's procedure on that basis. Phase 05 does not modify anything here and shares no output
> with it: 04 and 05 are two branches off `shiao.h5ad`, and every number and figure in this
> phase comes from the mixed epithelium described above. Read them as such, and read 05 for the
> version that knows which cells are tumour.

## Repository layout

```
04_drvi_epithelial/
├── README.md
├── signature_interpretation_all.sh   # phase-level driver for 04_3 -> 04_7
├── utils/
│   ├── signature_common.py           # paths, writers, DRVI accessor, embeddings, target region
│   └── sig_collections.py            # the scie and emt collections: lists, axes, target region
├── 04_1_subsetting/                  # from shiao.h5ad to the epithelial object
├── 04_2_drvi_run/                    # DRVI on that object, n_latent 64
├── 04_3_signatures/                  # the lists of one collection -> .gmt, coverage, Jaccard
├── 04_4_cytotrace2/                  # per-patient potency (own conda env, see below)
├── 04_5_cell_first/                  # Route A
├── 04_6_factor_first/                # Route B
├── 04_7_convergence/                 # Route C, the main result
├── 04_8_cycle_confound/              # how much of "stemness" is the cycle - a diagnostic
├── 04_9_embedding_control/           # Route A on Harmony instead of DRVI - a control
├── tables/{scie,emt}/                # every result table of 04_3 - 04_7, one folder per collection
└── figures/                          # one folder per step, then one per collection
```

`utils/` follows 00 and 02: helpers shared by several steps live there, imported with the idiom
`02_2_integration/run_integration.py` uses. `tables/` is phase-level for the same reason
`figures/` is - 04_7 reads what 04_5 and 04_6 wrote, so a per-step `tables/` would mean steps
reaching into each other's folders.

### The two collections

04_3 - 04_7 are **one procedure applied to two independent bodies of prior knowledge**. The step
folders are one per *method step*, not one per readout; which lists are being interpreted is a
flag, `--collection`, declared in `utils/sig_collections.py`:

| | `scie` (the default) | `emt` |
|---|---|---|
| question | is there an epithelial state that is stem-like **and** immune-evasive? | which cells sit in the **hybrid**, partial-EMT state? |
| lists | 11, on `immune` / `stemness`, plus CytoTRACE2 | 9, on `epithelial` / `hybrid` / `mesenchymal`, plus a derived E-to-M score per list version |
| primary | no primary stemness list; `IMMUNOGENIC_CONSENSUS` is the primary immune one | list **B**; A and C are robustness replicates of it |
| hybrid lists | - | scored and reported, but they **validate** the call rather than making it - see `sig_collections.py` |
| target region | stem-**high** x immunogenic-**low** | epithelial **high** x mesenchymal **high**, i.e. co-expression |
| named risks | cell cycle, sequencing depth | cell cycle, fibroblast ambient RNA / doublets |

They share no output. Every table goes to `tables/<collection>/<name>_<collection>_<run_id>.csv`
and every figure to `figures/<step>/<collection>/<name>_<collection>_<run_id>.png`, and **04_6
corrects its FDR inside one collection**, so running one cannot move a single number of the
other. `04_4_cytotrace2` is the one step that is not collection-scoped: it computes a
measurement, and it is the `scie` collection that declares it wants to use it.

Adding a third collection means appending a `Collection` to `sig_collections.py`. No step script
changes.

## 04_1_subsetting

A full re-run of the phase-01 pre-processing on the subset, exactly as in 03_1 and for the same
reason: scran size factors, HVGs, PCA, neighbours, UMAP and leiden were all computed on
populations of which the epithelium is a minority (12% of the full object, 43% of the 03
compartment), and none of them describe these cells. Only the raw counts and the metadata carry
over.

The split between notebooks and scripts follows 01 and 03: notebooks where decisions are taken by
looking at figures, scripts where the computation is deterministic and long.

### Execution order

| # | File | What it does | Where |
|---|------|--------------|-------|
| 1 | `subset_and_qc.ipynb` | Subset `cell_type` on the 11 epithelial labels; write `compartment = 'epi'`; restore `.X` = raw counts; drop every phase-01 derived slot (PCA/UMAP/neighbours/HVG/size factors/leiden/cell cycle) and the emptied categories; recompute QC metrics; re-apply the cell thresholds; gene filter `min_cells=3`; batch census and cohort drops. | local (notebook) |
| 2 | `subsetting_all.sh` | Driver for steps 3-6: runs them in sequence, resumes from the last completed one, logs to `logs/`. | local |
| 3 | `scran_norm_epi.py` | scran normalization re-estimated on the subset (scib: `quickCluster → computeSumFactors → logNormCounts`); adds `size_factors`, clears `.raw`, casts `.X` to float32. | local |
| 4 | `cell_cycle_score_epi.py` | Tirosh/Regev cell-cycle re-scoring → `S_score`, `G2M_score`, `phase`. | local |
| 5 | `reduce_data_epi.py` | Batch-aware HVG (`cohort`, 2000) + PCA(50) + neighbours + UMAP; writes the HVG list **and** the 2,000-gene DRVI input, plus the overlap with both the 01_5 and the 03_1 HVG lists. | local |
| 6 | `clustering_epi.py` | Leiden optimal-resolution sweep (0.1-2.0) vs `cell_type` NMI → `shiao_epi.h5ad`. | local |
| 7 | `visualization_epi.ipynb` | Figures on the final object: composition, cell cycle, normalization, HVG overlap, UMAPs, leiden sweep. Read-only. | local (notebook) |

> **Note on step 2.** The scripts are duplicated from 03_1 (which duplicated 01_3/01_4/01_5)
> rather than called with different paths, so this phase reads as a self-contained Materials &
> Methods section and the earlier scripts stay frozen. The parameters are identical except where
> stated below.

> **Note on where.** ~74k cells, less than half of 03: the whole chain runs locally in sequence
> and there is no SLURM wrapper.

### QC figures (`figures/04_1_subset_and_qc/`, from step 1)

Same panels as 01_2 (`n_genes_by_counts`, `total_counts`, `pct_counts_mt`) plus
`pct_counts_ribo`, at **three rounds** - `before_filter`, `after_cell_and_gene_filter`,
`after_cohort_drop` - and in **two versions** each:

- `violin_epi_qc_<round>_unintegrated.png` - clean, the distribution alone;
- `violin_epi_qc_<round>_unintegrated_jitter.png` - the same with the cells drawn on top.

The jittered version strips a random 20,000-cell subsample (fixed seed), not the whole
compartment: at 75k points the strip becomes a solid black band that hides the violin it is
supposed to annotate. Plus `scatter_epi_counts_vs_genes_unintegrated.png`, the complexity
scatter coloured by mitochondrial content.

### Usage

```bash
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets

# 1. the notebook, interactively (kernel: benchmark-py-r)
#    -> $DATA_DIR/04_epi/shiao_epi_raw.h5ad

# 2. the headless chain
cd 04_drvi_epithelial/04_1_subsetting
./subsetting_all.sh                 # resumes: a step whose output exists is skipped
./subsetting_all.sh --force         # re-run everything
./subsetting_all.sh --dry-run       # print what would run
./subsetting_all.sh reduce cluster  # only the named steps (norm, cc, reduce, cluster)

# 3. visualization_epi.ipynb, interactively
```

`CC_GENES` defaults to `$DATA_DIR/regev_lab_cell_cycle_genes.txt`; export it to override.
The active conda environment must be `benchmark-py-r` (scib + rpy2 + R scran): the driver does
not activate it.

## 04_2_drvi_run

DRVI on the object 04_1 produced: the question of 03_2, one lineage deeper.

A notebook, `drvi_epi.ipynb`, for the same reason DRVI is a notebook in 02_2 and 03_2: the
latent size is judged by eye, from how many dimensions vanish. **`N_LATENT = 64` is the run of
this phase**, the same size as 03_2. 32 was tried first - half the cells of 03 (74,441 vs
176,610) and 10 labels of a single lineage instead of 18 across four suggested a smaller space
would do - and **nothing vanished: 0 of 32 dimensions**, which is the signal that the space is
too tight, so the phase settled at 64. `N_LATENT` in the configuration cell is the only thing
to edit: the run id (`drvi_epi_<N>`), the model, the embedding and the figure folder all follow
from it, so the 32 run is still on disk next to this one, not overwritten.

Three differences from the 03_2 notebook:

- the input is `shiao_epi_hvg_2k.h5ad` (74,441 cells x 2,000 HVGs selected **on this
  compartment**), not the non-immune `shiao_nonimm_hvg_2k.h5ad`;
- `N_LATENT = 64`, the same as 03_2, after 32 saturated (above);
- `compartment` is constant (`epi`) after the subset and so is `fraction` (`non_imm`), so
  neither is plotted, and `dataset_origin` (the CD45 sort) is left out as everywhere in Part 2.

Nothing is inherited from 03 - not the model, not the HVGs, not the latent size. As in 04_1,
the subset comes from the **unintegrated** `shiao.h5ad` chain and not from the 03 object.

**Model.** `DRVI` with `n_latent=64`, `encoder_dims=[256, 128]`, `decoder_dims=[128, 256]`
(against the `[128, 128]` defaults) and `dispersion='gene-batch'`, trained on `layers['counts']`
with `batch_key='cohort'`, `SEED=123`, up to 400 epochs with early stopping after 50 without
improvement (scvi-tools' default 0.9/0.1 split). `dispersion='gene-batch'` fits a dispersion per
gene *and* per batch, which is why 04_1 dropped the cohorts under 200 cells. The lightest DRVI run
of the thesis by cell count: 74,441 cells against the 176,610 of 03_2, at the same 64 dimensions.

**Outputs.** `model_<run_id>.pt` and `embed_<run_id>.h5ad` (the latent space with the dimension
stats and the OOD/IND scores, what 04_3 reads), plus `shiao_epi_<run_id>.h5ad` - the definitive
04_1 object, all genes, with the latent space added as `obsm['X_drvi']`. The embedding carries
no genes by construction, so that third file is the compartment itself in the DRVI space, for
any downstream step that needs genes and latent coordinates in the same object.

`OVERWRITE = False` in the configuration cell reuses the model and embedding already on disk, so
re-running the notebook to redraw a figure does not retrain.

### The same run headless (`run_drvi_epi.py`, `submit_drvi_epi.slurm`)

Judging the latent size needs a pair of eyes; training at a chosen size does not, and it is the
long half of the phase. `run_drvi_epi.py` is the notebook without the kernel - same input, same
architecture, same `SEED`, same early stopping, the same three outputs and the same 25 figures -
so a size can be trained wherever it is convenient and the notebook re-opened afterwards with
`OVERWRITE = False`, which reads model and embedding from disk instead of recomputing them. As
in `subsetting_all.sh`, an output already there is reported `[have]` and reused, so a crash
after training never costs the training.

```bash
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
cd 04_drvi_epithelial/04_2_drvi_run

python run_drvi_epi.py                 # n_latent 64, the run of this phase
python run_drvi_epi.py --n-latent 32   # another size, side by side with it
python run_drvi_epi.py --overwrite     # retrain and rewrite everything

# on the HPC (CPU there, hence partition `long`; 48G is enough for 74k cells)
mkdir -p logs
sbatch --export=ALL,DATA_DIR=$DATA_DIR submit_drvi_epi.slurm --n-latent 64
```

`submit_drvi_epi.slurm` passes everything after the script path to the script unchanged,
activates `$DRVI_ENV` (default `catalano_env`, which pins the same `drvi-py==0.2.7` as the local
`benchmark-py-r`, so the model it writes loads back here) and logs to `04_2_drvi_run/logs/`. The
inputs to copy up first are `$DATA_DIR/04_epi/shiao_epi_hvg_2k.h5ad` and - only if the downstream
object is wanted there - `shiao_epi.h5ad`, by hand:
`02_integration_benchmark/utils/sync_to_cluster.sh` only walks the phase-02 grid. A missing
`shiao_epi.h5ad` is announced at the start and costs nothing else. Only `embed_<run_id>.h5ad` has
to come back, since that is all 04_3 reads.

The number to read in the log before deciding whether to re-run at another size is the vanished
count (`n vanished / 64 latent dimensions`). None vanished means the size is too tight and a
larger one is worth a run - which is exactly what the first attempt at 32 reported (0 / 32).

### Figures (`figures/04_2_drvi_epi_64/`, notebook or `run_drvi_epi.py`)

Same set as the DRVI figures of 03_2, each name suffixed with the run id, next to the `04_1_*`
folders:

- `umap_<key>_<run_id>.png` - one UMAP of the DRVI space per metadata/QC key, plus
  `umap_combined_` and `umap_per_cell_type_` (one panel per label: three of them are 81% of the
  compartment and dominate the single panel). Each has its direct counterpart in
  `figures/04_1_visualization/umap_*_epi_unintegrated.png` - same cells in the *unintegrated* PCA
  space, same keys, same palettes, so only the space changes.
- `latent_dimension_stats[_rmVanished]_<run_id>.png` - per-dimension reconstruction effect, with
  and without the vanished dimensions: the plot behind the latent-size choice, and a diagnostic
  of the run rather than a filter.
- `latent_dims_in_umap_` and `latent_dims_in_heatmap_<key>_` - each dimension on the UMAP, and how
  the dimensions respond to `cell_type` (also sorted by label), `cohort`, `treatment`, `response`,
  `phase`.
- `ood_*_` / `ind_linear_weighted_mean_` - interpretability scores. OOD comes from the decoder
  reconstructions (favours the genes *specific* to a dimension); IND averages each factor's effect
  over all cells, so a shared gene stays high in every dimension it belongs to.

## 04_3 - 04_7: signature interpretation

Second downstream stage of the phase, and the one that turns the 04_2 latent space into biology.
04_2 stopped at *which* genes each dimension moves; this stage asks what those dimensions
**mean**, against a curated collection of gene signatures, and - unlike 03_3, which was
enrichment alone - it asks in **both directions** and reports where the two agree.

DRVI returns dimensions that are data-driven and unnamed; the signature collection returns names
with no guarantee that the matching axis exists in this dataset. Interpretation is a mapping
between two sets **neither of which is ground truth**, and it can be traversed either way - but
the two directions answer different questions and fail differently.

| | Route A, cell-first (04_5) | Route B, factor-first (04_6) |
|---|---|---|
| unit of analysis | the **cell** | the **gene** |
| question | do the cells prior knowledge calls stem-like, immune-evasive or hybrid-EMT sit anywhere in particular along a latent dimension? | what program does this dimension encode, irrespective of what I was looking for? |
| why it is needed | the **only** route that assigns cells to states - the project's deliverable | the **only** discovery route, and the corrective for A's confirmation bias |
| how it fails | noisy per-cell scores, correlated with depth and cycle and with each other; confirmatory by construction | gene overlap with no cellular counterpart; long signatures under-enrich against a truncated top-gene list |

Route B is also what makes the choice of DRVI pay off: its additive decoder gives every dimension
a directly readable gene-level footprint. Without it any of the higher-scoring phase-02 methods
would have done.

**Agreement between the routes is the criterion for calling a dimension a genuine cell state.**
The two failure modes do not overlap - a dimension can pass A by coincidence among correlated
scores, or pass B by gene-set overlap with no cells behind it, but it is unlikely to pass both for
the wrong reason. Disagreement is informative rather than a failure, and **all three categories
are reported separately. Nothing is promoted on a single route.**

### Execution order

| # | Step | What it does | Where |
|---|------|--------------|-------|
| 1 | `04_3_signatures/build_signatures_epi.py` | Ingest the collection's text files; write `signatures_<collection>.gmt` with the provenance in the description field; coverage table; pairwise Jaccard. | local |
| 2 | `04_4_cytotrace2/cytotrace2_epi.py` | CytoTRACE2 potency, one patient at a time, from raw counts; concatenate. **Not in the default chain** - needs the `cytotrace2-py` env. | local |
| 3 | `04_5_cell_first/cell_first_epi.py` | Route A: scoring, within-patient standardisation, derived readouts, confounder table, target-region definitions, dimension x signature association. | local |
| 4 | `04_3_signatures/signature_composition_epi.py` | Which genes actually carry each score, and whether they are detectable: per-gene contribution, effective gene count, flags. Feeds nothing - a diagnostic. Lives in 04_3 but runs after 04_5 so it can use the real scores. | local |
| 5 | `04_6_factor_first/factor_first_epi.py` | Route B: top genes per dimension **and direction**, offline ORA, signed significance matrix. | local |
| 6 | `04_7_convergence/convergence_epi.py` | Route C: the convergence table and the comparative figures. | local |
| 7 | `04_8_cycle_confound/cycle_confound_epi.py` | How much of "stemness" is the cell cycle: gene content of each list, variance it explains, what the consensus vote does to it, and the cycle loading of every dimension. Feeds nothing - a diagnostic. | local |

```bash
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
cd 04_drvi_epithelial

./signature_interpretation_all.sh                  # steps 1, 3, 4, 5, 6, 7 on scie - resuming
./signature_interpretation_all.sh --collection emt # the same six steps on the EMT lists
./signature_interpretation_all.sh --force          # re-run everything
./signature_interpretation_all.sh --dry-run        # print what would run
./signature_interpretation_all.sh cellfirst convergence   # only the named steps
PYTHON=~/miniconda3/envs/cytotrace2-py/bin/python \
  ./signature_interpretation_all.sh cytotrace      # the one step under a different interpreter
```

Step names: `signatures`, `cytotrace`, `cellfirst`, `factorfirst`, `convergence`. A step whose
output already exists is reported `[have]` and skipped, as in `04_1/subsetting_all.sh`; the check
is existence only, so delete a file truncated by a crash before resuming. Logs go to
`04_drvi_epithelial/logs/`. Environment `benchmark-py-r`; the driver does not activate it.

### Inputs (verified, not assumed)

| Object | Shape | What was checked |
|---|---|---|
| `04_epi/shiao_epi.h5ad` | 74,441 x 26,371 | `.X` log-normalised (`uns['log1p']` present, asserted `< 50` at runtime); `layers['counts']` raw integers; all required `.obs` keys present |
| `04_epi/shiao_epi_hvg_2k.h5ad` | 74,441 x 2,000 | the DRVI training feature set, read for its `var_names` alone - it is the ORA background |
| `04_epi/embed_drvi_epi_64.h5ad` | 74,441 x **64** | latent space in `.X`; per-dimension stats in `.var`; per-gene scores in `.varm` |

**`.obs` keys used.** `cohort` (patient, 29), `cell_type` (CellTypist, 10 labels), `treatment`
(`BASE` / `PD1` / `RTPD1`), `response` (`NR` / `R1` / `R2`), `n_genes_by_counts`, `pct_counts_mt`,
`S_score`, `G2M_score`, `phase`.

> **There is no separate `timepoint` column.** The timepoint axis is carried by `treatment`, with
> the per-biopsy detail in `sample` (`P02_A_P`, `P02_B_P`, …). Nothing in this stage needed a
> timepoint, so nothing was invented.

**Cell cycle.** `S_score` and `G2M_score` are **reused from 04_1** as they stand: not recomputed,
and the cycle is **not** regressed out. They serve only as covariates in the confounder table and
as the stratifier of the G1-only check.

**Latent space.** The embedding's `.X` **is** the latent representation (74,441 x 64) - there is no
`.obsm` key on it. The same coordinates also live as `obsm['X_drvi']` in
`shiao_epi_drvi_epi_64.h5ad`, which this stage does not need.

**Vanished dimensions are not pruned.** `signature_common.PRUNE_VANISHED = False`: 04_5 - 04_7
run on **all 64 dimensions** and **all 128 dimension-directions**. `03_2`'s
`plot_pruned_umap_nonimm.py` measured what pruning does to a DRVI space of this kind and found it
inert, so pruning buys no cleanliness while deciding, ahead of the analysis, which axes are
allowed to mean something; the cost is a larger FDR denominator in Route B. The flags are still
read programmatically from `var['vanished']` and `var['vanished_*_direction']` - never from a plot
- and reported instead: `dimension_row_order_*.csv` carries a `vanished` column and
`convergence_*.csv` a `dimension_vanished` one. Set the constant to `True` to prune again.

**Cell identity.** All three objects hold the same 74,441 cells **in the same order**; asserted at
runtime in both routes rather than assumed. Nothing needed reindexing.

## 04_3_signatures

Ten files, one gene symbol per line, in `$DATA_DIR/signatures/` (override with
`SIGNATURE_DIR`). Ingested into `tables/<collection>/signatures_<collection>.gmt`, whose **description field carries the
provenance**, so the collection that feeds both routes is also the Appendix table. The `.gmt`
holds the **mapped** genes, i.e. the collection as actually used, so the Appendix and the tested
sets cannot drift apart.

Several of the files are **CRLF** and at least one carries a **UTF-8 BOM**; the reader opens
with `encoding='utf-8-sig'` and strips whitespace. Without that, `B2M\r` matches no `var_name` and
the immune signatures silently score zero.

The registry - name, file on disk, axis, provenance - is `signature_common.SIGNATURES`:

| Signature | Axis | Provenance |
|---|---|---|
| `HALLMARK_IFNA` | immune | Interferon Alpha response, MSigDB Hallmark |
| `HALLMARK_IFNG` | immune | Interferon Gamma response, MSigDB Hallmark |
| `ISDS` | immune | IFN-Stem Cell-Down signature, PMC5481166 |
| `KEGG_APM` | immune | Antigen Presentation Machinery, KEGG |
| **`IMMUNOGENIC_CONSENSUS`** | immune | curated by the lab from the APM signature, immunogenic genes and immunomodulators only |
| `BENPORATH_ES1` | stemness | MSigDB `BENPORATH_ES_1` |
| `ESC_ASSOU` | stemness | PMC1906587, Table S3 |
| `ESC_WONG` | stemness | MSigDB `WONG_EMBRYONIC_STEM_CELL_CORE` |
| `LIM_STEM` | stemness | MSigDB `LIM_MAMMARY_STEM_CELL_UP` |
| `FMASC` | stemness | fetal mammary stem cells, PMC3277444, Suppl. Table 2 |

Two file names differ from the signature name: `Immune consensus.txt` → `IMMUNOGENIC_CONSENSUS`,
`fMaSC.txt` → `FMASC`.

### Why `EMP` was dropped

The stemness axis carried six lists until now. The sixth, `EMP` (Embryonic Multipotent
Progenitors, PMID 29784918, `EMP.txt`, 15 symbols), has been removed; five remain. `EMP.txt` is
still on disk and the restoring change is one line in `utils/sig_collections.py`, where the same
reasoning is repeated for whoever reads the code rather than this file.

Every number below was measured on the six-list run, and the re-run has overwritten its tables -
**this section and that code comment are the record of it.** (The last *committed* with-`EMP`
tables are the earlier 104-row ones from before `PRUNE_VANISHED` was turned off, so they are not
a like-for-like comparison.)

**It is not a coverage failure.** 13 of the 15 symbols map, 87 %, comfortably above the 60 %
floor. The problem is what those 13 measure on *this* data.

*Three genes carry 92 % of the score.* Of the total EMP signal across the 74,441 cells, `RPSA`
carries **64.8 %**, `FN1` 16.7 % and `MFAP2` 10.6 %; the remaining ten together carry 7.9 %.
`RPSA` is a 40S ribosomal protein, ubiquitous and depth-sensitive, and the score follows it -
Spearman **0.70** with `RPSA` expression against **0.11** with `SOX11`, the progenitor marker the
list is actually named for. `FN1` and `MFAP2` are ECM genes, which is why the EMP score tracked
the mesenchymal axis of the `emt` collection (ρ 0.27-0.31 on `EMT_A/B/C_MESENCHYMAL`) as closely
as it tracked the other stemness lists (ρ 0.21-0.27) - an unwanted coupling between two
collections that are supposed to be independent bodies of prior knowledge.

*The genes that make it an EMP list are below the droplet noise floor.* Detection rates: `NDNF`
0.2 % of cells, `IGF2BP1` 0.3 %, `FRAS1` 1.1 %, `EPHA7` 1.2 %, `UNC5B` 1.7 %, `GPC3` 2.0 %,
`CCND2` 2.3 %, `LEF1` 2.4 %, `SOX11` 3.3 %. **20.9 % of cells have zero counts across all 13
genes.** Its `hvg_fraction_of_mapped` of 0.46 was the highest of the stemness lists and looked
like a point in its favour; it is 6 genes out of 13, too few for the Route B hypergeometric test
to be informative, which is exactly what happened.

*It contributed no evidence.* Across all 128 dimension-directions of
`convergence_scie_drvi_epi_64.csv`, `EMP` was **never significant on either route**: best-of-row
3× on Route A (ρ 0.019 / 0.054 / 0.056, bar 0.20) and 3× on Route B (FDR 0.11 / 0.87 / 0.11,
bar 0.05). Its maximum |ρ| over all dimensions was **0.104**, against 0.30-0.35 for every other
stemness readout and 0.29 for CytoTRACE2. It produced no `convergent` verdict and changed none.

*And it cost stability.* In `quadrant_stability_*.csv` it sat at Jaccard 0.11-0.14 with the ESC
block and 0.14 with CytoTRACE2. Removing it raises the **median pairwise Jaccard of the stemness
quadrants from 0.184 to 0.255**. 374 cells entered the seven-readout consensus on the strength of
the EMP vote alone.

*What is deliberately not claimed.* `LIM_STEM` is the other readout that disagrees with the ESC
block - it holds the minimum pairwise Jaccard (0.073, against `ESC_WONG`) with or without EMP -
and it is **kept**. It maps 465 of 479 symbols, reaches |ρ| 0.31 and is the best Route A readout
of 26 dimension-directions: its disagreement is a documentable biological one, not an artefact of
what happens to be measurable. Only EMP fails on measurability.

### What is inside each score (`signature_composition_epi.py`)

`coverage_*.csv` answers *how many symbols map*. It does not answer *what the score then
measures*, and the two come apart badly. `sc.tl.score_genes` averages the mapped genes, so a gene
expressed an order of magnitude above the rest of its list sets the average no matter how many
other genes are in it, and a gene under the droplet detection floor contributes nothing but its
zeros. A list can map 87 % of its symbols and still produce a score that is, on this object, a
proxy for one ribosomal gene - that is exactly what `EMP` was, and it is why this step exists.

Two tables come out of it. `signature_gene_contribution_*.csv`, one row per (signature, gene):

| Column | What it says |
|---|---|
| `detection_rate` | fraction of cells with a non-zero count. Under 1 % the gene is at the noise floor and measures nothing per cell |
| `share_of_mean_signal` | the gene's share of its list's summed expression - how much of the score's **level** it sets. Sums to 1 |
| `share_of_score_variance` | `Cov(x_g / n, m) / Var(m)` with `m` the list mean - how much of the score's **spread across cells** it sets. Sums to 1, and is negative for a gene that moves against its own list |
| `spearman_with_score` | the gene against the score Route A actually uses, joined from 04_5 |

`share_of_score_variance` is the sharper of the two: every quadrant and every Spearman downstream
is built on the *ranking* of cells, not on the level.

And `signature_concentration_*.csv`, one row per signature, whose headline column is
**`effective_n_genes`** - inverse Simpson on the contribution shares, `1 / Σ share²`, i.e. the
number of genes the score behaves as if it had. A list of equally-contributing genes scores its
own length; `EMP` scored **2.17 of 13**.

Three flags, all reported and none enforced - this step drops nothing and stops nothing:

| Flag | Threshold | Reading |
|---|---|---|
| `dominated_by_one_gene` | top gene > 30 % of the score's variance | the other genes are decoration |
| `low_effective_n` | `effective_n_genes` < 10 | too few effective genes to average away any single one |
| `undetectable_in_many_cells` | > 10 % of cells with zero counts across the whole list | for those cells the score is its control set, not a measurement |

**`effective_n_fraction` is reported and deliberately not flagged.** Thresholding it looked like
the right way to measure the *shape* of a list independently of its length, and it does not work:
every long list on this object sits at 0.18-0.28 - `FMASC` at 0.27 with 457 effective genes,
`ESC_ASSOU` at 0.18 with 146 - because expression is heterogeneous across any few hundred genes,
not because those lists are defective. A threshold there flags the collection's most robust
scores. What separated `EMP` was the **absolute** effective count, 2.17, together with a single
dominant gene, which the two flags above already catch. Read the fraction next to the absolute
number; do not read it alone.

**A flag is not a verdict.** The thresholds are chosen here, not derived from anything, and a
concentrated list can still be the right instrument - `IMMUNOGENIC_CONSENSUS` is 25 curated genes
and is *supposed* to be dominated by the HLA locus. Acting on a flag means editing the registry in
`utils/sig_collections.py` with the reasoning written next to it, which has happened exactly once.

**What it says about the two collections.** On `scie` only `ISDS` is flagged, and mildly: 8.3
effective genes of 27, `CD74` at 13.5 % of its variance. Everything else is clean, the five
stemness lists at 77-457 effective genes. The immune four share their top genes - `CD74`, `HLA-B`,
`HLA-A` in every one of them - which is the Jaccard matrix's nesting story arriving by a second
route. `IMMUNOGENIC_CONSENSUS` behaves as ~10 genes of 25 and is *not* flagged, which is the
intended shape for a curated list built around the HLA locus.

`emt` is the weaker collection by this measure. All nine lists are short, so `low_effective_n`
fires on eight of them and is not informative on its own; what is informative is that
**`EMT_B_HYBRID` and `EMT_A_MESENCHYMAL` trip all three flags** - `YBX1` at 31.9 % and `VIM` at
38.0 % of their variance, with 15.7 % and 25.6 % of cells scoring on nothing - and that
`EMT_A_HYBRID` has the worst detectability in the phase at **36.1 % of cells with zero counts
across the whole list**. That is an independent reason for something 04_5 had already found the
hard way: the first, quantile-band definition of the hybrid state was unstable across list
versions (Jaccard 0.08-0.13) and had to be replaced by co-expression. The hybrid lists are the
least reliable instruments here, and the co-expression definition leans on them least.

For the record, `EMP` scored **2.17 effective genes of 13** and tripped all three flags - the only
list in either collection to do so. Re-deriving that on the current code is a `sc.tl.score_genes`
call on `$DATA_DIR/signatures/EMP.txt`, which is still on disk; it is not in the registry and no
step reads it.

How much of each list is actually measured is in `tables/<collection>/coverage_*.csv`, against two universes:
all genes of the epithelial object (what Route A scores on) and the 2,000 HVGs (the Route B ORA
background). The step **stops** below a 60 % mapping floor (`--allow-low-coverage` to override):
these lists date from 2007-2012 and carry deprecated symbols, and a gene that does not map is
**not measured**, which is not the same as not expressed.

**The immune axis is expected to be LOW** in the states of interest - the project looks for immune
*evasion*, not immunogenicity. `IMMUNOGENIC_CONSENSUS` is the **primary** immune readout, the
other four secondary confirmation. **No stemness signature is primary** and there is deliberately
**no stemness consensus**: the intersection of the three embryonic lists was tested by the lab and
captures proliferation only, which `S_score` / `G2M_score` already cover. Each stemness list is
used on its own and their disagreement is reported rather than averaged away. That axis is
provisional on the lab's side; an update is expected.

**The signatures are not independent** - the immune four are largely nested and the embryonic
stemness lists overlap substantially - so **they must not be treated as independent tests, and the
FDR of 04_6 must not be presented as if they were.** The pairwise overlap is published
(`jaccard_*.csv`, `shared_genes_*.csv`, `jaccard_signature_overlap_*.png`) so the reader can see
how much of any apparent agreement between two readouts is just shared genes.

## 04_4_cytotrace2

CytoTRACE2 is a per-cell predictor of differentiation potency, **not a gene set**. It belongs to
Route A only and has no Route B counterpart. Its value is that it is **independent of the lab's
lists**, so it is the only non-circular evidence on the stemness axis - which matters more than
usual here, because the stemness consensus was deliberately not built and no single stemness
signature is primary.

**Run per patient.** The score is computed within the object it is given, so scoring all 29
cohorts at once would let the ranking absorb batch and "stemness" would become a batch axis. Each
cohort is scored on its own and the results concatenated, so the output cannot be confounded with
`cohort` by construction. Scores are then comparable within a patient and only ordinally across
them - which is how 04_5 uses them anyway, every readout there being standardised within
`(cohort, cell_type)`.

The package's API is not what an `.h5ad` pipeline expects - tab-delimited genes x cells input,
non-log-transformed counts, `species` defaulting to `"mouse"` - and the script's header documents
every point of it; the parameters used are in *Key parameters* below.

**This is its own step because of the environment, not because of the machine.** The dependency
cannot live in `benchmark-py-r`, so it runs under a different interpreter; it was tried on the
cluster first and that turned out to be unnecessary (see *Runtime, measured* below).

### The environment, which is not optional

`cytotrace2-py` 1.1.0.4 declares **`numpy<2.0.0`** as a hard requirement. `pip install
cytotrace2-py` inside `benchmark-py-r` therefore does not add a package, it rolls the analysis
stack back - numpy 2.4 → 1.26, pandas 3.0 → 2.3, scipy 1.18 → 1.17, zarr 3.2 → 3.1, anndata
0.13 → 0.12, scanpy 1.12 → 1.11 - and leaves `fast-array-utils` and `tifffile` unsatisfiable.
It gets its own environment, exactly as scGen does in phase 02:

```bash
conda env create -f environments/cytotrace2-py.yml     # Python 3.11, numpy 1.26, CPU-only torch
conda activate cytotrace2-py
```

The isolation costs nothing here: this step reads `shiao_epi.h5ad`, writes tab-delimited
matrices and writes one `.csv`, and anndata 0.12 reads the 0.13-written file without complaint.
Two details the `.yml` pins for a reason - `setuptools<81`, because `cytotrace2_py` still
imports `pkg_resources`, which setuptools 81 removed (install succeeds, import fails); and the
CPU torch index, because the predictor runs on the CPU and the CUDA wheels would duplicate
~2.5 GB for nothing.

On first use the package downloads **17 model files (~8 MB each) from Google Drive** into its own
`site-packages/cytotrace2_py/resources/models/`, so the first call needs network access.

```bash
conda activate cytotrace2-py
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
cd 04_drvi_epithelial/04_4_cytotrace2
python cytotrace2_epi.py --cores 12          # or --dry-run to export the matrices only
```

### Runtime, measured

The full run, 29 cohorts / 74,441 cells, on 12 cores: **47 min of scoring + 6 min of export**,
0 cells unscored. The cost is almost entirely *fixed per call* - the package reloads its 17
models every time - so it is ~68 s for a 209-cell cohort and 196 s for the 9,243-cell one,
about 0.012 s per cell on top of the constant. Peak memory is a few GB (one cohort densified
at a time, the largest being 9,243 x 26,371) and the temporary `.txt` matrices total ~4 GB and
are deleted at the end unless `--keep-work`. Nothing here needs a scheduler, and the SLURM
wrapper this step used to carry was removed once the local run measured it.

### Result

Run on all 29 cohorts; `cytotrace2_<run_id>.csv` exists. 74,441 cells scored, none missing, none
unknown. Gene mapping is stable across cohorts: 14,043 input gene names map to mouse orthologs
and 14,042 of them are in the model's feature set.

| Potency | cells |
|---|---|
| Differentiated | 53,025 |
| Unipotent | 17,576 |
| Oligopotent | 3,379 |
| Multipotent | 461 |
| Pluripotent / Totipotent | 0 |

**The cycle check passes.** The worry written into the script was that potency would just be
proliferation: it is not. `basal` has the highest mean score (0.224), and `Lumsec-prol` -
the cycling population and by far the largest cell type at 24,442 cells - ranks **5th of 10**
(0.169), below `Lumsec-HLA`, `Lumsec-myo` and `Lumsec-basal`. Full table in
`tables/scie/cytotrace2_by_cell_type_scie_<run_id>.csv`.

> **But it is not depth-free.** CytoTRACE2's own premise is that the number of expressed genes
> tracks potency, so it inherits the coupling the lab's lists have: rho **+0.445** against
> `n_genes_by_counts` in `confounders_<run_id>.csv`, and **+0.39 median rho even within
> `(cohort, cell_type)`** (IQR [+0.16, +0.61], 42% of the 108 groups of >=50 cells above 0.5),
> which is to say the within-stratum standardisation 04_5 applies does **not** remove it.
> 18 of the 29 cohorts also raised the package's `<500 genes per cell` warning, up to 54.2% of
> cells in `Patient64`. CytoTRACE2 is therefore evidence independent *of the lab's lists* -
> which is what the non-circularity claim actually says - but not independent *of sequencing
> depth*: a second opinion on the stemness axis, not an arbiter.
>
> The ranking is the useful part, though: at +0.445 CytoTRACE2 is **less** depth-coupled than
> `ESC_ASSOU` (+0.642), `ESC_WONG` (+0.503) and `BENPORATH_ES1` (+0.484), and its cell-cycle
> coupling (+0.262 on `S_score`) is below theirs too. It is the least confounded of the four
> readouts that agree with each other, not another equally bad one.
>
> Per-patient scoring worked on the batch side but not completely: `cohort` still explains
> R^2 = 0.16 of the raw score, which is exactly what the 04_5 standardisation is for.

### What it changed downstream

> **Historical.** Every number in this subsection was measured on the **seven-readout**
> configuration, i.e. with `EMP` still in the registry. It is the record of what CytoTRACE2
> changed at the time and is left as measured. `EMP` has since been dropped (see *Why `EMP` was
> dropped* under 04_3, and *What dropping `EMP` changed* below); the current tables carry six
> readouts and a ≥ 3-of-6 bar.

04_5, 04_6 and 04_7 were re-run in that order with the `.csv` in place. In short: CytoTRACE2
changed the stemness consensus and changed **no** Route C verdict.

*Route A.* CytoTRACE2 lands in the middle of the readouts on every axis. Its own quadrant is
4,651 cells (6.25% of the compartment), squarely inside the 4,295-5,206 range the six lists
span. On the pairwise Jaccard it sits with the ESC block - 0.41 with `ESC_ASSOU`, 0.39 with
`ESC_WONG`, 0.35 with `BENPORATH_ES1` - and far from `LIM_STEM` (0.14) and `EMP` (0.14); the
median pairwise Jaccard across the 7 is 0.184. It agrees with 57 of the 61 cells the six lists
had called unanimously, which is the strongest single piece of support the stemness axis has.

*The consensus quadrant shrank from 4,017 to 2,854 cells, and not because CytoTRACE2
disagreed.* The rule is `>= ceil(n_readouts/2)`, so the bar moved from 3 votes of 6 to 4 of 7.
At the old bar of 3, the seven readouts would call 4,918 cells. Worth knowing before reading
the drop as a tightening of evidence: it is a change of threshold, and it is arithmetic.

*The cycle check got slightly worse, and that is the honest cost.* The stricter consensus is
more cycling, not less: G1 share of the target quadrant 28.8% → 21.0%, and the Jaccard between
the quadrant and the quadrant recomputed inside G1 alone 0.598 → **0.538**. Raising the vote bar
concentrates on cells that many readouts agree on, and those are disproportionately cycling.

*Route C is unchanged where it can be compared.* Of the 104 rows both the old and the new
`convergence_<run_id>.csv` have, **0 changed verdict**. CytoTRACE2 has no Route B counterpart by
construction, so it enters Route C only through the effect sizes, and it did not move any across
a bar.

### What dropping `EMP` changed

04_3, 04_5, 04_6 and 04_7 were re-run in that order after the removal (see *Why `EMP` was
dropped* under 04_3 for the reasoning). Every table and figure under `tables/scie/` and
`figures/*/scie/` in the repo is from that run. In short: it **changed no Route C verdict**, and
it improved the two things the stemness axis was weakest on.

*Route A - the axis got more coherent.* The **median pairwise Jaccard of the stemness quadrants
rose from 0.184 to 0.255**, and the range tightened to 0.073-0.541. The consensus quadrant grew
from 2,854 to **4,145 cells**, but read that as arithmetic before biology: the bar is
`>= ceil(n_readouts/2)`, so six readouts move it from 4-of-7 to **3-of-6**. Held at a matched bar
of 4 the two consensus sets agree at Jaccard 0.87, and 374 cells had been entering the old
consensus on the EMP vote alone. 12,559 cells are now called by at least one definition, 181 by
all six.

*The cycle check improved, which partly reverses the cost CytoTRACE2 imposed.* G1 share of the
target quadrant **21.0% → 28.5%**, and the Jaccard between the quadrant and the quadrant
recomputed inside G1 alone **0.538 → 0.594**. That is still not a clean pass - the target is
enriched for G2M and S - but it is back above where the seven-readout consensus had pushed it.
The depth check on the immune axis is untouched, as it must be: median `n_genes_by_counts` 1,007
in the `IMMUNOGENIC_CONSENSUS`-low group vs 1,335 in the rest, AUROC 0.560.

*Route B - one fewer gene set in the denominator.* BH now runs over 128 directions x 60 gene sets
= **7,680 pairs** (was 7,808), of which 6,813 return a non-empty overlap and **693 are significant
at global FDR < 0.05**. Loosening the denominator by 128 pairs cannot make an already-significant
result go away, and none did.

*Route C - nothing moved.* The verdict distribution is **identical**: 54 `neither`, 43
`factor_only_candidate_patient_or_technical`, 18 `convergent`, 9 `cell_only_state_not_on_one_axis`,
4 `both_routes_different_family`, and 9 of the 18 convergent rows still carry a confounder flag.
`EMP` had been best-of-row in 6 of the 128 dimension-directions (DR 13-, 19-, 21+, 45+, 54-, 55+)
and significant in none of them; in the new table the second-best readout steps in and all six
keep their verdict. The `stem_rho AND immunogenic_low_rho` test still finds no single axis
carrying both at |ρ| ≥ 0.20.

## 04_5_cell_first (Route A)

Scored on the **unintegrated, all-genes** object. A 150-gene signature reduced to whatever
survived HVG selection is no longer that signature, so `shiao_epi.h5ad` is read here and never
`shiao_epi_hvg_2k.h5ad`.

**Scoring.** `sc.tl.score_genes`, `use_raw=False`, `ctrl_size=len(signature)` (matched per
signature, so a 25-gene and an 800-gene list are not scored against the same control size),
`n_bins=25`, `random_state=0`. The control gene set is sampled **at random**, so the seed is not
optional. A signature under 10 mapped genes is skipped and reported.

**Standardisation.** Every score is z-scored **within `(cohort, cell_type)`**. Absolute scores are
not comparable across patients sequenced at different depths, nor across cell types with different
baseline expression. A stratum with no spread - a single cell, say - gets z = 0, the only honest
value for a stratum of one.

**Quadrant.** stem-high (z ≥ q0.75) and immunogenic-low (z ≤ q0.25) against
`IMMUNOGENIC_CONSENSUS`. Quantile cutoffs rather than fixed z values: the scores are not normal,
and a fixed z would give wildly different group sizes across readouts and make the stability
comparison meaningless. The quadrant is defined **once per stemness readout** (five signatures
plus CytoTRACE2, now present) and the consensus cell set is the majority vote across them; the
**stability of that set across the definitions is itself a reported result**, in
`quadrant_stability_*.csv`, not something to be averaged away.

**The three checks that qualify everything downstream**, all written to `tables/`:

- *confounders* - Spearman of every raw score against depth, mito, `S_score` and `G2M_score`. Long
  embryonic lists largely report how many genes were detected; this is where that is measured
  rather than assumed.
- *is stem-high just cycling?* - the quadrant recomputed **inside G1 alone**, the cycle held out
  rather than regressed out, and the two cell sets compared by Jaccard.
- *is immunogenic-low just shallow?* - the depth of the immunogenic-low group against the rest, by
  Mann-Whitney and by the AUROC of "shallower" predicting the group. This one matters more here
  than anywhere else: immune evasion is defined by the **absence** of signal, which low sequencing
  depth mimics perfectly.

Then the latent space enters: every dimension against every standardised readout
(`dim_signature_spearman_*.csv`), and the consensus quadrant's effect size on each dimension both
ways, AUROC and standardised mean difference (`dim_target_effect_size_*.csv`). The row order those
tables establish goes to `dimension_row_order_*.csv` and is reused by Route B, so the two heatmaps
are directly comparable.

## 04_6_factor_first (Route B)

**Reused from `03_3_enrichment/enrichment_nonimm.ipynb` unchanged:** the `interpretability_scores`
accessor that rebuilds DRVI's genes x dimension-directions table from `embed.varm` alone (no
model, no GPU, no `scvi-tools`); reading each dimension in its **two directions separately**; the
**HVG background**; `N_TOP_GENES = 200`; and the `OOD_combined` scores, which favour the genes
*specific* to a dimension. Both directions and the declared background were therefore already
03_3's behaviour. What changed here:

- **the gene sets are the lab's collection**, with Hallmark 2020 alongside as a sanity-check
  collection so a dimension enriching for nothing custom can still be named;
- **every test runs offline** (`gp.enrich`, hypergeometric) against an explicit background, never
  Enrichr's implicit all-human-genes universe;
- **BH is applied once across all dimension-direction / gene-set pairs**, not per query. 128
  directions x 60 gene sets = **7,680 pairs**, the pairs with no overlap included: they are p = 1
  and cannot become significant, but they belong in the denominator. With 100+ directions a
  per-query FDR is far too permissive;
- **no direction is pruned.** 03_3 let DRVI's accessor drop the directions it had marked vanished;
  here all 2 x 64 are tested (see *Vanished dimensions are not pruned* above);
- the output is a matrix on the **same row order as Route A**, read from
  `tables/<collection>/dimension_row_order_*.csv`, so the two heatmaps are directly comparable.

**The gene universe, for the Methods.** DRVI was trained on the 2,000 batch-aware HVGs of 04_1, so
a gene outside that set could never have entered a top-gene list: the ORA background is the
**training feature set, not the transcriptome**. The consequence is stated rather than hidden -
**each signature is effectively tested in its HVG-restricted form**, which is why
`n_in_hvg_background` is a column of the coverage table, and why a short list can be scored in
full by Route A and still be untestable here. That is a property of ORA, not a contradiction of
the all-genes rule, which governs scoring.

## 04_7_convergence (Route C)

One row per dimension **and direction**, 128 rows. Route A is computed per dimension, so its
direction is the **sign of the Spearman correlation**: a signature correlating positively with
DR 7 is a statement about `DR 7+`. That is what makes the two routes joinable.

Bars: Route A `|ρ| ≥ 0.20`, Route B global FDR < 0.05. Every row gets one verdict:

| verdict | reading |
|---|---|
| **`convergent`** | both routes, same signature family - **a candidate cell state** |
| `both_routes_different_family` | both routes fire, on unrelated signatures |
| `factor_only_candidate_patient_or_technical` | the program is on the axis, no coherent cell group is |
| `cell_only_state_not_on_one_axis` | the cells separate; this axis is not the description |
| `neither` | - |

Every row also carries the confounder flags of the signature it claims (depth or cycle coupling at
`|ρ| ≥ 0.30`, and `immune_low_is_absence_of_signal` where the claim rests on low immunogenicity),
so a convergent row that cannot be read as clean says so in its own line.

The step ends on the project's actual target: whether any single dimension-direction carries
stemness **and** immunogenic-low at the same time, above the same bar. If any does, they go to
`target_axes_*.csv`; if none does, that is a result rather than a gap - the target state is then an
**intersection of two axes** in the latent space rather than a direction of it, which is what the
Route A quadrant assumed by crossing two scores in the first place.

> **Run 04_5 → 04_6 → 04_7 in that order, always.** 04_6 does not merely follow 04_5, it *reads*
> `dimension_row_order_<run_id>.csv` and builds its signed matrix on exactly the dimensions listed
> there. Re-run 04_5 alone and 04_6's matrix is left indexed on the previous dimension list, after
> which 04_7 raises a bare `KeyError` on the first dimension the two disagree about - not a
> warning, a crash mid-table. This had already happened in this phase: the row order moved to 64
> dimensions while 04_6's matrix stayed at the 52 of the run before it, so `convergence_*.csv` on
> disk covered 104 of the 128 rows and could not be regenerated until 04_6 was re-run. The driver's
> `--force` chain does the three in order for this reason; the failure mode is only reachable by
> asking for steps by name.

**Coverage note.** Route C now covers all 64 dimensions (128 rows); the previous table covered 52
(104 rows). The 12 that were missing are the dimensions flagged `vanished`, and they are not
inert: the 24 new rows contribute 3 `convergent` verdicts and all 4 `both_routes_different_family`
ones - a category that literally could not appear before, since no row had both routes to compare.
`PRUNE_VANISHED = False` is doing real work.

### Result tables (`tables/<collection>/`)

One `.csv` per file, the caveat in the `#` header, the collection and the run id in the name:
`tables/scie/convergence_scie_drvi_epi_64.csv`, `tables/emt/convergence_emt_drvi_epi_64.csv`.
The `*` below stands for `_<collection>_<run_id>`.

| File | Contents |
|---|---|
| `signatures_<collection>.gmt` | the collection, provenance in the description field |
| `coverage_*.csv` | signature, axis, primary/robustness, provenance, n genes, n mapped, fraction, n in HVG background |
| `jaccard_*.csv`, `shared_genes_*.csv` | pairwise overlap of the collection |
| `signature_concentration_*.csv` | per signature: effective gene count, dominant gene, detectability, flags |
| `signature_gene_contribution_*.csv` | per (signature, gene): detection rate, share of the score's level and of its variance, ρ with the score |
| `cycle_confound_by_readout_*.csv` | per readout: cell-cycle genes it contains, the variance they carry, and the R² of the cycle on its score |
| `cycle_confound_by_vote_*.csv` | per vote threshold: AUROC of the cycle predicting membership, % G1, % G2M |
| `cycle_confound_by_dimension_*.csv` | per dimension: ρ with S and G2M, cycle loading, strongest readout association, convergent flag |
| `confounders_*.csv` | Spearman of each raw score vs depth, mito, S, G2M (ρ and p) |
| `confounder_checks_*.csv` | the named risks of the collection as scalars: G1, and depth (`scie`) or doublet (`emt`) |
| `quadrant_stability_*.csv` | Jaccard of the called cell set across the collection's definitions of the target |
| `quadrant_vote_distribution_*.csv` | how many definitions call each cell |
| `quadrant_per_patient_*.csv`, `quadrant_per_cell_type_*.csv` | group sizes, per patient and per label |
| `dim_signature_spearman_*.csv` | 64 dimensions x every readout, Route A |
| `dim_target_effect_size_*.csv` | AUROC and standardised mean difference per dimension |
| `dimension_row_order_*.csv` | the row order both routes share, with the `vanished` flag |
| `dim_geneset_signed_significance_*.csv` | 64 dimensions x the collection's gene sets, signed −log10 FDR, Route B |
| `factor_first_significant_*.csv`, `factor_first_hallmark_significant_*.csv` | every significant pair |
| `target_axes_*.csv` | axes carrying both halves of the target, when there are any |
| **`convergence_*.csv`** | **the main result: 128 rows, one per dimension-direction** |

### Figures (`figures/<step>/<collection>/`)

`04_3_signatures/`: `signature_coverage`, `jaccard_signature_overlap`, `signature_composition`.
`04_5_cell_first/`: `confounder_heatmap`, `quadrant_stability`, `dim_signature_heatmap`, and the
plane, whose name is the collection's: `stemness_immunogenicity_plane` for `scie`,
`emt_coexpression_plane` for `emt`.
`04_6_factor_first/`: `dim_geneset_signed_heatmap`.
`04_8_cycle_confound/`: `cycle_behind_stemness`.
`04_7_convergence/`: `routes_side_by_side`, `convergence_scatter`.
`04_9_embedding_control/`: `dim_signature_heatmap_side_by_side`, `association_concentration`,
`max_abs_rho`, `information_vs_alignment` - plus, outside the collection folders, **one folder
per method** (`harmony/`, `pca/`) holding that space's own UMAPs and its own Route A heatmap.

### Rules that keep the two routes independent

- **Signature genes are never added to the DRVI training feature set.** The HVG list of 04_1 is
  used exactly as 04_1 wrote it. Injecting the signatures and then reporting a dimension enriched
  for them would be circular and would destroy the evidential value of convergence.
- **No dimension is dropped before the analysis.** The vanished flags are read programmatically
  from `var['vanished']` / `var['vanished_*_direction']` - never from a plot - and reported
  alongside the results instead of filtering them.
- **The signatures are not independent** and the BH correction is not presented as that many
  independent tests. See the Jaccard matrix - on the EMT collection the block structure is
  extreme, list B vs list C reaching 0.76 on the mesenchymal axis.
- **The two collections are corrected separately.** 04_6 applies BH inside one collection, so a
  SCIE p-value cannot move because the EMT lists were added, and neither run may be read as
  having been corrected for the other.
- **Nothing is promoted on a single route.** All three categories are reported separately.

## 04_9_embedding_control (Route A on a space that is not DRVI's)

A **control**, not a result, and the one the README's own justification for DRVI invites:
Route B "is also what makes the choice of DRVI pay off (…) without it any of the
higher-scoring phase-02 methods would have done". That sentence is about Route B. It says
nothing about **Route A**, which is the route that assigns the cells - the deliverable. This
step asks Route A of **Harmony** and reports the two answers side by side.

Either outcome is worth writing down:

- if Harmony's corrected PCs carry the readouts as strongly as DRVI's dimensions, then Route
  A is **not** evidence for DRVI and the whole weight of the choice sits on Route B;
- if the association is **concentrated** on DRVI (one dimension carrying a readout) and
  **smeared** on Harmony (ten PCs each carrying a bit), that is the disentanglement claim
  measured on this dataset instead of cited from the DRVI paper.

### ⚠ What is under test here, and what is not

**DRVI is the hypothesis; Harmony is the reference level, not a rival.** The asymmetry is
deliberate and every output of this step has to be read through it: **Harmony has never
claimed axis-level interpretability.** It claims batch correction with the biological signal
preserved, and it was benchmarked on exactly that, against nine other methods, in phase 02.
A low concentration for Harmony is therefore **not a defect of Harmony**, and nothing in this
step may be quoted as ranking the two methods - phase 02 is where they are ranked, on what
they both promise. What is tested here is whether **DRVI delivers on its own promise** for
these readouts, with Harmony as the level a space that makes no such promise happens to
reach. The figures of this step carry that sentence as a subtitle for the same reason.

That splits the comparison into **two questions, and only the first is fair to both**:

| | question | how it is measured | fair to both? |
|---|---|---|---|
| 1 | does the space **carry** the state? | AUROC of the consensus target region from **all** the dimensions, L2 logistic regression, 5-fold CV **grouped by patient** | **yes** - it asks Harmony for nothing Harmony does not offer |
| 1b | is it on one axis, as a separation? | `auroc_best_single_dimension`, **oriented**: `0.5 + abs(AUROC - 0.5)`. A latent dimension has no privileged sign - the phase says so everywhere else, *the sign of rho IS the direction* - so an axis separating the target at 0.315 separates it exactly as well as one at 0.685. The signed value is kept in the column next to it | yes |
| 2 | is the state on **one axis**? | `max abs rho` per readout, and how concentrated the association is across dimensions | no - this is DRVI's own claim, and Harmony is the floor |

The **gap between them** is the quantity that makes the useful statement possible: it is the
state a space *represents but does not isolate*. "Harmony carries the same information but on
no single axis, DRVI puts it on one" is a defensible conclusion; neither number can reach it
alone. And if the multivariate AUROCs match **and** the concentrations match, Route A is not
an argument for DRVI at all - which is a result, not a failure of the control.

Grouped by patient rather than shuffled, because the question a multivariate readout has to
answer is whether the state is recoverable in a patient the fit has never seen; a random
split would put cells of one patient on both sides and report the patient effect as a
success. The classifier is deliberately **linear**: a non-linear one would measure its own
capacity as much as the geometry of the space, and the geometry is the question.

### Two fairness caveats, reported rather than corrected

**`effective_n_dims` cannot be read without the space's own redundancy next to it.** It mixes
two things - how aligned a readout is with a single axis, and how correlated the axes are
with each other, because correlated axes *share* an association. The two spaces differ
structurally on the second, and the difference is measured, not assumed:

| | dimensions | effective rank of the space | mean abs r between dimensions |
|---|---|---|---|
| DRVI (`drvi_epi_64`) | 64 | **30.7** | 0.094 (max 0.70) |
| Harmony (`harmony_epi_50`) | 50 | **43.6** | 0.037 |
| the uncorrected PCA | 50 | 50 exactly - orthogonal by construction | 0 |

DRVI's 64 dimensions behave as ~31 independent ones, Harmony's 50 as ~44 - the correction
breaks the PCA's exact orthogonality, but only slightly. Redundant axes spread any readout
over more of them, so **that bias runs against DRVI, not against Harmony.** Hence
`effective_rank_of_space` and `mean_abs_dim_correlation` are columns of the comparison table
and a subtitle of the concentration figure: the number never travels alone.

**The maximum of |ρ| over 64 draws is slightly inflated against the same maximum over 50.**
Small, and it cannot be removed without retraining DRVI at 50, which would change the run the
whole phase is built on. It is stated instead, and it is one more reason the counts are
reported as a percentage of the space rather than as a count.

**The coordinate system is an argument of Route A, not part of it.** A1 - A5 of
`cell_first_epi.py` - scoring, within-`(cohort, cell_type)` standardisation, the target
region, the consensus vote, the confounder and G1 and depth checks - are computed from
`shiao_epi.h5ad` and never see an embedding; only A6 does. So there is **no second copy of
Route A**: `cell_first_epi.py` takes `--embedding {drvi,harmony,pca}` exactly as it takes
`--collection`, and a control run writes only the three outputs that depend on the space -
`dim_signature_spearman`, `dim_target_effect_size`, `dimension_row_order` - plus its own
heatmap. Everything describing the cells stays owned by the `drvi_epi_64` run and is
reported `[skip]`. That is the point rather than an economy: **the cells being placed are
the same cells in both spaces**, so the only thing that varies is the axes.

**Harmony is re-run on the subset**, never reused from phase 02, for the reason 04_1 re-did
the pre-processing and 04_2 retrained DRVI: the phase-02 Harmony ran on 619,693 cells over
HVGs dominated by immune and stromal genes. It runs here on `shiao_epi_hvg_2k.h5ad`, **the
same 2,000 batch-aware HVGs DRVI was trained on**, so the two spaces differ by the method and
not by the genes.

**Route B and Route C are DRVI-only and are not parameterised.** They read the additive
decoder off `embed.varm`, which no other method has; the linear analogue would be the PCA
loadings, and Harmony corrects the embedding and *not* the loadings, so a "gene programme"
for a harmonised PC would be the programme of the PC **before** correction. The embeddings
written here carry no `varm` at all, so 04_6 fails on them immediately rather than quietly
producing something that reads like a Route B result. This is also why 04_9 has its own
driver instead of a flag on `signature_interpretation_all.sh`.

### The call, and what is pinned

`scib.integration.harmony`'s two lines, reproduced rather than imported the way
`02_2_integration/integration_methods.py` reproduces scVI's:

```python
sc.tl.pca(adata)                                             # 50 components
harmonize(adata.obsm["X_pca"], adata.obs, batch_key="cohort")
```

Everything that affects the result is kept identical to scib - the same **50** components
(scanpy's default, i.e. the size the phase-02 Harmony ran at), `batch_key='cohort'`, and
harmony-pytorch's own defaults for `theta`, `sigma`, `ridge_lambda` and the two iteration
caps. The two additions are `svd_solver='arpack'` and `random_state`, pinned so the run is
reproducible; scib leaves both to scanpy and a thesis table should not depend on that.

**50 against DRVI's 64 is not a handicap on any per-readout maximum** - the PCA is nested, so
the extra dimensions of a larger space can only add - but it *is* one on anything that
**counts** dimensions, which is why the comparison reports the effective count and the
percentage rather than the raw number.

**The uncorrected PCA comes free.** Harmony corrects the PCA it has just computed, so writing
that PCA out as well costs one extra UMAP and gives Route A a null arm: `--embedding pca`
separates "the integration matters for this readout" from "any 50-dimensional linear space
would have done". `run_harmony_epi.py` writes it by default (`--no-pca-arm` to skip);
`embedding_control_all.sh` asks for it only under `--embeddings "harmony pca"`, so the
default chain neither writes nor analyses an arm nobody requested.

### The three numbers, per readout

| Column | What it says |
|---|---|
| `max_abs_rho` | the strongest association any single dimension carries. This is the bar 04_7 applies (\|ρ\| ≥ 0.20), so it is what decides whether a readout has an axis at all |
| `pct_dims_abs_rho_ge_0.2` | how much of the space clears that bar, as a percentage - a percentage because 50 and 64 dimensions must not be compared by a count |
| **`effective_n_dims`** | inverse Simpson on the shares of Σρ², i.e. **the number of dimensions the association behaves as if it were spread over**. 1.0 is a readout living on one axis; 12 is a readout smeared across the space |
| `effective_rank_of_space`, `mean_abs_dim_correlation` | the redundancy of the space itself, repeated on every row so the column above is never read alone (see the caveats) |

The third is the disentanglement claim made measurable, and it is the same statistic
`signature_concentration` applies to the *genes inside* a score, applied here to the
*dimensions outside* it. Two spaces can reach the same `max_abs_rho` and still say very
different things about whether a dimension **is** a cell state: a strong association on a
space with an effective count of 10 is a programme the space represents but does not isolate.

### Execution order

| # | File | What it does | Where |
|---|------|--------------|-------|
| 1 | `run_harmony_epi.py` | PCA(50) + Harmony on `shiao_epi_hvg_2k.h5ad` → `embed_harmony_epi_50.h5ad` (and the uncorrected `embed_pca_epi_50.h5ad`), in DRVI's embedding format; the 04_2 figure set minus the decoder-only panels | local |
| 2 | `../04_5_cell_first/cell_first_epi.py --embedding harmony` | the same Route A, per collection | local |
| 3 | `compare_embeddings_epi.py` | the runs side by side: the comparison tables and the four figures. Re-derives the target region from the cached per-cell scores with `signature_common`'s own quadrant functions - not a second copy of the rule - and **checks it against the vote distribution 04_5 wrote**, stopping if they disagree: two spaces compared on two different cell sets is not a comparison. Give it the same `--high-q/--low-q/--mid-lo-q/--mid-hi-q` 04_5 was run with; `--no-multivariate` skips the CV fits, the only slow part | local |

```bash
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
cd 04_drvi_epithelial/04_9_embedding_control

./embedding_control_all.sh                            # harmony, both collections, resuming
./embedding_control_all.sh --collection scie          # one collection
./embedding_control_all.sh --embeddings "harmony pca" # with the uncorrected null arm
./embedding_control_all.sh --dry-run                  # print what would run
./embedding_control_all.sh cellfirst compare          # only the named steps
```

Step names: `harmony`, `cellfirst`, `compare`; a step whose output exists is `[have]` and
skipped, as in the other drivers. Environment `benchmark-py-r`, which already carries
`harmony-pytorch` 0.1.8 from phase 02 - **nothing new to install**. Everything is CPU and
local: the PCA is a couple of minutes on 74k x 2,000, `harmonize` a handful, and the UMAPs
the rest; no scheduler.

**Adding a third method** - scVI, scANVI, anything with an embedding - is one entry in the
`EMBEDDINGS` registry of `utils/signature_common.py` plus a writer that produces
`.X` = cells x dimensions with `var` carrying `title` / `order` / `vanished`. No step script
changes, exactly as adding a third collection does not change one.

### Outputs

    tables/<coll>/dim_signature_spearman_<coll>_harmony_epi_50.csv    # Route A on Harmony
    tables/<coll>/dim_target_effect_size_<coll>_harmony_epi_50.csv
    tables/<coll>/dimension_row_order_<coll>_harmony_epi_50.csv       # `harmony_order`, not `drvi_order`
    tables/<coll>/embedding_comparison_<coll>_epi_embeddings.csv      # one row per (space, readout)
    tables/<coll>/embedding_target_effect_<coll>_epi_embeddings.csv   # the quadrant on each space

    figures/04_9_embedding_control/
    ├── pca_variance_ratio_epi_50.png     the elbow of the PCA BOTH arms come from, hence
    │                                     at the root and not inside either method
    ├── harmony/                          one folder per METHOD, named for the method
    │   ├── umap_*_harmony_epi_50.png     its own UMAPs: 04_2's set, same panel names, so a
    │   ├── dims_in_heatmap_*             Harmony figure and its DRVI counterpart differ only
    │   ├── umap_dims_harmony_epi_50.png  by the run id in the filename
    │   └── <coll>/dim_signature_heatmap_<coll>_harmony_epi_50.png
    │                                     its Route A heatmap. NOT in figures/04_5_cell_first/,
    │                                     which is the phase's own run - a Harmony heatmap next
    │                                     to the DRVI one invites the two to be read as one
    │                                     result set - and not loose among the comparison
    │                                     figures either
    ├── pca/                              the same for the null arm
    └── <coll>/                           what is genuinely cross-method, and only that:
        dim_signature_heatmap_side_by_side, association_concentration, max_abs_rho,
        information_vs_alignment - question 1 against question 2, the one the control is for

A third method added to the registry gets its own folder and reorganises nothing.

The cross-run files carry `epi_embeddings` where a run id would be: they are about several
runs, so naming them after one would be a lie. Their `#` header falls back to the bare run id
for the same reason, while a control run's tables say `Harmony run harmony_epi_50` where the
phase's own say `DRVI run drvi_epi_64`.

### Result: the control does not support Route A as an argument for DRVI

Run on 2026-09-01, both collections, DRVI (`drvi_epi_64`, 64 dims) against Harmony
(`harmony_epi_50`, 50 dims), same 74,441 cells, same scores, same target region - re-derived
and checked against the vote distribution 04_5 wrote (4,145 cells on `scie`, 3,266 on `emt`).

**Question 1, does the space carry the state - Harmony, slightly.**

| | `scie` (4,145 cells) | `emt` (3,266 cells) |
|---|---|---|
| DRVI, whole space, patient-grouped CV | 0.850 ± 0.034 | 0.774 ± 0.025 |
| **Harmony**, same | **0.868 ± 0.029** | **0.815 ± 0.017** |
| DRVI, best single dimension (oriented) | **0.732** (DR 27) | 0.647 (DR 45) |
| Harmony, same | 0.685 (HD 12) | **0.676** (HD 4) |

Harmony carries the target region **as well as or better than** DRVI on both collections. On
`scie` DRVI isolates it better on one axis (0.732 vs 0.685); on `emt` it does not (0.647 vs
0.676). Both spaces show a large gap between the whole space and their best axis - DRVI
+0.118 and +0.127, Harmony +0.183 and +0.139 - so **neither space puts this state on a single
dimension**.

**Question 2, is it on one axis - not more so on DRVI.**

- `max abs rho`: DRVI is the stronger of the two on **8 of 11** readouts on `scie` and **9 of
  12** on `emt`, with a median paired difference of **-0.032** and **-0.025** in Harmony's
  disfavour. That is DRVI's one clear advantage, and it is small. Both spaces clear the 04_7
  bar on all 11 `scie` readouts; on `emt` Harmony clears it on 9 and DRVI on 7.
- **concentration, which is the disentanglement claim, goes the other way.** Median
  `effective_n_dims` is **16.8 on DRVI against 10.1 on Harmony** (`scie`) and **20.1 against
  13.3** (`emt`): the association is spread over *more* DRVI dimensions, not fewer. Correcting
  for the redundancy of each space widens the gap rather than closing it - 16.8 of 30.7
  effective DRVI directions (55%) against 10.1 of 43.6 Harmony ones (23%).

**What this licenses, and what it does not.** It does not say Harmony is the better method:
Harmony is here as the reference level and is ranked in phase 02, on what it promises. It says
that **for these readouts, on this compartment, Route A gives no evidence for DRVI's
axis-level interpretability** - a batch-corrected PCA reaches the same conclusions about the
same cells, and concentrates them on fewer axes. The justification for DRVI in this phase
therefore rests **entirely on Route B**, exactly as the README's own sentence said it did,
and that sentence should now be read as the whole of the argument rather than the main part
of it. Route C's convergent verdicts are unaffected: they are joint statements about a DRVI
dimension and its decoder, and nothing here touches either.

The honest way to present it: *"the interpretation is not an artefact of the embedding - it
reproduces on an independent, non-disentangled space - and the case for DRVI is the decoder,
not the axes."* The first half is a robustness result worth having on its own.

## Data location (`DATA_DIR`)

As in every phase, the heavy objects live **outside** the repo and every script resolves
`DATA_DIR` from the environment. This phase writes into its own subfolder, `04_epi/`, so the
epithelial chain is never confused with the non-immune one of 03:

    $DATA_DIR/
    ├── shiao.h5ad                             # input, from 01_5
    └── 04_epi/
        ├── shiao_epi_raw.h5ad                 # 1  (subset + refilter, raw counts)
        ├── shiao_epi_norm.h5ad                # 3  (scran log-norm)
        ├── shiao_epi_norm_cc.h5ad             # 4  (+ cell cycle)
        ├── shiao_epi_reduced.h5ad             # 5  (+ HVG/PCA/neighbours/UMAP)
        ├── shiao_epi_hvg_2k_list.csv          # 5  (selected HVG symbols)
        ├── shiao_epi_hvg_2k.h5ad              # 5  (DRVI input for 04_2)
        ├── shiao_epi.h5ad                     # 6  (+ leiden) - definitive epithelial object
        ├── shiao_epi_leiden_resolution_profile.csv   # 6  (resolution vs NMI)
        ├── embed_harmony_epi_50.h5ad          # 04_9 (Harmony on the same 2k HVGs, the control)
        ├── embed_pca_epi_50.h5ad              # 04_9 (that PCA uncorrected, the null arm)
        ├── model_drvi_epi_64.pt               # 04_2 (the trained DRVI model, one flat file)
        ├── embed_drvi_epi_64.h5ad             # 04_2 (latent space + stats + OOD/IND scores)
        ├── shiao_epi_drvi_epi_64.h5ad         # 04_2 (the 04_1 object + obsm['X_drvi'])
        #  the discarded 32 run keeps the same three names with _32
        │
        ├── signature_scores_<coll>_drvi_epi_64.csv    # 04_5 per cell, raw + within-stratum z scores
        ├── cytotrace2_drvi_epi_64.csv                 # 04_4 per cell, 74,441 rows (not collection-scoped)
        ├── factor_first_top200_genes_drvi_epi_64.tsv  # 04_6 ranked list per dimension-direction,
        │                                              #    read off the decoder, so shared by both collections
        ├── factor_first_top200_<coll>_drvi_epi_64.tsv # 04_6 every tested pair (the cache), per collection
        └── msigdb_hallmark_2020.json          # 04_6 cached library, so a re-run needs no network

`$DATA_DIR/signatures/` holds the lab's `.txt` files, the input to 04_3 (override the
location with `SIGNATURE_DIR`). The small result tables do **not** live here: they are versioned
in `04_drvi_epithelial/tables/`, because they are Appendix material rather than heavy objects.

## Object conventions (carried through every step)

- `.X` = **raw integer counts** in `shiao_epi_raw.h5ad`, **scran log1p-normalized** from step 3
  onward (already in log space - do **not** re-`log1p`). `.layers['counts']` = raw integer counts,
  kept unchanged throughout, and the layer DRVI trains on.
- `var_names` = gene symbols (Ensembl in `.var['gene_ids']`).
- `batch_key = 'cohort'` (patient), `label_key = 'cell_type'` (CellTypist, 10 labels observed).
- `compartment` is **constant** (`epi`) by construction, and so is `fraction` (`non_imm`, since
  every epithelial label is non-immune): neither is usable as a covariate. `compartment` exists
  precisely because `fraction` is constant in the 03 object too and could not tell the two apart.
- `dataset_origin` (the technical CD45 sort) is carried along untouched but plays no role here.

## Key parameters (verbatim, ready for Materials & Methods)

**Subsetting and re-filtering (step 1)**
- Subset: `cell_type` in the 11 epithelial labels of `Cells_Adult_Breast.pkl` - `LummHR-SCGB`,
  `LummHR-active`, `LummHR-major`, `Lumsec-HLA`, `Lumsec-KIT`, `Lumsec-basal`, `Lumsec-lac`,
  `Lumsec-major`, `Lumsec-myo`, `Lumsec-prol`, `basal`. Verbatim the "Epithelial" block of
  `01_4/fraction_reassignment.py`. 10 observed; `Lumsec-lac` has 0 cells in this dataset.
- `.X` reset to `.layers['counts']`; every derived slot and the emptied categories dropped.
- Cell filter (keep iff): `pct_counts_mt < 10` **and** `n_genes_by_counts > 100`, the thresholds of
  01_2 and 03_1 re-applied rather than assumed. Gene filter `sc.pp.filter_genes(min_cells=3)`, as
  in 03_1 and stricter than 01_2's `min_cells=1`.
- Scrublet is **not** re-run: doublet detection was done per biopsy on the full object, where the
  mixture of lineages is what makes a doublet detectable at all.
- Cohorts: `DROP_INCOMPLETE_COHORTS = True` **and** `DROP_SMALL_COHORTS = True` with
  `MIN_CELLS_PER_COHORT = 200`. **The one parameter that differs from 03_1**, where the size
  criterion was off. See below.

**Why the size threshold is on here (decided 2026-08-18)**

`cohort` is the batch key of everything downstream: scib's `hvg_batch` selects genes per patient
before merging, and DRVI fits a dispersion per gene *and* per batch. In 03 the smallest cohort
carried 314 cells, so the criterion never bound; on the epithelium the per-patient counts span
**14 to 9,243**, and without a threshold Patient06 would be a dispersion column fitted on 14
cells. `200` puts the floor back at 209 - the regime 03 worked in - for 1.3% of the cells:

| dropped | reason | cells |
|---|---|---|
| Patient01 | 0 cells in PD1 | 542 |
| Patient30 | 0 cells in RTPD1 | 87 |
| Patient06 | 14 cells | 14 |
| Patient23 | 176 cells | 176 |
| Patient04 | 192 cells | 192 |
| **total** | | **1,011 of 75,452 (1.3%)**, 34 → 29 cohorts |

**Normalization (step 3)** - identical to 01_3 and 03_1
- `scib.preprocessing.normalize`: `min_mean=0.1`, `log=True`, `precluster=True` (leiden
  `quickCluster`), `sparsify=False`; `SEED=0`.

**Cell cycle (step 4)** - identical to 01_4 and 03_1
- `sc.tl.score_genes_cell_cycle`, Tirosh/Regev 97-gene signature (43 S + 54 G2M), `random_state=0`.
- Re-scored rather than inherited: the reference gene sets are sampled from expression bins built
  on the current object. This matters more here than in 03: `Lumsec-prol` is a *proliferating*
  label and one of the two largest in the compartment, so the cycle bins of this population are
  not those of any wider one.

**Feature selection + reduction (step 5)** - identical to 01_5 and 03_1
- `scib.preprocessing.reduce_data`: `batch_key='cohort'`, `flavor='cell_ranger'`,
  `n_top_genes=2000`, `n_bins=20`, PCA 50 comps (`svd_solver='arpack'`), neighbours + UMAP on
  `X_pca`. HVG selection is per-patient then merged (`overwrite_hvg=True`).
- The script prints the overlap with **both** the 01_5 and the 03_1 HVG lists. The second is the
  number that justifies this phase over the previous one.

**Clustering (step 6)** - identical to 03_1
- `scib.clustering.cluster_optimal_resolution`, `label_key='cell_type'`, leiden
  (`flavor='igraph'`, `n_iterations=2`), selected by max NMI; every per-resolution column
  (`optscib_epi_leiden_<res>`) is kept. Resolution grid **0.1-2.0**, the full scib range: the 10
  observed labels are sub-labels of a single lineage, so the structure leiden has to resolve is
  finer than that count suggests.

**Signature interpretation (04_3 - 04_7)**

*Environment.* python 3.12.13, scanpy 1.12.1, anndata 0.13.2, numpy 2.4.6, pandas 3.0.3,
scipy 1.18.0, gseapy 1.3.0, statsmodels 0.14.6, scikit-learn 1.9.0, matplotlib 3.11.0,
seaborn 0.13.2, drvi-py 0.2.7. Conda env `benchmark-py-r`.

*Seeds.* `signature_common.SEED = 0` for the `score_genes` control sampling and the 20,000-cell
plot subsample; `CT2_SEED = 14` (the CytoTRACE2 default, pinned explicitly). The DRVI run itself
used `SEED = 123` in 04_2 and is not retrained.

*Latent dimensions.* `signature_common.PRUNE_VANISHED = False`: all 64 dimensions and all 128
dimension-directions enter both routes; `var['vanished']` is reported, not applied.

*04_5, Route A.* `sc.tl.score_genes(use_raw=False, ctrl_size=len(signature), n_bins=25,
random_state=0)` on the all-genes object. Standardisation `groupby(['cohort','cell_type'])`,
z-score, `ddof=0`. Quadrant `high_q=0.75`, `low_q=0.25`; consensus = majority of the stemness
definitions (≥ 3 of 6: `BENPORATH_ES1`, `ESC_ASSOU`, `ESC_WONG`, `LIM_STEM`, `FMASC`, CytoTRACE2).
Minimum signature size 10 mapped genes; minimum mapping fraction 0.60.

*04_6, Route B.* `OOD_combined` scores, `N_TOP_GENES=200`, both directions of every dimension.
`gp.enrich` (offline hypergeometric), background = the 2,000 HVGs of `shiao_epi_hvg_2k.h5ad`. BH
(`fdr_bh`) across all 7,680 dimension-direction / gene-set pairs at once, α = 0.05. Sanity-check
collection `MSigDB_Hallmark_2020`.

*04_7, Route C.* Route A bar `|ρ| ≥ 0.20`; Route B bar global FDR < 0.05; confounder flags at
`|ρ| ≥ 0.30` against depth and against `max(|ρ_S|, |ρ_G2M|)`.

*04_4, CytoTRACE2.* `species="human"` (**not** the package default `"mouse"`), `seed=14`,
`batch_size=10000`, `smooth_batch_size=1000`, `disable_plotting=True`, run per `cohort` from
`layers['counts']`. Input written as tab-delimited genes x cells `.txt`, the only format the
package reads.

## Reproducibility gaps (04_3 - 04_7)

Values that were not specified anywhere and had to be chosen. Every one is a module-level constant
or a CLI flag, so it can be changed and the chain re-run.

| # | Gap | Chosen | Where |
|---|---|---|---|
| 1 | Quadrant cutoffs - "stem-high / immunogenic-low" defines no threshold | q0.75 / q0.25 on within-stratum z | `cell_first_epi.py --high-q/--low-q` |
| 2 | How to combine the stemness definitions into one cell set | majority vote, `n_defs >= max(2, ceil(n_readouts/2))` - **≥ 3 of 6** (5 lists + CytoTRACE2). It was ≥ 4 of 7 while `EMP` was in the registry and ≥ 3 of 6 before CytoTRACE2; the bar is arithmetic, so a change in the number of readouts moves it | `cell_first_epi.py` |
| 3 | Route A significance bar - no effect size was specified | `|ρ| ≥ 0.20` | `convergence_epi.py --rho-min` |
| 4 | Confounder flag thresholds | `|ρ| ≥ 0.30` for depth and for cycle | `convergence_epi.py` |
| 5 | "Effect size (AUROC **or** standardised mean difference)" - requested as alternatives | both computed and reported | `dim_target_effect_size_*.csv` |
| 6 | `N_TOP_GENES` for Route B | 200, inherited from 03_3 rather than re-derived | `factor_first_epi.py --n-top-genes` |
| 7 | Cell-cycle conditioning method - "survives conditioning on phase" names no test | quadrant recomputed within G1 only, compared by Jaccard; the cycle is **not** regressed out | `cell_first_epi.py` |
| 8 | Whether to prune the vanished dimensions before the routes | not pruned: all 64 dimensions, all 128 directions, the flag reported instead | `signature_common.PRUNE_VANISHED` |
| 9 | Where small tables live - 03_3 wrote its tables to `$DATA_DIR` | versioned in `tables/`, per-cell matrices in `$DATA_DIR`; they are Appendix material | `signature_common.write_table` |
| 10 | Signature name ↔ file name for two files | `Immune consensus.txt` → `IMMUNOGENIC_CONSENSUS`, `fMaSC.txt` → `FMASC` | `signature_common.SIGNATURES` |
| 11 | `timepoint` `.obs` key | does not exist; `treatment` carries it. Nothing used it, so nothing was invented | above |
| 12 | CytoTRACE2 version | `cytotrace2-py` 1.1.0.4, pinned in its **own** env (`environments/cytotrace2-py.yml`) because it requires `numpy<2` and cannot coexist with `benchmark-py-r` | `environments/cytotrace2-py.yml` |

**Not a gap but a limitation to carry:** CytoTRACE2 has now been run on all 29 cohorts, so the
stemness axis no longer rests on the lab's five lists alone - but the sixth readout is not the
clean arbiter it was hoped to be. It is independent of the lists (which is what non-circularity
claims) and it clears the cycle check, `Lumsec-prol` ranking 5th of 10 cell types; it is **not**
independent of sequencing depth, rho +0.39 median even within `(cohort, cell_type)`, so it shares
the confounder that makes the five lists hard to read rather than adjudicating between them. The
stemness axis is therefore better evidenced than before and still weaker than the immune axis,
which has neither problem. See *04_4_cytotrace2 → Result*, and `confounders_*.csv` once 04_5 is
re-run with the `.csv` in place.

## Established numbers

Measured on `shiao.h5ad` before this phase runs:

- Epithelial cells: **75,452** of 619,693 (**12.2%**); 74,823 of them also survive the 03_1
  filters, i.e. the epithelium is **43%** of the non-immune compartment.
- **10** of the 11 epithelial labels are observed (`Lumsec-lac`: 0 cells). Three of them are 81%
  of the compartment: Lumsec-prol (24,647), LummHR-major (18,401), Lumsec-basal (17,729), then
  LummHR-SCGB (6,000) and basal (3,317); the remaining five are under 1,500 cells each and
  Lumsec-HLA has 97.

- All **34** cohorts carry epithelial cells; **29** survive the two cohort criteria, for
  **74,441** cells entering the gene filter. Smallest cohort after the drops: Patient02 (209),
  then Patient66 (214) and Patient20 (218); largest Patient52 (9,243).
- The cell filter removes **no cell**: every epithelial cell already passed `pct_mt < 10` and
  `n_genes > 100` in 01_2, and the gene set is still complete when it is re-applied.
- Genes after `min_cells=3`: **26,427** of 30,869, and **26,371** after the cohort drops (a whole
  patient leaving takes a few more genes under the threshold). The object step 3 receives is
  therefore **74,441 cells x 26,371 genes**.
- 04_2 at `n_latent = 64`: **12 dimensions vanish, 52 do not** - reported, and not used as a
  filter anywhere in 04_3 - 04_7.

Measured by reproducing step 1 headlessly on `shiao.h5ad`; the notebook prints the same numbers
when it runs. To be filled in from the chain: HVG overlap with 01_5 and with 03_1, selected leiden
resolution and cluster count.

## Critical methodological notes

- The subset is taken from the **unintegrated** `shiao.h5ad`, not from a phase-02 integration
  output and not from the 03 object. DRVI performs its own batch correction from raw counts, so an
  already-corrected matrix would be the wrong input, and re-using the 03 object would carry the
  QC/cohort decisions of that phase into this one.
- CellTypist is **not** re-run: `cell_type` is inherited from 01_4. Which means the epithelial
  definition is exactly as good as that annotation - a cell mislabelled `Fibro-major` by CellTypist
  is not in this object, and one mislabelled `LummHR-major` is. The leiden clusters of step 6 and
  the marker plots of step 7 are where such misassignments become visible.
- Malignant cells are **not** separated out. The atlas labels are normal-breast identities, so
  tumour epithelium is distributed across the luminal labels rather than carrying one of its own;
  no CNV inference (inferCNV / copyKAT) is run in this phase.
- Removing the stroma removes the largest source of variance left in the 03 compartment.
  Everything downstream of that - size factors, HVGs, PCA, neighbours, clusters - has to be
  recomputed for this reason and no other.
