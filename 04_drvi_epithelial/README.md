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
question "is this a tumour cell state?" cannot be answered from these data. The same caveat
carries over unchanged into 05.

## Repository layout

```
04_drvi_epithelial/
├── README.md
├── signature_interpretation_all.sh   # phase-level driver for 04_3 -> 04_7
├── utils/
│   └── signature_common.py           # paths, signature registry, writers, DRVI accessor
├── 04_1_subsetting/                  # from shiao.h5ad to the epithelial object
├── 04_2_drvi_run/                    # DRVI on that object, n_latent 64
├── 04_3_signatures/                  # the lab's lists -> one .gmt, coverage, Jaccard
├── 04_4_cytotrace2/                  # per-patient potency (the one HPC step)
├── 04_5_cell_first/                  # Route A
├── 04_6_factor_first/                # Route B
├── 04_7_convergence/                 # Route C, the main result
├── tables/                           # every result table of 04_3 - 04_7
└── figures/                          # one folder per step: 04_1_*, 04_2_<run_id>, 04_3_*, 04_5_* ...
```

`utils/` follows 00 and 02: helpers shared by several steps live there, imported with the idiom
`02_2_integration/run_integration.py` uses. `tables/` is phase-level for the same reason
`figures/` is - 04_7 reads what 04_5 and 04_6 wrote, so a per-step `tables/` would mean steps
reaching into each other's folders.

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
| question | do the cells prior knowledge calls stem-like or immune-evasive sit anywhere in particular along a latent dimension? | what program does this dimension encode, irrespective of what I was looking for? |
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
| 1 | `04_3_signatures/build_signatures_epi.py` | Ingest the 11 text files; write `lab_signatures.gmt` with the provenance in the description field; coverage table; pairwise Jaccard. | local |
| 2 | `04_4_cytotrace2/cytotrace2_epi.py` | CytoTRACE2 potency, one patient at a time, from raw counts; concatenate. **Not in the default chain.** | HPC |
| 3 | `04_5_cell_first/cell_first_epi.py` | Route A: scoring, confounder table, within-patient standardisation, quadrant definitions, dimension x signature association. | local |
| 4 | `04_6_factor_first/factor_first_epi.py` | Route B: top genes per dimension **and direction**, offline ORA, signed significance matrix. | local |
| 5 | `04_7_convergence/convergence_epi.py` | Route C: the convergence table and the comparative figures. | local |

```bash
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
cd 04_drvi_epithelial

./signature_interpretation_all.sh                  # steps 1, 3, 4, 5 - resuming
./signature_interpretation_all.sh --force          # re-run everything
./signature_interpretation_all.sh --dry-run        # print what would run
./signature_interpretation_all.sh cellfirst convergence   # only the named steps
./signature_interpretation_all.sh cytotrace        # the HPC step, locally
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

Eleven files, one gene symbol per line, in `$DATA_DIR/signatures/` (override with
`SIGNATURE_DIR`). Ingested into `tables/lab_signatures.gmt`, whose **description field carries the
provenance**, so the collection that feeds both routes is also the Appendix table. The `.gmt`
holds the **mapped** genes, i.e. the collection as actually used, so the Appendix and the tested
sets cannot drift apart.

Five of the eleven files are **CRLF** and at least one carries a **UTF-8 BOM**; the reader opens
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
| `EMP` | stemness | Embryonic Multipotent Progenitors, PMID 29784918 |
| `LIM_STEM` | stemness | MSigDB `LIM_MAMMARY_STEM_CELL_UP` |
| `FMASC` | stemness | fetal mammary stem cells, PMC3277444, Suppl. Table 2 |

Two file names differ from the signature name: `Immune consensus.txt` → `IMMUNOGENIC_CONSENSUS`,
`fMaSC.txt` → `FMASC`.

How much of each list is actually measured is in `tables/coverage_*.csv`, against two universes:
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

**This is the only part of the chain with a SLURM wrapper**, and the reason it is its own step:
the dependency is not in `benchmark-py-r`, and the runtime is minutes per patient across 29
patients rather than seconds.

```bash
cd 04_drvi_epithelial/04_4_cytotrace2 && mkdir -p logs
sbatch --export=ALL,DATA_DIR=$DATA_DIR submit_cytotrace2_epi.slurm --cores 8
```

> **Not yet run.** `cytotrace2-py` (1.1.0.4 on PyPI, `pip install cytotrace2-py`) is **not** in
> `benchmark-py-r` and is not pinned in `environments/benchmark-py-r.yml`. The script checks for
> it and stops with that instruction rather than failing halfway; `--dry-run` exports the
> per-patient matrices without it. Re-run 04_5 afterwards and it picks the `.csv` up automatically,
> adding CytoTRACE2 as a seventh stemness readout. Without it Route A runs on the six signature
> lists alone **and says so in its output**.

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
comparison meaningless. The quadrant is defined **once per stemness readout** (six signatures,
plus CytoTRACE2 when present) and the consensus cell set is the majority vote across them; the
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
  directions x 61 gene sets = **7,808 pairs**, the pairs with no overlap included: they are p = 1
  and cannot become significant, but they belong in the denominator. With 100+ directions a
  per-query FDR is far too permissive;
- **no direction is pruned.** 03_3 let DRVI's accessor drop the directions it had marked vanished;
  here all 2 x 64 are tested (see *Vanished dimensions are not pruned* above);
- the output is a matrix on the **same row order as Route A**, read from
  `tables/dimension_row_order_*.csv`, so the two heatmaps are directly comparable.

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

### Result tables (`tables/`)

One `.csv` per file, the caveat in the `#` header, the run id in the name.

| File | Contents |
|---|---|
| `lab_signatures.gmt` | the collection, provenance in the description field |
| `coverage_*.csv` | signature, axis, provenance, n genes, n mapped, fraction, n in HVG background |
| `jaccard_*.csv`, `shared_genes_*.csv` | pairwise overlap of the collection |
| `confounders_*.csv` | Spearman of each raw score vs depth, mito, S, G2M (ρ and p) |
| `confounder_checks_*.csv` | the G1 and depth checks as scalars |
| `quadrant_stability_*.csv` | Jaccard of the called cell set across stemness definitions |
| `quadrant_vote_distribution_*.csv` | how many definitions call each cell |
| `quadrant_per_patient_*.csv`, `quadrant_per_cell_type_*.csv` | group sizes, per patient and per label |
| `dim_signature_spearman_*.csv` | 64 dimensions x every readout, Route A |
| `dim_target_effect_size_*.csv` | AUROC and standardised mean difference per dimension |
| `dimension_row_order_*.csv` | the row order both routes share, with the `vanished` flag |
| `dim_geneset_signed_significance_*.csv` | 64 x 11 signed −log10 FDR, Route B |
| `factor_first_significant_*.csv`, `factor_first_hallmark_significant_*.csv` | every significant pair |
| `target_axes_*.csv` | axes carrying both halves of the target, when there are any |
| **`convergence_*.csv`** | **the main result: 128 rows, one per dimension-direction** |

### Figures

`figures/04_3_signatures/`: `signature_coverage`, `jaccard_signature_overlap`.
`figures/04_5_cell_first/`: `confounder_heatmap`, `stemness_immunogenicity_plane`,
`quadrant_stability`, `dim_signature_heatmap`.
`figures/04_6_factor_first/`: `dim_geneset_signed_heatmap`.
`figures/04_7_convergence/`: `routes_side_by_side`, `convergence_scatter`.

### Rules that keep the two routes independent

- **Signature genes are never added to the DRVI training feature set.** The HVG list of 04_1 is
  used exactly as 04_1 wrote it. Injecting the signatures and then reporting a dimension enriched
  for them would be circular and would destroy the evidential value of convergence.
- **No dimension is dropped before the analysis.** The vanished flags are read programmatically
  from `var['vanished']` / `var['vanished_*_direction']` - never from a plot - and reported
  alongside the results instead of filtering them.
- **The signatures are not independent** and the BH correction is not presented as eleven
  independent tests. See the Jaccard matrix.
- **Nothing is promoted on a single route.** All three categories are reported separately.

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
        ├── model_drvi_epi_64.pt               # 04_2 (the trained DRVI model, one flat file)
        ├── embed_drvi_epi_64.h5ad             # 04_2 (latent space + stats + OOD/IND scores)
        ├── shiao_epi_drvi_epi_64.h5ad         # 04_2 (the 04_1 object + obsm['X_drvi'])
        #  the discarded 32 run keeps the same three names with _32
        │
        ├── signature_scores_drvi_epi_64.csv   # 04_5 per cell, raw + within-stratum z scores
        ├── cytotrace2_drvi_epi_64.csv         # 04_4 per cell (once that step has run)
        ├── factor_first_top200_genes_drvi_epi_64.tsv  # 04_6 ranked list per dimension-direction
        ├── factor_first_top200_drvi_epi_64.tsv        # 04_6 every tested pair (the cache)
        └── msigdb_hallmark_2020.json          # 04_6 cached library, so a re-run needs no network

`$DATA_DIR/signatures/` holds the lab's eleven `.txt` files, the input to 04_3 (override the
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
definitions (≥ 3 of 6). Minimum signature size 10 mapped genes; minimum mapping fraction 0.60.

*04_6, Route B.* `OOD_combined` scores, `N_TOP_GENES=200`, both directions of every dimension.
`gp.enrich` (offline hypergeometric), background = the 2,000 HVGs of `shiao_epi_hvg_2k.h5ad`. BH
(`fdr_bh`) across all 7,808 dimension-direction / gene-set pairs at once, α = 0.05. Sanity-check
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
| 2 | How to combine the six stemness definitions into one cell set | majority vote (≥ 3 of 6) | `cell_first_epi.py` |
| 3 | Route A significance bar - no effect size was specified | `|ρ| ≥ 0.20` | `convergence_epi.py --rho-min` |
| 4 | Confounder flag thresholds | `|ρ| ≥ 0.30` for depth and for cycle | `convergence_epi.py` |
| 5 | "Effect size (AUROC **or** standardised mean difference)" - requested as alternatives | both computed and reported | `dim_target_effect_size_*.csv` |
| 6 | `N_TOP_GENES` for Route B | 200, inherited from 03_3 rather than re-derived | `factor_first_epi.py --n-top-genes` |
| 7 | Cell-cycle conditioning method - "survives conditioning on phase" names no test | quadrant recomputed within G1 only, compared by Jaccard; the cycle is **not** regressed out | `cell_first_epi.py` |
| 8 | Whether to prune the vanished dimensions before the routes | not pruned: all 64 dimensions, all 128 directions, the flag reported instead | `signature_common.PRUNE_VANISHED` |
| 9 | Where small tables live - 03_3 wrote its tables to `$DATA_DIR` | versioned in `tables/`, per-cell matrices in `$DATA_DIR`; they are Appendix material | `signature_common.write_table` |
| 10 | Signature name ↔ file name for two files | `Immune consensus.txt` → `IMMUNOGENIC_CONSENSUS`, `fMaSC.txt` → `FMASC` | `signature_common.SIGNATURES` |
| 11 | `timepoint` `.obs` key | does not exist; `treatment` carries it. Nothing used it, so nothing was invented | above |
| 12 | CytoTRACE2 version | not installed; `cytotrace2-py` 1.1.0.4 is current on PyPI and **not** pinned in `benchmark-py-r.yml` | `environments/` |

**Not a gap but a limitation to carry:** CytoTRACE2 has **not been run**. Until it is, the
stemness axis rests entirely on the lab's six lists, which are heavily depth- and cycle-coupled
and disagree with each other - both measured in `confounders_*.csv` and `quadrant_stability_*.csv`.
There is currently **no non-circular evidence on the stemness axis**, and this is the single
largest weakness of the stage. The immune axis does not have this problem.

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
