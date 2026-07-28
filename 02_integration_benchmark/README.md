# 02_integration_benchmark

Third phase of the thesis: the **technical benchmark**. Ten batch-correction methods are run
on the unintegrated object produced by `01_pre_processing`, each output is scored with the
13 `scib` metrics, and the results are collected into a single summary table.

Input from phase 01: `shiao.h5ad` (619,693 cells x 30,869 genes) and
`shiao_hvg_2k_unintegrated_list.csv` (2,000 batch-aware HVGs).
Integration `batch_key = 'cohort'` (34 patients), biological `label_key = 'cell_type'`
(48 CellTypist labels).

## Repository layout

```
02_integration_benchmark/
├── README.md
├── benchmark_grid.tsv            # the run matrix: one row per run, drives everything
├── 02_1_prepare/                 # from shiao.h5ad to the benchmark inputs
├── 02_2_integration/             # the ten methods (local or SLURM)
├── 02_3_plot_method_umap/        # per-method UMAP QC panels (after 02_2, before metrics)
├── 02_4_metrics/                 # the 13 metrics (SLURM) + the final summary table
├── figures/                      # its output: per-method UMAPs + the summary table
└── utils/                        # smoke tests, shared helpers, vendored scIB plotting code
```

## Execution order

| # | Folder | Script | What it does | Where |
|---|--------|--------|--------------|-------|
| 0 | `utils` | `smoke_test_metrics.py` | Validates the metrics stack (scib, rpy2, R kBET, all 13 metrics) on a 5k-cell fixture before any real job is launched. | local + HPC |
| 0 | `utils` | `submit_smoke_test.slurm` | SLURM wrapper for the smoke test. | HPC (SLURM) |
| 0 | `utils` | `make_smoke_input.py` | Builds the tiny input for the *integration* smoke test (`smoke_hvg.h5ad` + `smoke_hvg_list.csv`) from `smoke_fixture.h5ad`. See "Smoke test of the integration step" below. | local |
| 0 | `utils` | `smoke_test_metrics_pipeline.sh` | Local smoke of the metrics *pipeline*: runs the 02_4 scripts over the 5k `smoke_out/` objects (one command) to validate the code path before any cluster job. Distinct from `smoke_test_metrics.py`, which validates the *stack*. See "Smoke test of the metrics step" below. | local |
| 1 | `02_1_prepare` | `subset_hvg.py` | Subsets `shiao.h5ad` to the 2,000 HVGs, asserts `layers['counts']` is still raw integers → `shiao_hvg_2k.h5ad`. | local |
| 2 | `02_1_prepare` | `scale_batch.py` | Per-batch z-scoring into a preallocated array indexed by row position, so cell order cannot change; PCA recomputed on the scaled matrix → `shiao_hvg_2k_scaled.h5ad`. | local |
| 3 | `02_1_prepare` | `h5ad_to_rds.R` | Converts each variant to a Seurat v3 `.rds` via `zellkonverter` (called once per variant), for the three R methods. | local |
| 3 | `02_1_prepare` | `hvg_csv_to_rds.R` | Converts the HVG symbol list to an `.rds` character vector, read by the Seurat anchor calls. | local |
| 4 | `02_2_integration` | `run_all.sh` | **The integration step.** Walks `benchmark_grid.tsv` and calls the right dispatcher for each row, **locally or on SLURM**. Default: runs each integration here, in sequence. `--slurm`: submits `submit_integration.slurm` for the matching rows, throttled to 3 concurrent tasks in partition `normal` (CPU). `--method` / `--scaling` / run_id filters, `--dry-run`. **Resumes by default**: a run whose `output` already exists is reported `[have]` and skipped (locally, and left out of the array spec on SLURM), so re-running the same command after a crash picks up where it stopped; `--force` overwrites. The R methods' `.rds` intermediate is deleted once converted, `--keep-rds` to keep it. | local + HPC |
| 4a | `02_2_integration` | `run_integration.py` + `integration_methods.py` | Python dispatcher and method bodies: BBKNN, Scanorama, Harmony, scVI, scANVI. Called per row by `run_all.sh`; runnable on its own. | local + HPC |
| 4b | `02_2_integration` | `run_integration.R` + `integration_methods.R` | R dispatcher and method bodies: fastMNN, Seurat CCA, Seurat RPCA. Called per row by `run_all.sh`; runnable on its own. | local + HPC |
| 4c | `02_2_integration` | `run_scgen.py` | scGen, in its own `scgen-py` environment (calls `SCGEN` directly; that env has no scib). On the cluster it needs an equivalent env, `$SCGEN_ENV` (see `submit_integration.slurm`). | local + HPC |
| 4d | `02_2_integration` | `run_drvi.ipynb` **(PLACEHOLDER — not written yet)** | DRVI, run interactively; must use the same cells and the same 2,000 HVGs as every other method. Not in the SLURM array (a notebook); run by hand. | local (GPU) |
| 4e | `02_2_integration` | `rds_to_h5ad.py` | Converts the R outputs back to `.h5ad`. Python via rpy2/anndata2ri, **not** `zellkonverter::writeH5AD`, which has no native R writer and would provision CPython through basilisk. | local + HPC |
| 4f | `02_2_integration` | `submit_integration.slurm` | The SLURM array task script `--slurm` submits: one array task per grid row, partition `normal` (CPU), throttled `%3`, no `--time`. Every method uses `catalano_env` except scGen (`$SCGEN_ENV`). | HPC (SLURM) |
| 5 | `02_3_plot_method_umap` | `plot_all.sh` | **The UMAP QC step.** Walks `benchmark_grid.tsv` and plots every integrated run into `figures/<run_id>/`, so one command covers the whole grid. Runs after 02_2 and **before** the metrics: a cheap visual check that a method mixed the cohorts without destroying the cell types (see "Figures" below). **Resumes by default**: a run whose five panels all exist is `[have]` and skipped; `--force` redraws. | local |
| 5a | `02_3_plot_method_umap` | `plot_methods_umaps.py` | Five UMAP QC panels for a single run. Called per row by `plot_all.sh`; runnable on its own. | local |
| 6 | `02_4_metrics` | `run_all_metrics.sh` | **The metrics step.** Walks `benchmark_grid.tsv`, expands the `types` column (17 runs → 21 jobs) and scores every (run, type), **locally or on SLURM**. Default: runs the check + metrics + kBET here, in sequence, then merges. `--slurm`: submits the two arrays (6f) restricted to the matching indices, throttled to 3 concurrent tasks, kBET `afterok`-dependent on the metrics array. `--method` / `--scaling` / run_id filters, `--no-kbet`, `--no-check`, `--dry-run`. **Resumes by default**, per (run, type): a task whose CSV already carries a kBET value is `[have]` and skipped, so a `full,embed` row half-scored resumes on the missing half; `--force` re-scores. The kBET row is the completeness marker because `metrics.py` writes the CSV and `metrics_kbet.py` fills that row afterwards. | local + HPC |
| 6a | `02_4_metrics` | `check_integrations.py` | Pre-flight check on one integrated object: cell count and order, batch/label keys, expected `obsm`/graph, finiteness. Cheap, and catches the failures that would otherwise surface hours into a metrics job. Called per (run, type); runnable on its own. | local + HPC |
| 6b | `02_4_metrics` | `metrics.py` | **The 13 metrics** minus kBET, one invocation per (method, scaling, output type). Called per (run, type); runnable on its own. | local + HPC |
| 6c | `02_4_metrics` | `metrics_kbet.py` | kBET alone, in a separate job/array; patches its value into the CSV `metrics.py` wrote. `--max-cells` optionally caps it (default: no cap, kBET runs on the full 620k in its own job). | local + HPC |
| 6d | `02_4_metrics` | `merge_metrics.py` | Merges the per-run CSVs into the single table consumed by the plotting code. Run once, after the jobs finish. | local + HPC |
| 6e | `02_4_metrics` | `make_summary_table.R` | Produces the final summary table (overall = 0.6 bio + 0.4 batch, min-max scaled). By default reproduces the published scIB figure with the vendored plotting code in `utils/` (`plotSingleTaskRNA.R` + `knit_table.R` + `img/`); falls back to a built-in scorer if its R packages (`dynutils`, `Hmisc`, `ggimage`) are missing. | local |
| 6f | `02_4_metrics` | `submit_metrics.slurm` + `submit_kbet.slurm` | The SLURM array task scripts `--slurm` submits: one array task per (run, type), in partition `normal`, throttled to 3 concurrent (`%3`) and with no `--time` (the node's default limit applies). Runnable directly with `sbatch` too. | HPC (SLURM) |
| - | `utils` | `metrics_shared.py` | Shared preparation for `metrics.py` and `metrics_kbet.py`: strict cell-order alignment, reference PCA/graph, the per-type `reduce_data` options and metric flags — so the two jobs score identically. | — |
| - | `utils` | `scib_compat.py` | Restores the numpy/pandas APIs scib 1.1.7 still expects. Imported by every script in this phase that touches scib. | — |
| - | `utils` | `h5ad_compat.py` | Writes `.h5ad` files that `anndata 0.10` can also read, and verifies the result. Used by steps 1 and 2. | — |

> **Note on step 0.** The smoke test is not optional bookkeeping. `scib 1.1.7` predates the
> numpy/pandas versions in the current environment and fails *halfway through* the metric
> computation, not at the start. 

> **Notes on steps 4, 5 and 6.** Each step has one wrapper that walks `benchmark_grid.tsv` and
> drives the whole grid with a single command; the lettered rows (`4a`–`4f`, `5a`, `6a`–`6f`) are
> the underlying per-run scripts, also runnable on their own for a single method or run. The
> integration (step 4, `run_all.sh`) and the metrics (step 6, `run_all_metrics.sh`) both run the
> grid **either locally or on SLURM** (`--slurm`): locally they loop in sequence; on the cluster
> they submit an array throttled to 3 concurrent tasks in partition `normal` (CPU), with no
> `--time` (the node's own limit applies) — for the metrics, the kBET array is `afterok`-dependent
> on the metrics one. The UMAP QC (step 5) is local only. On the cluster every method uses
> `catalano_env`, except scGen (its own `$SCGEN_ENV`) and the metrics stack; DRVI is a notebook
> and is always run by hand. The real runs are kept separate from the 5k smoke tests
> (`utils/smoke_test_metrics_pipeline.sh` for the metrics half).

## Smoke test of the integration step

There are **two** smoke tests in this phase, one per half of the pipeline. Step 0 above
(`smoke_test_metrics.py`) exercises the *metrics* stack. This one exercises the *integration*
stack: it runs every 02_2 dispatcher on a tiny object before the real grid is launched.

**Why.** The ten integration methods span three environments (`benchmark-py-r`, `scgen-py`,
and the R stack via `rpy2`), two languages, a GPU, and a Seurat↔AnnData round-trip. A single
run on the full 620k-cell object takes from minutes (BBKNN) to hours (scANVI, kBET-scale
Seurat), so a broken dispatcher, a missing dependency, a Seurat-version API change or a
silent cell reordering is expensive to discover late. The smoke test surfaces all of those in
a couple of minutes on 5,252 cells. It validates the **code path only** — *does the method
run, is the output the expected type (`knn` graph / `embed` / `full`), and is the input cell
order preserved* — **not** biological integration quality, which is meaningless at this size.
Cell order matters because the metrics compare objects cell by cell while checking names only
as a set, so a reordering passes validation and silently corrupts every score; each dispatcher
asserts the order survived, and the smoke test is where that assertion first runs for real.

**The fixture.** `smoke_fixture.h5ad` (5,252 cells, 4 cohorts, 7 cell types, `layers['counts']`
intact). `utils/make_smoke_input.py` restricts it to the HVGs it contains and writes
`smoke_hvg.h5ad` (through `h5ad_compat`, so `scgen-py`'s anndata 0.10 can read it too) plus
`smoke_hvg_list.csv`. The fixture holds 1,944 of the 2,000 HVGs; the Seurat methods must be
given exactly the HVGs that are present, because `FindIntegrationAnchors` fails with
`subscript out of bounds` on an anchor feature that is not in the object. On the real data all
2,000 are present, so this restriction is a property of the fixture, not of the benchmark.

**Usage.**

```bash
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
conda activate benchmark-py-r     # so rpy2 can find R

# 1. Build the small inputs (once)
python utils/make_smoke_input.py
Rscript 02_1_prepare/h5ad_to_rds.R   -i $DATA_DIR/smoke_hvg.h5ad -o $DATA_DIR/smoke_hvg.rds
Rscript 02_1_prepare/hvg_csv_to_rds.R -i $DATA_DIR/smoke_hvg_list.csv \
        -o $DATA_DIR/smoke_hvg_list.rds --n-expected 1944

O=$DATA_DIR/smoke_out            # throwaway output dir

# 2. Python methods
for m in bbknn harmony scvi scanvi scanorama; do
  python 02_2_integration/run_integration.py -m $m \
      -i $DATA_DIR/smoke_hvg.h5ad -o $O/$m.h5ad --emb-out $O/$m.npy
done

# 3. R methods, then the .h5ad conversion the metrics read
for m in fastmnn seurat_cca seurat_rpca; do
  Rscript 02_2_integration/run_integration.R -m $m \
      -i $DATA_DIR/smoke_hvg.rds -o $O/$m.rds -v $DATA_DIR/smoke_hvg_list.rds
done
python 02_2_integration/rds_to_h5ad.py -i $O/fastmnn.rds     -o $O/fastmnn.h5ad \
       --types full,embed --emb-out $O/fastmnn.npy
python 02_2_integration/rds_to_h5ad.py -i $O/seurat_cca.rds  -o $O/seurat_cca.h5ad  --types full
python 02_2_integration/rds_to_h5ad.py -i $O/seurat_rpca.rds -o $O/seurat_rpca.h5ad --types full

# 4. scGen, in its own environment
conda run -n scgen-py python 02_2_integration/run_scgen.py \
      -i $DATA_DIR/smoke_hvg.h5ad -o $O/scgen.h5ad
```

Each dispatcher prints an `[out]` line reporting the `obsm`/graph it produced and asserts the
cell order internally; a non-zero exit or a failed assertion is the signal to fix the script
before touching the grid. On the tested PC: T1000 (4 GB) run scVI and scANVI one at a time, and run scGen
with `CUDA_VISIBLE_DEVICES=""` if other GPU processes are resident. The smoke outputs in
`smoke_out/` are disposable and can be deleted once every method has come back green.

## Smoke test of the metrics step

The counterpart of the integration smoke test, for the *scoring* half. It is a **local**
runner in `utils/` — `smoke_test_metrics_pipeline.sh` — kept separate from the real,
cluster-only metrics step (the two SLURM arrays) so no smoke logic ever leaks into the real
run. It walks the same grid but scores the tiny `smoke_out/` objects against `smoke_hvg.h5ad`
instead of the real integrations, running the whole metrics code path — the pre-flight check,
both metric jobs, the CSV tree, the merge — in a few minutes on 5,252 cells before a single
cluster job is submitted.

**Why, on top of step 0.** `smoke_test_metrics.py` (step 0) validates that the *environment*
can compute every metric on an identity integration. This validates the *pipeline*: that
`metrics.py` / `metrics_kbet.py` read the real integration outputs (embed / full / knn, from
five different writers) correctly, that the per-type metric selection is right, that kBET
patches cleanly into the CSV `metrics.py` wrote, and that the tree merges into the shape the
plotting expects. It checks the **code path only** — every expected metric produces a finite
value — not biological quality, which is meaningless at this size.

```bash
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
conda activate benchmark-py-r          # so rpy2 can find R (kBET)
cd 02_integration_benchmark
utils/smoke_test_metrics_pipeline.sh --dry-run   # preview the 9 runs / 11 jobs
utils/smoke_test_metrics_pipeline.sh             # -> $DATA_DIR/smoke_metrics/ + smoke_metrics_merged.csv
utils/smoke_test_metrics_pipeline.sh --method harmony   # one method; also --no-kbet / --no-check
```

It scores the nine unscaled smoke outputs (DRVI is a notebook, so it has no smoke output),
which expand into 11 metrics jobs; every expected metric per type must come back finite
(13 for `full`, 12 for `embed`, 7 for the `knn` graph output). The outputs under
`smoke_metrics/` are disposable.

## The run grid (`benchmark_grid.tsv`)

Every run is one row. The file is the single source of truth: the local driver and the
SLURM arrays both read it, so the benchmark matrix stays a readable piece of data instead of
logic scattered across scripts.

Every run is one row (17 rows). A `types` value with two entries expands into two metrics
jobs downstream, which is how the 17 runs become 21 metrics jobs.

| Column | Meaning |
|---|---|
| `run_id` | `<method>_<scaling>`, e.g. `scanorama_scaled` - stable key for the output file and the job name |
| `method` | bbknn, scanorama, harmony, fastmnn, seurat_cca, seurat_rpca, scgen, scvi, scanvi, drvi |
| `language` | `python`, `R` or `notebook` - selects the dispatcher |
| `env` | conda environment (`benchmark-py-r`, `scgen-py`) |
| `scaling` | `unscaled` / `scaled` |
| `input` | the prepared object the method reads: `.h5ad` for Python, `.rds` for the R methods |
| `output` | integrated object under `$DATA_DIR/02_integration/`, named `<run_id>.h5ad` |
| `types` | scib output type(s): `knn`, `embed`, `full`, `full,embed` or `embed,full` |
| `reference` | the `-u` object for the metrics - **must match the scaling variant**. This is an *unintegrated .h5ad*, not the reference *batches* of the anchor-based Seurat methods: those are `run_integration.R`'s `--reference` (default `Patient53,Patient16,Patient43`, the three largest patients), overridable for a whole run with `SEURAT_REFERENCE=...` |
| `hvgs` | `--hvgs` value. Uniformly `0` here: every reference is already HVG-subset, so scib must not re-select HVGs (fatal on the scaled, negative-valued data). The `2000` case does not arise because the benchmark subsets to HVGs before integrating. |

The `method` values are the canonical names, also used by the R dispatcher
`run_integration.R` (`fastmnn` / `seurat_cca` / `seurat_rpca`); the `figures/` subfolders are
named per run as `<method>_<scaling>` (the `run_id`). So the grid needs
no name-translation layer.

### Methods

| Method | Language | Output type(s) | Scaled variant |
|---|---|---|---|
| BBKNN | Python | knn | yes |
| Scanorama | Python | full + embed | yes |
| Harmony | Python (`harmony-pytorch`) | embed | yes |
| fastMNN | R (`batchelor`) | embed + full | yes |
| Seurat CCA | R | full | yes |
| Seurat RPCA | R | full | yes |
| scGen | Python (`scgen-py`) | full | yes |
| scVI | Python | embed | no |
| scANVI | Python | embed | no |
| DRVI | Python (notebook) | embed | no |

**17 integration runs, 21 metrics jobs.** Scanorama and fastMNN write both a corrected matrix
and an embedding in a *single* execution, so they are scored twice from the same file by
changing `--type` only. scVI, scANVI and DRVI read `layers['counts']` and ignore `.X`
entirely: a scaled variant would be a duplicate run, so it is not in the grid.

Harmony runs through `scib.integration.harmony` (Python), not through R — the R `harmony`
package is not required anywhere in this phase.

## Metrics

Thirteen metrics, in the two groups used for the summary score.

| Batch correction (5) | Bio conservation (8) |
|---|---|
| PCR batch | NMI cluster/label |
| batch ASW (`ASW_label/batch`) | ARI cluster/label |
| graph iLISI | cell type ASW (`ASW_label`) |
| graph connectivity | isolated label F1 |
| kBET | isolated label ASW |
| | graph cLISI |
| | cell cycle conservation |
| | HVG conservation |

Trajectory conservation is present in the original paper, but **excluded** from this benchmark.

Not every metric exists for every output type, and this is structural rather than a defect:
`knn` output (BBKNN) supports only the graph-based metrics — silhouette, PCR, cell cycle and
HVG conservation are unavailable; `embed` output additionally loses HVG conservation. 
The summary table will therefore have empty cells by construction.

## Data location (`DATA_DIR`)

As in phases 00 and 01, the heavy objects live **outside** the repo and every script resolves
`DATA_DIR` from the environment. The repo holds code and lightweight figures only.

```
$DATA_DIR/
├── shiao.h5ad                             # from 01_5, unintegrated reference (all genes)
├── shiao_hvg_2k_unintegrated_list.csv          # from 01_5, the 2,000 HVG symbols
├── shiao_hvg_2k.h5ad                      # 02_1  (778 MB, sparse, counts preserved)
├── shiao_hvg_2k_scaled.h5ad               # 02_1  (1.1 GB, dense by construction)
├── shiao_hvg_2k.rds                       # 02_1  (2.9 GB, Seurat v3, data slot sparse)
├── shiao_hvg_2k_scaled.rds                # 02_1  (11 GB, Seurat v3, data slot dense)
├── shiao_hvg_2k_unintegrated_list.rds                  # 02_1  (the HVG symbols, for Seurat anchors)
├── 02_integration/<method>_<scaling>.h5ad # 02_2  integrated objects
│                                          #       the R methods' <run_id>.rds intermediate is
│                                          #       deleted after conversion (--keep-rds to keep it)
├── 02_embeddings/<method>_<scaling>.npy   # 02_2  latent spaces, kept for the figures
├── 02_metrics/shiao/metrics/<scaling>/hvg/<method>_<type>.csv   # 02_4  per-run metric CSVs
├── 02_metrics_merged.csv                  # 02_4  merged table (input to make_summary_table.R)
│
│   # smoke-test artifacts, all disposable (see the two smoke-test sections)
├── smoke_fixture.h5ad                     # utils (37 MB, 5,252 cells)
├── smoke_hvg.h5ad / smoke_hvg.rds         # utils + 02_1  the tiny benchmark inputs
├── smoke_hvg_list.csv / smoke_hvg_list.rds  # utils + 02_1  the 1,944 HVGs the fixture has
├── smoke_out/<method>.h5ad|.rds|.npy      # 02_2  smoke integrations
├── smoke_metrics/                         # 02_4  smoke metric CSVs + summary figure
└── smoke_metrics_merged.csv               # 02_4  merged smoke table
```

## Where things run

**Both integration and metrics can run locally or on SLURM** (`run_all.sh --slurm`,
`run_all_metrics.sh --slurm`). The usual split is integration **locally** — the local
environment is verified and has a GPU — and metrics **on the cluster**, with integrated objects
transferred as they are produced. The `--slurm` path exists for running integrations on the
cluster too: all methods on partition `normal` (CPU), so scVI/scANVI/scGen are slower there than
locally, and scGen needs its own `$SCGEN_ENV` (an equivalent of the local `scgen-py`); DRVI stays
a by-hand notebook. On the cluster every SLURM array is throttled to 3 concurrent tasks in
`normal`.

Embeddings are also exported on their own (~120 MB each — 619,693 cells x 50 dims, float32 —
against several GB for a full object),
which keeps the transfer cheap for the `embed` methods and gives the figures step something
small and durable to read after the integrated objects are cleaned up.

## Figures

### Per-method UMAP panels (`02_3_plot_method_umap`)

`02_3_plot_method_umap/plot_methods_umaps.py` writes `figures/<method>_<scaling>/` with five
panels per run, from the integrated object and the unintegrated reference:

- `cohort`, integrated (single panel)
- `cohort`, integrated vs unintegrated (two side by side)
- `cell_type`, integrated (single panel)
- `cell_type`, integrated vs unintegrated (two side by side)
- `cohort` + `cell_type`, both integrated (two side by side, no unintegrated)

These are a visual QC run **after** `run_all.sh` and **before** the metrics: a broken
integration (cohorts unmixed, or cell types destroyed) is obvious here in minutes, well before
the multi-hour metric jobs. The integrated embedding is derived from the output type
(`embed` → UMAP on `X_emb`; `full` → PCA on the corrected `.X` then UMAP; `knn` → UMAP on the
corrected graph); the unintegrated UMAP is the reference's `obsm['X_umap']` from 02_1, not
recomputed. `cell_type` colours are inherited from `uns['cell_type_colors']` and `cohort` from a 
fixed palette, so a category keeps its colour across every panel and method.

`plot_all.sh` is the driver (twin of `run_all.sh`): it walks `benchmark_grid.tsv` and
plots every run whose integrated `.h5ad` already exists, into `figures/<run_id>/`. It picks the
`--type` from the grid's `types` column (`embed` wins when a run has one, so scanorama/fastmnn
are shown on their corrected embedding), matches each run to its scaling reference, and skips
rows not yet integrated. It takes the same `run_id` and `--scaling` filters as
`run_all.sh`, plus `--dry-run` to preview.

```bash
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
cd 02_integration_benchmark/02_3_plot_method_umap
./plot_all.sh --dry-run          # preview what would be plotted
./plot_all.sh                    # every integrated run
./plot_all.sh --scaling unscaled # only the unscaled half
```

### Final summary table (`02_4_metrics`)

The final summary table is produced by the official scIB plotting code
(`utils/plotSingleTaskRNA.R` + `utils/knit_table.R` + `utils/img/`, vendored from
theislab/scib-reproducibility), driven by `make_summary_table.R -i <merged.csv> -o figures`.

Overall score = 0.6 x bio conservation + 0.4 x batch correction, on min-max scaled metrics.

## Environments

| Where | Environment | Notes |
|---|---|---|
| local | `benchmark-py-r` | everything except scGen; must be **activated**, otherwise rpy2 cannot find R |
| local | `scgen-py` | scGen only; old pinned stack that conflicts with the main one |
| HPC | `catalano_env` | metrics, and every integration except scGen when `run_all.sh --slurm` is used (different name for `environments/benchmark-hpc.yml`) |
| HPC | `$SCGEN_ENV` | scGen only on the cluster; an equivalent of the local `scgen-py` (`submit_integration.slurm` defaults to that name) |

One dependency is **not** declared in the `.yml` files and must be installed explicitly:

```bash
Rscript -e "remotes::install_github('theislab/kBET')"     # required by the kBET metric
```

The R packages the funky heatmap needs (`dynutils`, `Hmisc`, `ggimage` and the rest of the
plotting stack) **are** declared in `benchmark-py-r.yml` and `benchmark-hpc.yml`.

## Critical methodological notes

- **Every Python script in this phase imports `utils/scib_compat.py` before `scib`.** scib
  1.1.7 uses `np.in1d` (removed in numpy 2.4), `pd.value_counts` (removed in pandas 3.0) and
  the positional `Series[0]` fallback (removed in pandas 3.0). Without the shims the jobs die
  partway through the metric computation, not at import.
- **The `-u` reference must match the preprocessing variant** (scaled with scaled). Scoring a
  scaled output against an unscaled reference attributes part of the z-scoring effect to the
  integration method and makes the two columns non-comparable.
- **`--hvgs 0` whenever the reference is already restricted to the HVGs.** With reference and
  integrated object at the same gene count, scib would otherwise re-run HVG selection, which
  on scaled (negative) values fails.
- **Cell order is load-bearing.** The metrics compare the two objects cell by cell while only
  checking names as a set, so a reordering passes validation and silently produces wrong
  numbers. `scale_batch.py` therefore never splits the object: it fills a preallocated array
  by row position, so the order cannot change rather than being restored afterwards.
  `check_integrations.py` enforces the same property on the integration outputs.
- **The scaling is not `scib.preprocessing.scale_batch`.** That function stitches its 34
  per-batch pieces back together with `anndata.AnnData.concatenate`, removed in anndata 0.13,
  and raises before scaling anything; recovering it would mean pinning anndata below 0.11,
  which conflicts with scanpy 1.12. `scale_batch.py` reimplements the same operation —
  `sc.pp.scale` per batch, scanpy defaults, no `max_value` clipping — so the result is
  numerically what scIB would have produced. This is an implementation deviation from
  Luecken et al. 2022, not a methodological one.
- **The scaled reference gets a fresh PCA, and this is not cosmetic.** `scib.metrics.pcr()`
  reuses `obsm['X_pca']` together with `uns['pca']['variance']` whenever both are present,
  and it is called with `recompute_pca=False`. On the unscaled object the PCA inherited from
  01_5 is exactly right, because it was computed on these same 2,000 genes, and keeping it
  gives all 21 jobs an identical baseline for free. Carrying that same PCA into the scaled
  object would instead hand PCR batch and cell-cycle conservation a baseline computed on a
  different matrix, so it is recomputed there with the parameters scib itself would use
  (50 components, arpack). The neighbour graph and the UMAP are dropped from the scaled
  object rather than recomputed: no metric reads them on the reference.
- **The benchmark inputs are written so that `anndata 0.10` can read them too.** scGen runs
  in `scgen-py`, pinned to anndata 0.10.8, while the inputs are written by anndata 0.13 under
  pandas 3, which emits three encodings older readers reject — `nullable-string-array` (used
  for the *index*, i.e. barcodes and gene symbols), `nullable-boolean`, and `null` for `None`
  values inside `.uns`. `utils/h5ad_compat.py` downcasts them before writing and re-opens the
  file afterwards to verify, so no separate legacy export is needed and there is only ever one
  copy of each input.
- **Seurat runs on the legacy anchor API** (`FindIntegrationAnchors` + `IntegrateData`), not
  on v5 `IntegrateLayers`. The v5 API returns an embedding rather than a corrected matrix,
  which would turn Seurat from a `full` method into an `embed` one and change which metrics
  can be computed. CCA and RPCA are both **reference-based** (largest patients as reference):
  with 34 batches and 620k cells, full pairwise anchoring is not tractable. Both points are
  deviations from Luecken et al. 2022 and are declared as such in Materials & Methods.
- **kBET is scored in its own job.** `metrics.py` never computes it (`kBET_=False` for every
  output type); it is by far the slowest metric, so isolating it in `metrics_kbet.py` means a
  timeout costs one metric instead of thirteen. `metrics_kbet.py --max-cells` can cap it, but
  the default is `0`: no cap, kBET runs on the full 620k cells in its own job.
- **The leiden resolution sweep is capped at 1.0**, as in `01_5_scib_pp/scib_clustering.py`.
  NMI and ARI are computed on an optimal-resolution clustering, recomputed for each of the 21
  metrics jobs: it is the single most expensive part of the phase.
