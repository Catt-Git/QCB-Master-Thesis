# 02_integration_benchmark

Third phase of the thesis (the big one): the **technical benchmark**. Ten batch-correction methods are run
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
├── figures/                      # one folder per run: UMAP QC panels (+ DRVI latent-space
│                                 # figures for the DRVI runs), plus the summary table
└── utils/                        # smoke tests, cluster upload, shared helpers, vendored scIB
                                  # plotting code
```

## Execution order

One line per script; the details live in the sections linked from the last column.

| # | Folder | Script | What it does | Where | More |
|---|--------|--------|--------------|-------|------|
| 0 | `utils` | `smoke_test_metrics.py` | Validates the metrics *stack* (scib, rpy2, R kBET, all 13 metrics) on a 5k fixture. | local + HPC | [smoke, metrics](#smoke-test-of-the-metrics-step) |
| 0 | `utils` | `submit_smoke_test.slurm` | SLURM wrapper for the stack smoke test. | HPC | |
| 0 | `utils` | `make_smoke_input.py` | Builds the tiny inputs for the *integration* smoke test. | local | [smoke, integration](#smoke-test-of-the-integration-step) |
| 0 | `utils` | `smoke_test_metrics_pipeline.sh` | Runs the whole 02_4 *pipeline* over the 5k `smoke_out/` objects. | local | [smoke, metrics](#smoke-test-of-the-metrics-step) |
| 1 | `02_1_prepare` | `subset_hvg.py` | Subsets `shiao.h5ad` to the 2,000 HVGs → `shiao_hvg_2k.h5ad`. | local | |
| 2 | `02_1_prepare` | `scale_batch.py` | Per-batch z-scoring, then a fresh PCA **and UMAP** on the scaled matrix → `shiao_hvg_2k_scaled.h5ad`. | local | [notes](#critical-methodological-notes) |
| 3 | `02_1_prepare` | `h5ad_to_rds.R` | Converts each variant to a Seurat v3 `.rds` (`zellkonverter`), for the R methods. | local | |
| 3 | `02_1_prepare` | `hvg_csv_to_rds.R` | Converts the HVG symbol list to an `.rds` character vector, for the Seurat anchors. | local | |
| 4 | `02_2_integration` | `run_all.sh` | **Driver of the integration step**: the whole grid, locally or on SLURM. | local + HPC | [drivers](#the-three-drivers) |
| 4a | `02_2_integration` | `run_integration.py` + `integration_methods.py` | Python dispatcher and method bodies: BBKNN, Scanorama, Harmony, scVI, scANVI. | local + HPC | [models](#trained-models) |
| 4b | `02_2_integration` | `run_integration.R` + `integration_methods.R` | R dispatcher and method bodies: fastMNN, Seurat CCA, Seurat RPCA. | local + HPC | [notes](#critical-methodological-notes) |
| 4c | `02_2_integration` | `run_scgen.py` | scGen, in its own `scgen-py` / `$SCGEN_ENV` environment. | local + HPC | [models](#trained-models) |
| 4d | `02_2_integration` | `shiao_drvi_128.ipynb` | DRVI, run interactively: the latent size is chosen by eye. | local (GPU) | [latent size](#choosing-the-latent-size-drvi) |
| 4e | `02_2_integration` | `run_drvi.py` | The same DRVI run headless, once the size is chosen: same parameters, same outputs, same figures. | local + HPC | [latent size](#choosing-the-latent-size-drvi) |
| 4f | `02_2_integration` | `submit_drvi.slurm` | SLURM wrapper for `run_drvi.py`: one job, partition `long`, CPU. | HPC | [latent size](#choosing-the-latent-size-drvi) |
| 4g | `02_2_integration` | `rds_to_h5ad.py` | Converts the R outputs back to `.h5ad` for the metrics. | local + HPC | [notes](#critical-methodological-notes) |
| 4h | `02_2_integration` | `submit_integration.slurm` | The SLURM array task script for step 4. | HPC | [drivers](#the-three-drivers) |
| 5 | `02_3_plot_method_umap` | `plot_all.sh` | **Driver of the UMAP QC step**: five panels per integrated run. | local | [drivers](#the-three-drivers) |
| 5a | `02_3_plot_method_umap` | `plot_methods_umaps.py` | The five UMAP QC panels for a single run; caches the integrated layout into the object. | local | [figures](#per-method-umap-panels-02_3_plot_method_umap) |
| 5b | `utils` | `sync_to_cluster.sh` | **The local → cluster bridge**: uploads the files the grid names. | local | [cluster](#running-the-metrics-on-the-cluster) |
| 6 | `02_4_metrics` | `run_all_metrics.sh` | **Driver of the metrics step**: every (run, type), locally or on SLURM. | local + HPC | [drivers](#the-three-drivers) |
| 6a | `02_4_metrics` | `check_integrations.py` | Pre-flight check on one integrated object: cells, order, keys, `obsm`/graph, finiteness. | local + HPC | |
| 6b | `02_4_metrics` | `metrics.py` | **The 13 metrics** minus kBET, one invocation per (run, type). | local + HPC | [metrics](#metrics) |
| 6c | `02_4_metrics` | `metrics_kbet.py` | kBET alone, patched into the CSV `metrics.py` wrote. | local + HPC | [notes](#critical-methodological-notes) |
| 6d | `02_4_metrics` | `merge_metrics.py` | Merges the per-run CSVs into the single table the plotting reads. | local + HPC | |
| 6e | `02_4_metrics` | `make_summary_table.R` | The final summary table (scIB funky heatmap). | local | [summary table](#final-summary-table-02_4_metrics) |
| 6f | `02_4_metrics` | `submit_metrics.slurm` + `submit_kbet.slurm` | The two SLURM array task scripts for step 6. | HPC | [drivers](#the-three-drivers) |
| - | `utils` | `metrics_shared.py` | Shared preparation, so `metrics.py` and `metrics_kbet.py` score identically. | imported | |
| - | `02_2_integration` | `model_paths.py` | The single definition of the checkpoint layout. | imported | [models](#trained-models) |
| - | `utils` | `scib_compat.py` | Restores the numpy/pandas APIs scib 1.1.7 expects. | imported | [notes](#critical-methodological-notes) |
| - | `utils` | `h5ad_compat.py` | Writes `.h5ad` files `anndata 0.10` can read too. | imported | [notes](#critical-methodological-notes) |

A few of these deserve a line of their own:

- `run_integration.py` **reimplements scVI and scANVI** instead of calling `scib.integration`,
  which gives no access to the fitted model - same call sequence and same parametrisation, only
  the checkpoint added.
- `metrics_shared.py` holds the strict cell-order alignment, the reference PCA/graph and the
  per-type `reduce_data` options and metric flags.
- `model_paths.py` is stdlib only, so `scgen-py` can import it too.
- `h5ad_compat.py` is used by steps 1 and 2; `scib_compat.py` by every script here that touches scib.

## The three drivers

Steps 4, 5 and 6 each have **one wrapper that walks `benchmark_grid.tsv`** and drives the whole
grid with a single command; the lettered rows above are the per-run scripts underneath, also
runnable on their own for a single method or run.

The three share the same interface:

- **Filters** - `--method`, `--scaling`, or explicit run ids; `--dry-run` previews without doing
  anything. Unknown options, run ids and method names are rejected up front, so a typo cannot
  silently select nothing and exit 0.
- **Resume by default** - work already on disk is reported `[have]` and skipped, so re-running
  the same command after a crash picks up where it stopped; `--force` redoes it.

| Driver | Step | "Already done" means | Local | `--slurm` |
|---|---|---|---|---|
| `02_2_integration/run_all.sh` | integration | the row's `output` exists | each integration in sequence | `submit_integration.slurm`, one task per row |
| `02_3_plot_method_umap/plot_all.sh` | UMAP QC | the run's five panels exist | in sequence | - (local only) |
| `02_4_metrics/run_all_metrics.sh` | metrics | the (run, type) CSV carries a kBET value | check + metrics + kBET in sequence, then merge | `submit_metrics.slurm` + `submit_kbet.slurm`, one task per (run, type) |

Per-driver extras:

- `run_all.sh` deletes the R methods' `.rds` intermediate once converted; `--keep-rds` keeps it.
- `run_all_metrics.sh` expands the `types` column (17 runs → 21 jobs) and takes `--no-kbet` /
  `--no-check`. It resumes **per (run, type)**, so a half-scored `full,embed` row resumes on the
  missing half; the kBET row is the completeness marker because `metrics.py` writes the CSV and
  `metrics_kbet.py` fills that row afterwards. The merge and the summary table are run by hand
  once the arrays finish - the driver prints the two commands.
- `plot_all.sh` picks the `--type` from the grid's `types` column (`embed` wins when a run has
  one, so scanorama/fastmnn are shown on their corrected embedding), matches each run to its
  scaling reference, and skips rows not yet integrated. Both references must carry
  `obsm['X_umap']` (01_6 for the unscaled one, `scale_batch.py` for the scaled one) or the run
  aborts on the plotter's assertion.

**On SLURM.** Both arrays run in partition `normal` (CPU), throttled to 3 concurrent tasks (`%3`),
with no `--time` (the node's own limit applies); the kBET array is `afterok`-dependent on the
metrics one. Only the matching rows enter the array spec, so an incomplete grid (the scaled half
still integrating) costs nothing. Every method uses `catalano_env` except scGen
(`$SCGEN_ENV`); the `.slurm` scripts are runnable directly with `sbatch` too. DRVI is never in the
array - the latent size is a decision, not a grid row - but once that size is chosen the run
itself can go to the cluster as a single job (`submit_drvi.slurm`, partition `long`).

## Smoke test of the integration step

There are **two** smoke tests in this phase, one per half of the pipeline. This one exercises the
*integration* stack: it runs every 02_2 dispatcher on a tiny object before the real grid is
launched.

**Why.** The ten integration methods span three environments (`benchmark-py-r`, `scgen-py`,
and the R stack via `rpy2`), two languages, a GPU, and a Seurat↔AnnData round-trip. A single
run on the full 620k-cell object takes from minutes (BBKNN) to hours (scANVI, kBET-scale
Seurat), so a broken dispatcher, a missing dependency, a Seurat-version API change or a
silent cell reordering is expensive to discover late. The smoke test surfaces all of those in
a couple of minutes on 5,252 cells. It validates the **code path only** - *does the method
run, is the output the expected type (`knn` graph / `embed` / `full`), and is the input cell
order preserved* - **not** biological integration quality, which is meaningless at this size.
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
before touching the grid. On the tested PC (T1000, 4 GB) run scVI and scANVI one at a time; scGen
fits on the card once its decoder pass is chunked (see the notes), but fall back to
`CUDA_VISIBLE_DEVICES=""` if other GPU processes are resident. The smoke outputs in
`smoke_out/` are disposable and can be deleted once every method has come back green.

## Smoke test of the metrics step

The counterpart for the *scoring* half, and itself two tests - `smoke_test_metrics.py`
validates the **stack**, `smoke_test_metrics_pipeline.sh` the **pipeline**.

**The stack** (`smoke_test_metrics.py`, step 0). Neither test is optional bookkeeping:
`scib 1.1.7` predates the numpy/pandas versions in the current environment and fails *halfway
through* the metric computation, not at the start. This one checks that the environment can
compute every one of the 13 metrics on an identity integration, on a 5k fixture, before any real
job is launched. It has a SLURM wrapper (`submit_smoke_test.slurm`) so the cluster environment
can be validated the same way.

**The pipeline** (`smoke_test_metrics_pipeline.sh`). A **local** runner in `utils/`, kept
separate from the real cluster metrics step so no smoke logic ever leaks into the real run. It
walks the same grid but scores the tiny `smoke_out/` objects against `smoke_hvg.h5ad` instead of
the real integrations, running the whole metrics code path - the pre-flight check, both metric
jobs, the CSV tree, the merge - in a few minutes on 5,252 cells. What it proves, on top of the
stack test: that `metrics.py` / `metrics_kbet.py` read the real integration outputs (embed /
full / knn, from five different writers) correctly, that the per-type metric selection is right,
that kBET patches cleanly into the CSV `metrics.py` wrote, and that the tree merges into the
shape the plotting expects. Code path only - every expected metric produces a finite value -
not biological quality, which is meaningless at this size.

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

## Trained models

Three methods fit a model: **scVI**, **scANVI** and **scGen**. For those, training is the long
half of the run and what follows it is the fragile half - `batch_removal` for scGen, the latent
pass for the other two - so the fitted model is checkpointed in between and reloaded instead of
retrained on a re-run. A crash after training no longer costs the training.

The checkpoints do **not** go next to the integrated objects: `02_integration/` holds exactly
what the metrics read, and a model file among them would be neither an input nor an output of the
benchmark. Each method gets a flat directory of its own, matching the `02_drvi/` the DRVI notebook
already writes:

```
$DATA_DIR/02_scvi/scvi_unscaled_model.pt
$DATA_DIR/02_scanvi/scanvi_unscaled_scvi_model.pt      # stage 1
$DATA_DIR/02_scanvi/scanvi_unscaled_scanvi_model.pt    # stage 2
$DATA_DIR/02_scgen/scgen_unscaled_model.pt
$DATA_DIR/02_scgen/scgen_scaled_model.pt
```

The run id is in the file name, not in a subdirectory, so the two scGen runs share `02_scgen/`
without overwriting each other. scANVI has two files because it is two trainings: it initialises
from a full scVI fit, so a crash in the short second stage would otherwise redo the long first
one; a resume reloads whichever stages are already on disk.

`run_all.sh` and `submit_integration.slurm` pass the directory with `--model-dir`, since they own
the layout (as they do for `02_embeddings/`); `model_paths.py` holds the convention for running a
script by hand. **To force a retraining, delete the `.pt`** - `--force` alone re-runs the
integration but still reloads the saved model. The per-script `--retrain` flag does the same for a
single manual invocation.

## The run grid (`benchmark_grid.tsv`)

Every run is one row (17 rows). The file is the single source of truth: the local drivers and the
SLURM arrays both read it, so the benchmark matrix stays a readable piece of data instead of
logic scattered across scripts. A `types` value with two entries expands into two metrics jobs
downstream, which is how the 17 runs become 21 metrics jobs.

| Column | Meaning |
|---|---|
| `run_id` | `<method>_<scaling>`, e.g. `scanorama_scaled` - stable key for the output file, the figure folder and the job name. DRVI appends the latent size (`drvi_unscaled_128`), so runs at different `n_latent` do not overwrite each other; only the size kept for the benchmark is a row here (see "Choosing the latent size") |
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
| DRVI | Python (notebook) | embed | no (`n_latent=128`) |

**17 integration runs, 21 metrics jobs.** Scanorama and fastMNN write both a corrected matrix
and an embedding in a *single* execution, so they are scored twice from the same file by
changing `--type` only. scVI, scANVI and DRVI read `layers['counts']` and ignore `.X`
entirely: a scaled variant would be a duplicate run, so it is not in the grid.

Harmony runs through `scib.integration.harmony` (Python), not through R - the R `harmony`
package is not required anywhere in this phase.

### Choosing the latent size (DRVI)

DRVI is the only method with a free hyperparameter that changes what the benchmark scores, and
the criterion is visual: a *vanished* latent dimension carries no signal, so the number of
dimensions the model actually keeps says whether the latent space was large enough. That is why
DRVI is a notebook (`shiao_drvi_128.ipynb`), not a grid row driven by `run_all.sh`, and why it
is never in the SLURM array.

Three sizes were run (`N_LATENT` = 32, 64, 128), each with its own run id
(`drvi_unscaled_<N>`) so nothing overwrites anything. **128 is the run in the grid.** At 64 only
16 dimensions vanish (48 used); at 128, 63 vanish and 65 are used. The count of *used*
dimensions still grows from 48 to 65, so 64 was truncating the representation, while at 128 half
the space is left unused - the size is no longer the binding constraint, which is the signal to
stop there.

**The notebook.** Same cells and same 2,000 HVGs as every other method. `N_LATENT` at the top is
the only cell to edit: it drives the run id, every output and the figure folder. The model and
the latent space live in `$DATA_DIR/02_drvi/`; the benchmark cell writes the scored output
(`obsm['X_emb']` + `.npy`) exactly as `run_integration.py` does, and the notebook also draws the
DRVI-specific figures into `figures/<run_id>/`.

**The same run, headless (`run_drvi.py`).** Choosing the size needs a pair of eyes; training at a
chosen size does not, and it is the longest computation of the phase. `run_drvi.py` is the
notebook without the kernel - same input, same architecture, same seed, same early stopping, the
same four outputs and the same figures - so a size can be trained wherever it is convenient and
the notebook re-opened afterwards with `OVERWRITE = False`, which reads model and embedding from
disk instead of recomputing them. Resuming works the same way: an output already there is
reported `[have]` and reused, so a crash after training never costs the training.

```bash
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
python run_drvi.py                    # n_latent 128, the run in the grid
python run_drvi.py --n-latent 64      # another size, side by side with it
python run_drvi.py --overwrite        # retrain and rewrite everything

# on the cluster (CPU: no GPU there, hence partition `long`)
cd 02_integration_benchmark/02_2_integration && mkdir -p logs
sbatch --export=ALL,DATA_DIR=$DATA_DIR submit_drvi.slurm --n-latent 128
```

`submit_drvi.slurm` passes everything after the script path to `run_drvi.py` unchanged, activates
`$DRVI_ENV` (default `catalano_env`, which pins the same `drvi-py==0.2.7` as the local
environment) and writes its logs to `02_2_integration/logs/`. The input it needs is
`$DATA_DIR/shiao_hvg_2k.h5ad`: `utils/sync_to_cluster.sh --what inputs --method drvi` puts it
there. Since the outputs land in the layout the grid expects, 02_3 and 02_4 can then run on the
cluster without anything coming back first.

Nothing of the exploratory runs enters the benchmark: their integrated `.h5ad` and `.npy` are
deleted, so the grid cannot pick them up, and what is kept is `figures/drvi_unscaled_64/` alone,
as the evidence behind the choice. Their model and latent space under `$DATA_DIR/02_drvi/` are
disposable too - deleting them only means a retraining if that size is ever revisited.

DRVI on the **non-immune compartment** is a different question in a different phase, and the size
is re-chosen there: 176k cells and 18 labels are not this dataset, so `n_latent = 64` is the run of
`03_drvi_interpretation` (run ids `drvi_nonimm_<N>`, nothing to do with this grid row).

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
`knn` output (BBKNN) supports only the graph-based metrics - silhouette, PCR, cell cycle and
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
│                                          #       carries the 01_6 unintegrated UMAP
├── shiao_hvg_2k_scaled.h5ad               # 02_1  (1.1 GB, dense by construction)
│                                          #       fresh PCA + its own unintegrated UMAP,
│                                          #       both computed on the scaled matrix
├── shiao_hvg_2k.rds                       # 02_1  (2.9 GB, Seurat v3, data slot sparse)
├── shiao_hvg_2k_scaled.rds                # 02_1  (11 GB, Seurat v3, data slot dense)
├── shiao_hvg_2k_unintegrated_list.rds                  # 02_1  (the HVG symbols, for Seurat anchors)
├── 02_integration/<run_id>.h5ad           # 02_2  integrated objects
│                                          #       the R methods' <run_id>.rds intermediate is
│                                          #       deleted after conversion (--keep-rds to keep it)
│                                          # 02_3  appends obsm['X_umap'] in place, so a redraw
│                                          #       does not re-lay-out 620k cells
├── 02_embeddings/<run_id>.npy             # 02_2  latent spaces, kept for the figures
├── 02_scvi/<run_id>_model.pt              # 02_2  trained models, one flat file per run, kept
├── 02_scanvi/<run_id>_{scvi,scanvi}_model.pt  #      out of 02_integration/ (see "Trained models";
├── 02_scgen/<run_id>_model.pt             #      scANVI has one file per training stage)
├── 02_drvi/model_drvi_<N>.pt              # 02_2  DRVI only: the trained model, one flat file
├── 02_drvi/embed_drvi_unscaled_<N>.h5ad   # 02_2  DRVI only: latent space as .X, one var per
│                                          #       dimension with its stats and gene scores
│                                          #       (~730 MB at n_latent=128; read only by the
│                                          #       notebook, never by 02_3 / 02_4)
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

**Both integration and metrics can run locally or on SLURM** (see "The three drivers" for how).
The usual split is integration **locally** - the local environment is verified and has a GPU -
and metrics **on the cluster**, with integrated objects transferred as they are produced by
`utils/sync_to_cluster.sh`. The `--slurm` path exists for running integrations on the cluster
too: all methods land on partition `normal` (CPU), so scVI/scANVI/scGen are slower there than
locally, and scGen needs its own `$SCGEN_ENV` (an equivalent of the local `scgen-py`). DRVI is the
one method outside the array: the latent size is chosen in the notebook, and the run at that size
goes to the cluster - if it goes there at all - through `submit_drvi.slurm`.

Embeddings are also exported on their own (~120 MB each - 619,693 cells x 50 dims, float32 -
against several GB for a full object),
which keeps the transfer cheap for the `embed` methods and gives the figures step something
small and durable to read after the integrated objects are cleaned up.

### Running the metrics on the cluster

Three things have to be there: the **code**, the **objects** and the **environment**.

The code is a clone - `.gitignore` excludes `datasets/`, so it carries only scripts and figures.
The whole `02_integration_benchmark/` tree must stay intact: the wrappers resolve
`benchmark_grid.tsv` and `utils/` relative to their own location, and `metrics.py` /
`metrics_kbet.py` import `utils/metrics_shared.py` (which must import `utils/scib_compat.py`
before `scib`). `make_summary_table.R` additionally reads `utils/plotSingleTaskRNA.R`,
`knit_table.R` and `img/`.

```bash
# on the cluster
git clone git@github.com:Catt-Git/QCB-Master-Thesis.git ~/Tesi/QCB-Master-Thesis
```

The objects go to the cluster's `DATA_DIR` in the **same layout as the local one**, which is what
`utils/sync_to_cluster.sh` is for: it walks `benchmark_grid.tsv` and uploads the files it names
to the same relative path under the remote `DATA_DIR`, so nothing is typed by hand and adding a
grid row is enough to include it. `--what metrics` (default) sends the integrated objects plus
the references they are scored against - what 02_4 needs, so not the `.rds`, not `shiao.h5ad`,
not the checkpoints; `--what inputs` sends the 02_1 prepared inputs + HVG lists (what 02_2 needs
if the integrations run on the cluster); `--what all` both. Same filters, `--dry-run` and
resume-by-default as the other wrappers - a file already remote at the same byte size is
`[have]`, one still missing locally is `[skip]`, `--force` re-sends. Transfer is `rsync` (which
resumes a half-sent file in place), `--scp` for clusters without it; one SSH connection is opened
and reused, so authentication happens once.

```bash
# locally
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
export CLUSTER_HOST=<user>@<hostname>          # CLUSTER_DATA_DIR overrides the remote path
utils/sync_to_cluster.sh --scaling unscaled --dry-run
utils/sync_to_cluster.sh --scaling unscaled
```

Then submit. The environment is `catalano_env` (built from `environments/benchmark-hpc.yml`,
`module load Miniconda3/202411`), overridable with `$METRICS_ENV`:

```bash
# on the cluster
export DATA_DIR=/users/genomics/albertoc/Tesi/hopes_and_dreams/datasets
cd <clone>/02_integration_benchmark/02_4_metrics
mkdir -p logs
./run_all_metrics.sh --slurm --scaling unscaled --dry-run
./run_all_metrics.sh --slurm --scaling unscaled
```

> **`logs/` must exist before `sbatch`.** `#SBATCH --output=logs/...` is resolved against the
> directory `sbatch` was called from, and SLURM opens that file *before* the job body runs - the
> `mkdir -p logs` inside the two `.slurm` scripts is too late to help. Submit from
> `02_4_metrics/`, with the folder already created, or the array fails immediately.

## Figures

One folder per run, `figures/<run_id>/`: the five comparison panels every method gets from
`02_3`, plus - for DRVI only - the latent-space figures its notebook draws in the same folder.

### Per-method UMAP panels (`02_3_plot_method_umap`)

`plot_methods_umaps.py` writes `figures/<run_id>/` with five panels per run, from the integrated
object and the unintegrated reference:

- `cohort`, integrated (single panel)
- `cohort`, integrated vs unintegrated (two side by side)
- `cell_type`, integrated (single panel)
- `cell_type`, integrated vs unintegrated (two side by side)
- `cohort` + `cell_type`, both integrated (two side by side, no unintegrated)

These are a visual QC run **after** `run_all.sh` and **before** the metrics: a broken
integration (cohorts unmixed, or cell types destroyed) is obvious here in minutes, well before
the multi-hour metric jobs. The integrated embedding is derived from the output type
(`embed` → UMAP on `X_emb`; `full` → PCA on the corrected `.X` then UMAP; `knn` → UMAP on the
corrected graph); the unintegrated UMAP is the reference's `obsm['X_umap']`, not recomputed.
`cell_type` colours are inherited from `uns['cell_type_colors']` and `cohort` from a 
fixed palette, so a category keeps its colour across every panel and method.

**The integrated layout is cached back into the integrated object.** Laying out 619,693 cells
takes ~25 minutes, and recomputing it every time a colour or a panel changed was by far the most
expensive thing this step did. The first run writes `obsm['X_umap']` plus a
`uns['integrated_umap']` record of the `--type` and `--seed` it came from; later runs reuse it
when that record matches and recompute otherwise. `--recompute-umap` forces a fresh layout,
`--no-cache-umap` suppresses the write.

Three details make this safe rather than merely convenient:

- **Only `obsm['X_umap']` and the record are written** - never a PCA, never a neighbour graph.
  `metrics_shared.reduce_integrated` rebuilds both for itself, and persisting a graph would add
  hundreds of MB to objects that already reach 9.9 GB.
- **02_4 is unaffected**: every metrics path runs `reduce_data(umap=False)`, so no scib metric
  reads `obsm['X_umap']`. Adding it is inert.
- **The write is in place**, via h5py, not a read-modify-write through AnnData: rewriting a 9.9 GB
  object to add 5 MB would be minutes of I/O and a second copy on disk. `.X`, `.layers`, `.obs`,
  `.var` and `.obsp` are never touched.

The reuse test keys on the provenance record, not on the presence of `obsm['X_umap']`, because
several integrated objects already carry a layout of unknown origin from earlier exploratory work
(all six are unscaled Python methods). Those are recomputed and replaced rather than drawn.

**The `-u` reference must match the run's scaling variant here too**, exactly as in the metrics: a
scaled run is drawn against the scaled reference, so its "before" panel is the input the method
actually received, z-scoring included. This matters because per-batch scaling is itself a weak
batch correction. Plotting a scaled run against the *unscaled* baseline would show, as if the
method had produced it, mixing that the z-scoring had already done - and the two are not a small
difference: the scaled and unscaled `X_pca` are nearly unrelated spaces (PC1 of one against PC1
of the other correlates at |r| = 0.10, PC2 at 0.04), so the layouts they generate are genuinely
different pictures, not two renderings of one.

The scaled object holds that UMAP because `scale_batch.py` builds it - see below. A clean run of
the pipeline needs no extra step.

### The scaled object's unintegrated UMAP

`scale_batch.py` drops the `obsm['X_umap']`, `uns['neighbors']` and `.obsp` it inherits from
01_5/01_6, because all three describe the *unscaled* matrix and keeping them would attach a stale
layout to an object whose `.X` no longer matches. It then rebuilds the PCA **and the UMAP** on the
scaled matrix. The graph behind the layout is computed and discarded: nothing downstream reads a
neighbour graph off the reference, and on 619,693 cells it would add hundreds of MB to a file that
is already 1.1 GB.

The UMAP is there for 02_3 alone - `metrics.py` runs `reduce_data(umap=False)`, so no scib metric
reads it and 02_4 is indifferent. Its parameters are the ones 01_6 used (`n_neighbors=15`,
`metric='euclidean'`, `random_state=0`, on `X_pca`), taken from `shiao_hvg_2k.h5ad`'s
`uns['neighbors']['params']` rather than assumed, so the two unintegrated baselines differ by the
scaling alone and not by how they were laid out. Provenance is recorded in
`uns['unintegrated_umap']`.

Budget roughly an hour on 620k cells for the neighbours + layout, on top of the scaling itself.

> **If you already have a `shiao_hvg_2k_scaled.h5ad` from before this step existed**, do not
> re-run `scale_batch.py` to acquire the UMAP. Re-running rewrites `.X`, and that file is the
> exact input the scaled integrations consumed - deterministic or not, it must not be re-derived
> under finished runs. The layout was appended in place instead, reading only the stored
> `obsm['X_pca']` and leaving `.X`, `.layers`, `.obs` and `.var` byte-identical, so no integration
> and no metric had to be recomputed. Rebuilding from scratch is for a clean run, where nothing
> depends on the old file.

`plot_all.sh` drives the whole grid (see "The three drivers"):

```bash
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
cd 02_integration_benchmark/02_3_plot_method_umap
./plot_all.sh --dry-run          # preview what would be plotted
./plot_all.sh                    # every integrated run
./plot_all.sh --scaling unscaled # only the unscaled half
./plot_all.sh --force            # redraw, reusing the cached layouts (fast)
./plot_all.sh --force --recompute-umap   # redraw and lay the UMAPs out again (slow)
```

### DRVI latent-space panels (`02_2_integration`, notebook or `run_drvi.py`)

The DRVI run writes its own figures into the same `figures/<run_id>/` folder, each name
suffixed with the run id so the 64 and 128 panels stay distinguishable once pulled out of their
folder. They are *not* the method comparison - that is `02_3` above, identical for every method
- but the reading of the latent space itself:

- `umap_<key>_<run_id>.png` - one UMAP of the DRVI space per metadata/QC key, plus
  `umap_combined_<run_id>.png` (the 6-panel grid, same layout as the 01_6 unintegrated one) and
  `umap_per_cell_type_<run_id>.png` (one panel per CellTypist label, since ~50 labels in a
  single panel are unreadable).
- `latent_dimension_stats[_rmVanished]_<run_id>.png` - per-dimension reconstruction effect, max,
  mean, std, with and without the vanished dimensions: the plot behind the latent-size choice.
- `latent_dims_in_umap_<run_id>.png` and `latent_dims_in_heatmap_<key>_<run_id>.png` - each
  non-vanished dimension on the UMAP, and how the dimensions respond to `cell_type` (also sorted
  by label), `cohort`, `treatment`, `phase`.
- `ood_*_<run_id>.png` / `ind_linear_weighted_mean_<run_id>.png` - interpretability scores.
  OOD comes from the decoder reconstructions (fast, favours the genes *specific* to a dimension,
  and `OOD_min/max_possible` are its two halves); IND averages the effect of each factor over
  all cells (broader, a gene shared by several dimensions keeps a high score in all of them).

`figures/drvi_unscaled_64/` holds the same panels minus the interpretability scores and the five
`02_3` ones: it is the exploratory run kept only as evidence for the latent size.

### Final summary table (`02_4_metrics`)

The final summary table is produced by the official scIB plotting code
(`utils/plotSingleTaskRNA.R` + `utils/knit_table.R` + `utils/img/`, vendored from
theislab/scib-reproducibility), driven by `make_summary_table.R -i <merged.csv> -o figures`.
If those R packages (`dynutils`, `Hmisc`, `ggimage`) are missing it falls back to a built-in
scorer.

Overall score = 0.6 x bio conservation + 0.4 x batch correction, on min-max scaled metrics.

## Environments

| Where | Environment | Notes |
|---|---|---|
| local | `benchmark-py-r` | everything except scGen; must be **activated**, otherwise rpy2 cannot find R |
| local | `scgen-py` | scGen only; old pinned stack that conflicts with the main one. `run_scgen.py` calls `SCGEN` directly, since this env has no scib |
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
- **The `-u` reference must match the preprocessing variant** (scaled with scaled), in 02_3 as
  well as 02_4. Scoring - or plotting - a scaled output against an unscaled reference attributes
  part of the z-scoring effect to the integration method and makes the two columns
  non-comparable. Per-batch scaling is a weak batch correction in its own right, so the "before"
  has to be the input the method was actually handed.
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
  which conflicts with scanpy 1.12. `scale_batch.py` reimplements the same operation -
  `sc.pp.scale` per batch, scanpy defaults, no `max_value` clipping - so the result is
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
  pandas 3, which emits three encodings older readers reject - `nullable-string-array` (used
  for the *index*, i.e. barcodes and gene symbols), `nullable-boolean`, and `null` for `None`
  values inside `.uns`. `utils/h5ad_compat.py` downcasts them before writing and re-opens the
  file afterwards to verify, so no separate legacy export is needed and there is only ever one
  copy of each input.
- **scGen's decoder pass is chunked** (`run_scgen.py --decode-chunk`, default 16,384 cells).
  scGen decodes every cell in one forward pass inside `batch_removal`, which on 619k cells is a
  1.85 GB hidden activation plus a 4.7 GB output resident on the GPU at once - an immediate OOM
  on a 4 GB card, whatever the training batch size.
- **The R outputs are converted back with rpy2/anndata2ri** (`rds_to_h5ad.py`), **not** with
  `zellkonverter::writeH5AD`, which has no native R writer and would provision CPython through
  basilisk.
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
