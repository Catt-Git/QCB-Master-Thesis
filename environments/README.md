# environments

Conda/mamba environment files used across the thesis. Three environments, each with a
distinct role.

## The environments

| File | Env name | Where / when to use |
|---|---|---|
| `benchmark-py-r.yml` | `benchmark-py-r` | **Main local environment.** Python 3.12 + GPU (PyTorch, CUDA 12.4) + the R/Bioconductor stack (scran, batchelor, rpy2, anndata2ri) |
| `benchmark-hpc.yml` | `benchmark-hpc` (on the HPC used for the thesis: `catalano_env`) | **Cluster environment.** Same scientific stack but CPU-only (no CUDA pin) and with **pinned pip versions** (`scib==1.1.7`, `drvi-py==0.2.7`) for reproducible SLURM jobs |
| `scgen-py.yml` | `scgen-py` | **Isolated legacy env for scGen only.** Old pinned deps (Python 3.9, `anndata 0.10`, `scanpy 1.9`, `torch 2.0`, `scvi-tools 0.19.0`, `scgen 2.1.0`) that conflict with the main stack, so scGen gets its own env. |

## Usage

Create an environment from its file (use `mamba` if available, it's faster):

```bash
conda env create -f environments/benchmark-py-r.yml
conda activate benchmark-py-r
```

Then run the code as described in each phase's README (`export DATA_DIR=...` first).

- Local machine: `benchmark-py-r`.
- Cluster / SLURM: `benchmark-hpc`.
- scGen integration method only: `scgen-py`.

> **Name of the cluster environment.** `benchmark-hpc.yml` declares `name: benchmark-hpc`, but on
> the HPC used for the thesis the environment was built as **`catalano_env`**, and that is the name
> hardcoded as default in every `.slurm` wrapper and in the SLURM-facing docs. Either create it
> under that name (`conda env create -f environments/benchmark-hpc.yml -n catalano_env`) or keep
> `benchmark-hpc` and override it per phase via the environment variables the wrappers read:
> `DOWNLOAD_ENV`, `PREPROC_ENV`, `INTEGRATION_ENV`, `METRICS_ENV`.

## Required post-install step: kBET

`kBET` is one of the 13 benchmark metrics, but it is only distributed on GitHub,
so no `.yml` can declare it. Run this once per environment that computes metrics
(`benchmark-py-r` locally, `benchmark-hpc` on the cluster):

```bash
Rscript -e "remotes::install_github('theislab/kBET')"
```

Skipping it does not raise an error - the metrics jobs simply come back without
kBET. `02_integration_benchmark/utils/smoke_test_metrics.py` reports whether the
package is present before computing anything, which is the cheap way to find out.
