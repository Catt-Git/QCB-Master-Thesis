# environments

Conda/mamba environment files used across the thesis. Four environments, each with a
distinct role.

## The environments

| File | Env name | Where / when to use |
|---|---|---|
| `benchmark-py-r.yml` | `benchmark-py-r` | **Main local environment.** Python 3.12 + GPU (PyTorch, CUDA 12.4) + the R/Bioconductor stack (scran, batchelor, rpy2, anndata2ri) |
| `benchmark-hpc.yml` | `benchmark-hpc` (on the HPC used for the thesis: `catalano_env`) | **Cluster environment.** Same scientific stack but CPU-only (no CUDA pin) and with **pinned pip versions** (`scib==1.1.7`, `drvi-py==0.2.7`) for reproducible SLURM jobs |
| `scgen-py.yml` | `scgen-py` | **Isolated legacy env for scGen only.** Old pinned deps (Python 3.9, `anndata 0.10`, `scanpy 1.9`, `torch 2.0`, `scvi-tools 0.19.0`, `scgen 2.1.0`) that conflict with the main stack, so scGen gets its own env. |
| `cytotrace2-py.yml` | `cytotrace2-py` | **Isolated env for CytoTRACE2 only** (`04_drvi_epithelial/04_4_cytotrace2`). `cytotrace2-py` pins `numpy<2.0.0`, so it cannot live beside the main stack; CPU-only torch, Python 3.11. |

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
- CytoTRACE2 (04_4) only: `cytotrace2-py`.

> **Do not `pip install cytotrace2-py` into `benchmark-py-r`.** It declares `numpy<2.0.0` as a
> hard requirement, so pip does not refuse - it downgrades the main stack under you (numpy
> 2.4 → 1.26, pandas 3.0 → 2.3, scipy 1.18 → 1.17, zarr 3.2 → 3.1, anndata 0.13 → 0.12, scanpy
> 1.12 → 1.11) and leaves `fast-array-utils` and `tifffile` with requirements it cannot satisfy.
> Same class of conflict as scGen, same answer: its own environment. Should it happen anyway,
> the way back is
> `pip uninstall -y cytotrace2-py gdown` followed by
> `pip install numpy==2.4.6 scipy==1.18.0 pandas==3.0.3 zarr==3.2.1 anndata==0.13.2 scanpy==1.12.1`
> in one command, then `pip check`.

> **Name of the cluster environment.** `benchmark-hpc.yml` declares `name: benchmark-hpc`, but on
> the HPC used for the thesis the environment was built as **`catalano_env`**, and that is the name
> hardcoded as default in every `.slurm` wrapper and in the SLURM-facing docs. Either create it
> under that name (`conda env create -f environments/benchmark-hpc.yml -n catalano_env`) or keep
> `benchmark-hpc` and override it per phase via the environment variables the wrappers read:
> `DOWNLOAD_ENV`, `PREPROC_ENV`, `INTEGRATION_ENV`, `METRICS_ENV`.

## R inside the metrics environment

No metric of the 12 in `02_integration_benchmark` is computed through R, but the
metrics code imports `anndata2ri`, which starts an embedded R: the environment
still has to be *activated* (not just its `bin/python` invoked directly) so that
rpy2 can find R. `02_integration_benchmark/utils/smoke_test_metrics.py` prints
the R version it reaches, which is the cheap way to check this before launching
anything long.
