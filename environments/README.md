# environments

Conda/mamba environment files used across the thesis. Three environments, each with a
distinct role.

## The environments

| File | Env name | Where / when to use |
|---|---|---|
| `benchmark-py-r.yml` | `benchmark-py-r` | **Main local environment.** Python 3.12 + GPU (PyTorch, CUDA 12.4) + the R/Bioconductor stack (scran, batchelor, rpy2, anndata2ri) |
| `benchmark-hpc.yml` | `benchmark-hpc` | **Cluster environment.** Same scientific stack but CPU-only (no CUDA pin) and with **pinned pip versions** (`scib==1.1.7`, `drvi-py==0.2.7`) for reproducible SLURM jobs |
| `scgen-py.yml` | `scgen-py` | **Isolated legacy env for scGen only.** Old pinned deps (Python 3.9, `anndata=0.8`, `scanpy=1.9`, `scvi-tools 0.20.x`, `scgen==2.1.0`) that conflict with the main stack, so scGen gets its own env. |

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

To update an existing environment after editing a `.yml`:

```bash
conda env update -f environments/benchmark-py-r.yml --prune
```
