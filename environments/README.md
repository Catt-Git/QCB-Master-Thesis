# environments

Conda/mamba environment files used across the thesis. Three environments, each with a
distinct role.

## The environments

| File | Env name | Where / when to use |
|---|---|---|
| `benchmark-py-r.yml` | `benchmark-py-r` | **Main local environment.** Python 3.12 + GPU (PyTorch, CUDA 12.4) + the R/Bioconductor stack (scran, batchelor, rpy2, anndata2ri) |
| `benchmark-hpc.yml` | `benchmark-hpc` | **Cluster environment.** Same scientific stack but CPU-only (no CUDA pin) and with **pinned pip versions** (`scib==1.1.7`, `drvi-py==0.2.7`) for reproducible SLURM jobs |
| `scgen-py.yml` | `scgen-py` | **Isolated legacy env for scGen only.** Old pinned deps (Python 3.9, `anndata 0.10`, `scanpy 1.9`, `torch 2.0`, `scvi-tools 0.20.3`, `scgen 2.1.0`) that conflict with the main stack, so scGen gets its own env. |

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

## Required post-install step: kBET

`kBET` is one of the 13 benchmark metrics, but it is only distributed on GitHub,
so no `.yml` can declare it. Run this once per environment that computes metrics
(`benchmark-py-r` locally, `benchmark-hpc` on the cluster):

```bash
Rscript -e "remotes::install_github('theislab/kBET')"
```

Skipping it does not raise an error — the metrics jobs simply come back without
kBET. `02_integration_benchmark/utils/smoke_test_metrics.py` reports whether the
package is present before computing anything, which is the cheap way to find out.

## Rebuilding `scgen-py`

This one environment must be **rebuilt from scratch**, never updated in place:

```bash
conda env remove -n scgen-py
conda env create -f environments/scgen-py.yml
conda activate scgen-py
python environments/check_scgen_env.py     # must exit 0
```

The reason is recorded at the top of `scgen-py.yml`. In short: conda and pip both
managing the same Python packages is what broke the previous build, so here conda
provides only the interpreter and pip resolves everything else in one pass.
`conda env update` would reintroduce exactly the mixed ownership that caused the
failure.

`check_scgen_env.py` is not a formality. A partially removed numpy still satisfies
`import numpy`, so the environment can look healthy and fail hours into a run; the
script checks the installation itself, the numpy/pandas ABI, that this env can read
`.h5ad` files written by the main one, and that the full scGen call path runs.
