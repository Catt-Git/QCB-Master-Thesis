#!/usr/bin/env bash
#SBATCH --job-name=drvi
#SBATCH --output=logs/drvi_%j.out
#SBATCH --error=logs/drvi_%j.err
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#
# 02_2 integration: DRVI on the cluster (CPU).
#
# DRVI is not in submit_integration.slurm's array: the grid marks it `notebook`
# because the latent size is chosen by eye, and an array task cannot make that
# choice. Once the size IS chosen, the run underneath is an ordinary long
# computation, and this is it - a single job wrapping run_drvi.py, which does
# exactly what shiao_drvi_128.ipynb does with the same parameters.
#
# `long` rather than `normal`: the cluster has no GPU, so the 619,693-cell
# training runs on CPU and is much slower than the local one. No `--time`, so the
# partition's own limit (30 days) applies.
#
# Environment on the cluster (override with the var below):
#   DRVI_ENV   conda env carrying scvi-tools + drvi-py   [default: catalano_env]
#
# catalano_env (benchmark-hpc.yml) pins drvi-py==0.2.7, the same version as the
# local benchmark-py-r, so the model this job writes loads back in the notebook.
#
# Everything after the script path reaches run_drvi.py unchanged, so the latent
# size and the other switches are set there and not here:
#
#   cd 02_integration_benchmark/02_2_integration && mkdir -p logs
#   export DATA_DIR=/users/genomics/albertoc/Tesi/hopes_and_dreams/datasets
#   sbatch --export=ALL,DATA_DIR=$DATA_DIR submit_drvi.slurm                  # n_latent 128
#   sbatch --export=ALL,DATA_DIR=$DATA_DIR submit_drvi.slurm --n-latent 64    # another size
#
# The input is $DATA_DIR/shiao_hvg_2k.h5ad, so send it up first if the cluster
# does not have it yet:
#   utils/sync_to_cluster.sh --what inputs --method drvi
#
# Outputs land where the notebook puts them (02_drvi/, 02_integration/,
# 02_embeddings/) and the figures in the repo's figures/<run_id>/, so 02_3 and
# 02_4 can be run straight after, on the same cluster, without touching anything.

set -euo pipefail
mkdir -p logs

hostname; date
echo "Job name: ${SLURM_JOB_NAME:-}"; echo "Job ID: ${SLURM_JOB_ID:-}"
# Log the node, timestamp and job identifiers into the .out file for traceability.

: "${DATA_DIR:?set DATA_DIR (sbatch --export=ALL,DATA_DIR=... or export it first)}"
DRVI_ENV="${DRVI_ENV:-catalano_env}"

# SLURM copies this file to the node's spool directory before running it, so
# "$0"/BASH_SOURCE point there and NOT at the repo - they cannot locate
# run_drvi.py. SLURM_SUBMIT_DIR (where sbatch was called from) covers the
# documented usage above; INTEGRATION_DIR overrides it, as in
# submit_integration.slurm.
SCRIPT_DIR="${INTEGRATION_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
RUNNER="$SCRIPT_DIR/run_drvi.py"
[ -f "$RUNNER" ] || { echo "run_drvi.py not found: $RUNNER (set INTEGRATION_DIR, or sbatch from 02_2_integration/)" >&2; exit 1; }

# torch reads this at import; keep it on the cores the scheduler granted rather
# than the node's full count.
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"

# `set -u` off across the whole conda block, the shell hook included: the hook
# auto-activates `base`, which sources this Miniconda's activate.d/openjdk_activate.sh,
# and that script expands ${JAVA_HOME} before defining it - fatal under -u, before
# `conda activate` is even reached. Same guard as in submit_integration.slurm.
set +u
module load Miniconda3/202411
eval "$(/soft/system/software/Miniconda3/202411/bin/conda shell.bash hook)"
conda activate "$DRVI_ENV"
set -u

which python3
echo "env: $DRVI_ENV   OMP_NUM_THREADS=$OMP_NUM_THREADS   args: ${*:-<none>}"

# Spelled out rather than a bare "$@": older bash treats an empty "$@" as an
# unbound variable under `set -u`, and "${@:-}" would hand the script one empty
# argument, which argparse would reject.
if [ $# -gt 0 ]; then python3 "$RUNNER" "$@"; else python3 "$RUNNER"; fi

date
