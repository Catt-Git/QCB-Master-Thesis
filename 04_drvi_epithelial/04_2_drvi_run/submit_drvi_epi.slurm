#!/usr/bin/env bash
#SBATCH --job-name=drvi_epi
#SBATCH --output=logs/drvi_epi_%j.out
#SBATCH --error=logs/drvi_epi_%j.err
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#
# 04_2 DRVI run: the epithelial compartment on the cluster (CPU).
#
# A single job wrapping run_drvi_epi.py, which does exactly what drvi_epi.ipynb
# does with the same parameters. The notebook stays the place where the latent
# size is *judged*, from how many dimensions vanish; this is where a chosen size
# is trained when the local machine is not the right place for it.
#
# `long` rather than `normal`: the cluster has no GPU, so the training runs on
# CPU and is much slower than the local one. The lightest DRVI job of the thesis
# all the same - 74,441 cells against the 176,610 of 03_2 and the 619,693 of
# 02_2, and half the latent size - hence 48G. No `--time`, so the partition's own
# limit (30 days) applies.
#
# Environment on the cluster (override with the var below):
#   DRVI_ENV   conda env carrying scvi-tools + drvi-py   [default: catalano_env]
#
# catalano_env (benchmark-hpc.yml) pins drvi-py==0.2.7, the same version as the
# local benchmark-py-r, so the model this job writes loads back in the notebook.
#
# Everything after the script path reaches run_drvi_epi.py unchanged, so the
# latent size and the other switches are set there and not here:
#
#   cd 04_drvi_epithelial/04_2_drvi_run && mkdir -p logs
#   export DATA_DIR=/users/genomics/albertoc/Tesi/hopes_and_dreams/datasets
#   sbatch --export=ALL,DATA_DIR=$DATA_DIR submit_drvi_epi.slurm                # n_latent 64
#   sbatch --export=ALL,DATA_DIR=$DATA_DIR submit_drvi_epi.slurm --n-latent 32  # another size
#
# Two inputs, both written by 04_1 locally and both to be copied up before
# submitting (rsync/scp: the 02 helper utils/sync_to_cluster.sh only walks the
# phase-02 grid):
#   04_epi/shiao_epi_hvg_2k.h5ad   what DRVI trains on
#   04_epi/shiao_epi.h5ad          all genes, only read at the end to carry the
#                                  latent space into the downstream object
# The second is optional: without it the job says so at the start, trains, and
# writes everything except shiao_epi_<run_id>.h5ad.
#
# Outputs land where the notebook puts them (04_epi/model_<run_id>.pt,
# embed_<run_id>.h5ad, shiao_epi_<run_id>.h5ad) and the figures in the repo's
# figures/04_2_<run_id>/. Only the embedding has to come back for 04_3, which
# reads it and needs neither the model nor a GPU.

set -euo pipefail
mkdir -p logs

hostname; date
echo "Job name: ${SLURM_JOB_NAME:-}"; echo "Job ID: ${SLURM_JOB_ID:-}"
# Log the node, timestamp and job identifiers into the .out file for traceability.

: "${DATA_DIR:?set DATA_DIR (sbatch --export=ALL,DATA_DIR=... or export it first)}"
DRVI_ENV="${DRVI_ENV:-catalano_env}"

# SLURM copies this file to the node's spool directory before running it, so
# "$0"/BASH_SOURCE point there and NOT at the repo - they cannot locate
# run_drvi_epi.py. SLURM_SUBMIT_DIR (where sbatch was called from) covers the
# documented usage above; DRVI_RUN_DIR overrides it.
SCRIPT_DIR="${DRVI_RUN_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
RUNNER="$SCRIPT_DIR/run_drvi_epi.py"
[ -f "$RUNNER" ] || { echo "run_drvi_epi.py not found: $RUNNER (set DRVI_RUN_DIR, or sbatch from 04_2_drvi_run/)" >&2; exit 1; }

# torch reads this at import; keep it on the cores the scheduler granted rather
# than the node's full count.
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"

# `set -u` off across the whole conda block, the shell hook included: the hook
# auto-activates `base`, which sources this Miniconda's activate.d/openjdk_activate.sh,
# and that script expands ${JAVA_HOME} before defining it - fatal under -u, before
# `conda activate` is even reached. Same guard as in the phase-02 batch scripts.
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
