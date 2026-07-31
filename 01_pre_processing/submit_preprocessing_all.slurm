#!/usr/bin/env bash
#SBATCH --job-name=preprocessing_all
#SBATCH --output=logs/preprocessing_all_%j.out
#SBATCH --error=logs/preprocessing_all_%j.err
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=470G
#
# 01 pre-processing: the whole 01_3 -> 01_5 chain as ONE job.
#
# The memory is sized on the heaviest step (scran on 619k cells, the same 470G 
# submit_scran_norm.slurm asked for.
#
# Nothing about the pipeline lives here: this only sets up the environment and calls
# preprocessing_all.sh, which owns the step list, the [have]/[run] resume logic and
# the pre-flight input checks. Resuming therefore works exactly as it does locally:
# resubmit after a failure and the finished steps are skipped.
#
# Usually launched through the runner, which fills in every export and redirects the
# logs to 01_pre_processing/logs/ with absolute paths:
#   ./preprocessing_all.sh --slurm
#   ./preprocessing_all.sh --slurm --force
#   ./preprocessing_all.sh --slurm reduce cluster
#
# Direct use, submitted from 01_pre_processing/ so the relative log paths resolve.
# Anything after the script path reaches the runner unchanged (--force, step names):
#   sbatch --export=ALL,DATA_DIR=/users/genomics/albertoc/Tesi/hopes_and_dreams/datasets \
#          submit_preprocessing_all.slurm
#   sbatch --export=ALL,DATA_DIR=... submit_preprocessing_all.slurm --force reduce cluster
#
# Variables read (all optional except DATA_DIR):
#   DATA_DIR     the .h5ad chain's directory                    [required]
#   PREPROC_DIR  the repo's 01_pre_processing/ (see below)       [default: SLURM_SUBMIT_DIR]
#   CC_GENES     Tirosh/Regev list                              [default: $DATA_DIR/regev_lab_cell_cycle_genes.txt]
#   PREPROC_ENV  conda env with the rpy2/scran/scib/celltypist stack       [default: catalano_env]
#   FIG_DIR      where scran_norm.py drops its two diagnostics             [default: the repo's
#                figures/01_3_normalization; --export=ALL forwards it if you export it first]

set -euo pipefail

hostname; date
echo "Job name: ${SLURM_JOB_NAME:-}"; echo "Job ID: ${SLURM_JOB_ID:-}"
# Log the node, timestamp and job identifiers into the .out file for traceability.

: "${DATA_DIR:?set DATA_DIR (sbatch --export=ALL,DATA_DIR=... or use ./preprocessing_all.sh --slurm)}"

# SLURM copies this file to the node's spool directory before running it, so
# "$0"/BASH_SOURCE point there and NOT at the repo. preprocessing_all.sh exports
# PREPROC_DIR; SLURM_SUBMIT_DIR (where sbatch was called from) covers a direct
# `sbatch` issued from 01_pre_processing/, as phase 00 does.
SCRIPT_DIR="${PREPROC_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
RUNNER="$SCRIPT_DIR/preprocessing_all.sh"
[ -x "$RUNNER" ] || { echo "runner not found or not executable: $RUNNER (set PREPROC_DIR, or sbatch from 01_pre_processing/)" >&2; exit 1; }

PREPROC_ENV="${PREPROC_ENV:-catalano_env}"
export CC_GENES="${CC_GENES:-$DATA_DIR/regev_lab_cell_cycle_genes.txt}"
# scran's R side and scib's numpy/BLAS calls both read this; keep it on the cores the
# scheduler actually granted rather than the node's full count.
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

# `set -u` off across the whole conda block, the shell hook included: the hook
# auto-activates `base`, which sources this Miniconda's activate.d/openjdk_activate.sh,
# and that script expands ${JAVA_HOME} before defining it - fatal under -u, before
# `conda activate` is even reached. Same guard as in 02's batch scripts.
set +u
module load Miniconda3/202411
eval "$(/soft/system/software/Miniconda3/202411/bin/conda shell.bash hook)"
conda activate "$PREPROC_ENV"
set -u

which python3
echo "env: $PREPROC_ENV   OMP_NUM_THREADS=$OMP_NUM_THREADS   runner args: ${*:-<none>}"

# Spelled out rather than a bare "$RUNNER" "$@": older bash treats an empty "$@" as an
# unbound variable under `set -u`, and the usual "${@:-}" workaround would hand the
# runner one empty argument, which it would read as a (never matching) step name.
if [ $# -gt 0 ]; then "$RUNNER" "$@"; else "$RUNNER"; fi

date
