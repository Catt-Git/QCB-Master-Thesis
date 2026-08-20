#!/usr/bin/env bash
#SBATCH --job-name=ct2_epi
#SBATCH --output=logs/ct2_epi_%j.out
#SBATCH --error=logs/ct2_epi_%j.err
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#
# 04_4: CytoTRACE2 on the epithelial compartment, one patient at a time.
#
# This is the only part of 04_3 that gets a SLURM wrapper, and the reason it sits
# reason it is its own step rather than part of 04_5: its dependency is not in
# `benchmark-py-r`, its runtime is minutes-per-patient across 29 patients rather
# than seconds, and it is the one step here that is worth running somewhere else.
# Everything else in this step is a few minutes on the laptop.
#
# `long` rather than `normal`, as everywhere in this thesis, and no `--time`, so
# the partition's own limit (30 days) applies. 64G because the per-patient export
# densifies one cohort at a time (the largest is ~9k cells x 26,371 genes); the
# scoring itself is the lighter half.
#
# Environment on the cluster (override with the var below):
#   CT2_ENV    conda env carrying cytotrace2-py   [default: catalano_env]
#
# The package is NOT in benchmark-hpc.yml as of this writing. Install it once into
# whichever env this points at:
#   conda activate catalano_env && pip install cytotrace2-py
#
# `--cores` is worth setting: the package's own documentation recommends dropping
# to 1-2 cores on machines under 16 GB, and here the scheduler decides how many
# there are, so it is passed explicitly rather than left to the node's core count.
#
# Everything after the script path reaches cytotrace2_epi.py unchanged:
#
#   cd 04_drvi_epithelial/04_4_cytotrace2 && mkdir -p logs
#   export DATA_DIR=/users/genomics/albertoc/Tesi/hopes_and_dreams/datasets
#   sbatch --export=ALL,DATA_DIR=$DATA_DIR submit_cytotrace2_epi.slurm --cores 8
#   sbatch --export=ALL,DATA_DIR=$DATA_DIR submit_cytotrace2_epi.slurm --dry-run
#
# One input, to be copied up before submitting (rsync/scp: the 02 helper
# utils/sync_to_cluster.sh only walks the phase-02 grid):
#   04_epi/shiao_epi.h5ad     all genes; layers['counts'] is what is scored
#
# One output to bring back:
#   04_epi/cytotrace2_<run_id>.csv    per cell, what 04_5 reads
# The per-patient working matrices stay on the node and are deleted at the end
# unless --keep-work is passed.

set -euo pipefail
mkdir -p logs

hostname; date
echo "Job name: ${SLURM_JOB_NAME:-}"; echo "Job ID: ${SLURM_JOB_ID:-}"

: "${DATA_DIR:?set DATA_DIR (sbatch --export=ALL,DATA_DIR=... or export it first)}"
CT2_ENV="${CT2_ENV:-catalano_env}"

# SLURM copies this file to the node's spool directory, so "$0"/BASH_SOURCE point
# there and not at the repo. Same guard as submit_drvi_epi.slurm.
SCRIPT_DIR="${CT2_RUN_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
RUNNER="$SCRIPT_DIR/cytotrace2_epi.py"
[ -f "$RUNNER" ] || { echo "cytotrace2_epi.py not found: $RUNNER (set CT2_RUN_DIR, or sbatch from 04_4_cytotrace2/)" >&2; exit 1; }

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

# `set -u` off across the conda block: the shell hook auto-activates `base`, whose
# activate.d/openjdk_activate.sh expands ${JAVA_HOME} before defining it, which is
# fatal under -u. Same guard as in 04_2 and in the phase-02 batch scripts.
set +u
module load Miniconda3/202411
eval "$(/soft/system/software/Miniconda3/202411/bin/conda shell.bash hook)"
conda activate "$CT2_ENV"
set -u

which python3
echo "env: $CT2_ENV   OMP_NUM_THREADS=$OMP_NUM_THREADS   args: ${*:-<none>}"

if [ $# -gt 0 ]; then python3 "$RUNNER" "$@"; else python3 "$RUNNER"; fi

date
