#!/usr/bin/env bash
#SBATCH --job-name=soup_profile
#SBATCH --output=logs/soup_profile_%j.out
#SBATCH --error=logs/soup_profile_%j.err
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#
# 06_1: the empty-droplet profiles, on the cluster.
#
# The only step of phase 06 that has to run there, and for one reason: it reads
# `raw_feature_bc_matrix.h5`, the UNFILTERED Cell Ranger matrix. Those are millions of
# droplets per sample and tens of GB in total; they were never copied down to the laptop
# (`datasets/` holds objects built from the *filtered* matrices) and there is no reason to
# copy them now, because what this step produces is ~300 KB per sample.
#
# 16G is sized on the largest single raw matrix, which is read whole by `sc.read_10x_h5`
# and then summed over its droplet axis. That is far less than it sounds: a
# `raw_feature_bc_matrix.h5` of this dataset is ~6 MB on disk and a few hundred MB expanded,
# because the 737k-barcode whitelist is almost entirely empty. One sample at a time, so the
# peak is one matrix and not 148 of them; the loop is inside the Python script. `normal`
# rather than `long`: the whole batch is a matter of hours, not days. No `--time`, so the
# partition's own limit applies - the convention of every SLURM wrapper in this repo.
#
# CELLRANGER_DIR is worth setting explicitly. The default below follows what
# 00_5/run_cellranger_batch.sh documents ($DATA_DIR/cellranger_out), but on the cluster this
# thesis ran on the batch landed one level up, BESIDE datasets/:
#   /users/genomics/albertoc/Tesi/hopes_and_dreams/cellranger_out
#
# What to bring back afterwards (a few tens of MB for all 148 samples):
#     rsync -av <cluster>:$DATA_DIR/06_amb/soup/          $LOCAL_DATA_DIR/06_amb/soup/
#     rsync -av <cluster>:$DATA_DIR/06_amb/soup_census.csv $LOCAL_DATA_DIR/06_amb/
# Everything after this step is local.
#
# Environment (override with the vars below):
#   SOUP_ENV         conda env with scanpy    [default: catalano_env, as the other wrappers]
#   CELLRANGER_DIR   where <sample>/outs/ lives  [default: $DATA_DIR/cellranger_out]
#
# Usage:
#   cd 06_ambient_soupx/06_1_soup_profile && mkdir -p logs
#   export DATA_DIR=/users/genomics/albertoc/Tesi/hopes_and_dreams/datasets
#   sbatch --export=ALL,DATA_DIR=$DATA_DIR submit_soup_profile.slurm
#   sbatch --export=ALL,DATA_DIR=$DATA_DIR submit_soup_profile.slurm --force
#   sbatch --export=ALL,DATA_DIR=$DATA_DIR,CELLRANGER_DIR=$(dirname $DATA_DIR)/cellranger_out submit_soup_profile.slurm
# Everything after the script path reaches extract_soup_profile.py unchanged.

set -euo pipefail
mkdir -p logs

hostname; date
echo "Job name: ${SLURM_JOB_NAME:-}"; echo "Job ID: ${SLURM_JOB_ID:-}"

: "${DATA_DIR:?set DATA_DIR (sbatch --export=ALL,DATA_DIR=... or export it first)}"
SOUP_ENV="${SOUP_ENV:-catalano_env}"

# SLURM copies this file to the node's spool directory, so "$0" does not point at the repo.
# SLURM_SUBMIT_DIR is where sbatch was called from, which the usage above documents as
# 06_1_soup_profile/; SOUP_RUN_DIR overrides it.
SCRIPT_DIR="${SOUP_RUN_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
RUNNER="$SCRIPT_DIR/extract_soup_profile.py"
[ -f "$RUNNER" ] || { echo "extract_soup_profile.py not found: $RUNNER (set SOUP_RUN_DIR, or sbatch from 06_1_soup_profile/)" >&2; exit 1; }

CELLRANGER_DIR="${CELLRANGER_DIR:-$DATA_DIR/cellranger_out}"
[ -d "$CELLRANGER_DIR" ] || { echo "no Cell Ranger output at $CELLRANGER_DIR (set CELLRANGER_DIR)" >&2; exit 1; }
export CELLRANGER_DIR

# `set -u` off across the conda block, the shell hook included: the hook auto-activates
# `base`, which sources this Miniconda's activate.d/openjdk_activate.sh, and that script
# expands ${JAVA_HOME} before defining it - fatal under -u. Same guard as the 02/05 wrappers.
set +u
module load Miniconda3/202411
eval "$(/soft/system/software/Miniconda3/202411/bin/conda shell.bash hook)"
conda activate "$SOUP_ENV"
set -u

which python3
echo "env: $SOUP_ENV   CELLRANGER_DIR=$CELLRANGER_DIR   args: ${*:-<none>}"
echo

# Spelled out rather than a bare "$@": older bash treats an empty "$@" as unbound under
# `set -u`, and "${@:-}" would hand the script one empty argument argparse then rejects.
if [ $# -gt 0 ]; then python3 "$RUNNER" "$@"; else python3 "$RUNNER"; fi

date
