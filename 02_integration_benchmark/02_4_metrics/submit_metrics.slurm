#!/usr/bin/env bash
#SBATCH --job-name=scib_metrics
#SBATCH --output=logs/metrics_%A_%a.out
#SBATCH --error=logs/metrics_%A_%a.err
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --array=1-21%3
#
# 02_4 metrics: SLURM array, the 13 metrics minus kBET.
#
# One array task = one (run, output type). The task list is the grid expanded on
# its `types` column (17 runs -> 21 jobs), so --array=1-21. The `%3` throttle
# keeps at most three tasks running at once (at most three nodes in `normal`); the
# rest queue and start as those finish, i.e. effectively in sequence in batches of
# three. No `--time` is set: the partition's own default wall-clock limit applies.
# kBET is a separate array (submit_kbet.slurm) that runs afterwards and patches in.
#
#   sbatch submit_metrics.slurm
#
# Usually launched through 02_4_metrics/run_all_metrics.sh --slurm, which submits
# a filtered subset with the same %3 throttle and wires the kBET dependency.
# Requires DATA_DIR in the environment (or edit the default below) and the
# catalano_env conda environment on the node.

set -euo pipefail
mkdir -p logs

: "${DATA_DIR:?set DATA_DIR (sbatch --export=ALL,DATA_DIR=... or export it first)}"
ENV="${METRICS_ENV:-catalano_env}"
TASK="shiao"
BATCH_KEY="cohort"; LABEL_KEY="cell_type"; ORGANISM="human"

# SLURM copies this file to the node's spool directory before running it, so
# "$0"/BASH_SOURCE point there and NOT at the repo - they cannot locate
# benchmark_grid.tsv or the sibling scripts. run_all_metrics.sh exports
# METRICS_DIR; SLURM_SUBMIT_DIR (where sbatch was called from) covers a direct
# `sbatch` issued from 02_4_metrics/, as phase 00 does.
SCRIPT_DIR="${METRICS_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
GRID="$SCRIPT_DIR/../benchmark_grid.tsv"
[ -f "$GRID" ] || { echo "grid not found: $GRID (set METRICS_DIR, or sbatch from 02_4_metrics/)" >&2; exit 1; }
ROOT="$DATA_DIR/02_metrics"

# Expand the grid into one line per (run, type): the SLURM_ARRAY_TASK_ID-th line
# is this task. Columns: run_id method scaling output types reference hvgs.
JOB=$(awk -F'\t' 'NR>1 && $1!="" {
        n=split($8, ts, ",");
        for (i=1;i<=n;i++) print $1"\t"$2"\t"$5"\t"$7"\t"ts[i]"\t"$9"\t"$10
      }' "$GRID" | sed -n "${SLURM_ARRAY_TASK_ID}p")
[ -n "$JOB" ] || { echo "no job for array index $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
IFS=$'\t' read -r run_id method scaling output type reference hvgs <<< "$JOB"

int_path="$DATA_DIR/$output"
ref_path="$DATA_DIR/$reference"
out_dir="$ROOT/$TASK/metrics/$scaling/hvg"
csv="$out_dir/${method}_${type}.csv"
mkdir -p "$out_dir"

echo ">>> [$SLURM_ARRAY_TASK_ID] $run_id  ($method, $type)"

# `set -u` off across the whole conda block, the shell hook included: the hook
# auto-activates `base`, which sources this Miniconda's activate.d/openjdk_activate.sh,
# and that script expands ${JAVA_HOME} before defining it - fatal under -u, before
# `conda activate` is even reached. Same guard as conda_guarded() in run_all_metrics.sh.
set +u
module load Miniconda3/202411
eval "$(/soft/system/software/Miniconda3/202411/bin/conda shell.bash hook)"
conda activate "$ENV"
set -u

python "$SCRIPT_DIR/check_integrations.py" -i "$int_path" -u "$ref_path" \
  --types "$type" -b "$BATCH_KEY" -l "$LABEL_KEY"

python "$SCRIPT_DIR/metrics.py" -u "$ref_path" -i "$int_path" -o "$csv" \
  -m "$method" --type "$type" -b "$BATCH_KEY" -l "$LABEL_KEY" \
  --organism "$ORGANISM" --hvgs "$hvgs" -v

echo "[ok] $run_id $type -> $csv"
