#!/usr/bin/env bash
#SBATCH --job-name=scib_kbet
#SBATCH --output=logs/kbet_%A_%a.out
#SBATCH --error=logs/kbet_%A_%a.err
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --array=1-21%3
#
# 02_4 metrics: SLURM array, kBET alone.
#
# Same expansion as submit_metrics.slurm (17 runs -> 21 jobs), same `%3` throttle
# (at most three nodes in `normal` at once) and no `--time` (the partition default
# applies). Each task patches its kBET value into the CSV metrics.py already wrote,
# so run this AFTER the metrics array:
#
#   m=$(sbatch --parsable submit_metrics.slurm)
#   sbatch --dependency=afterok:$m submit_kbet.slurm

set -euo pipefail
mkdir -p logs

: "${DATA_DIR:?set DATA_DIR (sbatch --export=ALL,DATA_DIR=... or export it first)}"
ENV="${METRICS_ENV:-catalano_env}"
TASK="shiao"
BATCH_KEY="cohort"; LABEL_KEY="cell_type"

# SLURM copies this file to the node's spool directory before running it, so
# "$0"/BASH_SOURCE point there and NOT at the repo. See submit_metrics.slurm.
SCRIPT_DIR="${METRICS_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
GRID="$SCRIPT_DIR/../benchmark_grid.tsv"
[ -f "$GRID" ] || { echo "grid not found: $GRID (set METRICS_DIR, or sbatch from 02_4_metrics/)" >&2; exit 1; }
ROOT="$DATA_DIR/02_metrics"

JOB=$(awk -F'\t' 'NR>1 && $1!="" {
        n=split($8, ts, ",");
        for (i=1;i<=n;i++) print $1"\t"$2"\t"$5"\t"$7"\t"ts[i]"\t"$9
      }' "$GRID" | sed -n "${SLURM_ARRAY_TASK_ID}p")
[ -n "$JOB" ] || { echo "no job for array index $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
IFS=$'\t' read -r run_id method scaling output type reference <<< "$JOB"

int_path="$DATA_DIR/$output"
ref_path="$DATA_DIR/$reference"
csv="$ROOT/$TASK/metrics/$scaling/hvg/${method}_${type}.csv"

echo ">>> [$SLURM_ARRAY_TASK_ID] kBET $run_id ($method, $type)"

# `set -u` off across the whole conda block: see submit_metrics.slurm (openjdk_activate.sh).
set +u
module load Miniconda3/202411
eval "$(/soft/system/software/Miniconda3/202411/bin/conda shell.bash hook)"
conda activate "$ENV"
set -u

python "$SCRIPT_DIR/metrics_kbet.py" -u "$ref_path" -i "$int_path" -o "$csv" \
  -m "$method" --type "$type" -b "$BATCH_KEY" -l "$LABEL_KEY"

echo "[ok] kBET $run_id $type -> $csv"
