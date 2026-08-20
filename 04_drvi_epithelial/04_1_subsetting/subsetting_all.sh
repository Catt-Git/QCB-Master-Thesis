#!/usr/bin/env bash
#
# 04_1 subsetting: the headless half of the phase, in one command.
#
# Runs the four scripts that turn the re-filtered raw subset written by
# subset_and_qc.ipynb into the definitive epithelial object, in order:
#
#   1. scran_norm_epi.py       shiao_epi_raw.h5ad     -> shiao_epi_norm.h5ad
#   2. cell_cycle_score_epi.py shiao_epi_norm.h5ad    -> shiao_epi_norm_cc.h5ad
#   3. reduce_data_epi.py      shiao_epi_norm_cc.h5ad -> shiao_epi_reduced.h5ad
#                                                     +  shiao_epi_hvg_2k_list.csv
#                                                     +  shiao_epi_hvg_2k.h5ad
#   4. clustering_epi.py       shiao_epi_reduced.h5ad -> shiao_epi.h5ad
#
# Everything runs locally in sequence: ~74k cells, so unlike 01_3 scran does not
# need the cluster. Each step is a separate process, which also means the R
# session rpy2 opens for scran is torn down before the reduction step starts.
#
# Resuming is the default: a step whose output already exists is reported as
# [have] and skipped, so re-running after a crash picks up where it stopped.
# --force re-runs everything. The check is existence only, so delete an .h5ad
# truncated by a crash mid-write before resuming.
#
#
# Usage:
#   export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
#   ./subsetting_all.sh                 # every step, resuming
#   ./subsetting_all.sh --force         # re-run everything, overwriting
#   ./subsetting_all.sh --dry-run       # print what would run, do nothing
#   ./subsetting_all.sh reduce cluster  # only the named step(s), in file order
#
# Step names: norm, cc, reduce, cluster.
# Logs go to 04_1_subsetting/logs/subsetting_all_<timestamp>.log as well as to
# the terminal.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="$SCRIPT_DIR/logs"
PYTHON="${PYTHON:-python3}"

# One row per step: name, script, output that marks it done.
STEP_NAMES=(norm cc reduce cluster)
STEP_SCRIPTS=(
  scran_norm_epi.py
  cell_cycle_score_epi.py
  reduce_data_epi.py
  clustering_epi.py
)
STEP_OUTPUTS=(
  shiao_epi_norm.h5ad
  shiao_epi_norm_cc.h5ad
  shiao_epi_reduced.h5ad
  shiao_epi.h5ad
)

FORCE=0; DRY_RUN=0
STEP_FILTER=()
while [ $# -gt 0 ]; do
  case "$1" in
    --force|-f) FORCE=1; shift ;;
    --dry-run|-n) DRY_RUN=1; shift ;;
    -h|--help) sed -n '2,40p' "${BASH_SOURCE[0]}"; exit 0 ;;
    -*) echo "unknown option: $1" >&2; exit 1 ;;
    *) STEP_FILTER+=("$1"); shift ;;
  esac
done

# Checked after the options, so --help works without an environment.
: "${DATA_DIR:?set DATA_DIR to the directory holding the datasets (outside the repo)}"
EPI_DIR="$DATA_DIR/04_epi"

# The Tirosh/Regev list needed by the cc step. Same default as 01_4; override by
# exporting CC_GENES yourself.
export CC_GENES="${CC_GENES:-$DATA_DIR/regev_lab_cell_cycle_genes.txt}"

in_list() { local x="$1"; shift; local e; for e in "$@"; do [ "$e" = "$x" ] && return 0; done; return 1; }

# Validate the step names before doing any work, so a typo fails immediately
# rather than after the first hour of scran.
for requested in "${STEP_FILTER[@]:-}"; do
  [ -n "$requested" ] || continue
  in_list "$requested" "${STEP_NAMES[@]}" || {
    echo "unknown step '$requested'; valid steps: ${STEP_NAMES[*]}" >&2; exit 1; }
done

# The input of the whole chain comes from the notebook, not from this script.
INPUT_H5AD="$EPI_DIR/shiao_epi_raw.h5ad"
if [ ! -f "$INPUT_H5AD" ] && [ "$DRY_RUN" -eq 0 ]; then
  echo "missing $INPUT_H5AD" >&2
  echo "run subset_and_qc.ipynb first: it is the step that subsets the epithelial" >&2
  echo "cells, restores the raw counts and re-applies the gene/cell/cohort filters." >&2
  exit 1
fi

mkdir -p "$EPI_DIR" "$LOG_DIR"
LOG_FILE="$LOG_DIR/subsetting_all_$(date +%Y%m%d_%H%M%S).log"
[ "$DRY_RUN" -eq 0 ] && exec > >(tee -a "$LOG_FILE") 2>&1

echo "DATA_DIR  : $DATA_DIR"
echo "output dir: $EPI_DIR"
echo "CC_GENES  : $CC_GENES"
echo "python    : $($PYTHON --version 2>&1)"
[ "$DRY_RUN" -eq 0 ] && echo "log       : $LOG_FILE"
echo

n_run=0; n_have=0; n_skip=0
overall_start=$SECONDS

for i in "${!STEP_NAMES[@]}"; do
  name="${STEP_NAMES[$i]}"
  script="$SCRIPT_DIR/${STEP_SCRIPTS[$i]}"
  output="$EPI_DIR/${STEP_OUTPUTS[$i]}"

  if [ ${#STEP_FILTER[@]} -gt 0 ] && ! in_list "$name" "${STEP_FILTER[@]}"; then
    n_skip=$((n_skip + 1)); continue
  fi

  if [ "$FORCE" -eq 0 ] && [ -f "$output" ]; then
    echo "[have] $name: $output already exists, skipping"
    n_have=$((n_have + 1)); continue
  fi

  if [ "$DRY_RUN" -eq 1 ]; then
    echo "[dry] $name: $PYTHON $script -> $output"
    n_run=$((n_run + 1)); continue
  fi

  echo "==================================================================="
  echo ">>> $name  ($(basename "$script"))  $(date '+%Y-%m-%d %H:%M:%S')"
  echo "==================================================================="
  step_start=$SECONDS
  "$PYTHON" "$script"

  # A step that returns 0 without writing its output means the chain is broken
  # somewhere it did not raise; stop before the next step reads a stale file.
  [ -f "$output" ] || { echo "[error] $name finished but $output is missing" >&2; exit 1; }

  echo "[ok] $name in $(( (SECONDS - step_start) / 60 ))m$(( (SECONDS - step_start) % 60 ))s -> $output"
  n_run=$((n_run + 1))
done

echo "==================================================================="
echo "done in $(( (SECONDS - overall_start) / 60 ))m$(( (SECONDS - overall_start) % 60 ))s: \
$n_run step(s) run, $n_have already done, $n_skip filtered out"
[ "$n_have" -gt 0 ] && [ "$FORCE" -eq 0 ] && \
  echo "      (--force to re-run and overwrite the $n_have existing output(s))"
echo "next: visualization_epi.ipynb on $EPI_DIR/shiao_epi.h5ad"
exit 0
