#!/usr/bin/env bash
#
# 05_2 subsetting: the headless half of the step, in one command.
#
# Runs the four scripts that turn the re-filtered raw subset written by subset_and_qc.ipynb
# into the definitive object, in order:
#
#   1. scran_norm_tum.py       <prefix>_raw.h5ad     -> <prefix>_norm.h5ad
#   2. cell_cycle_score_tum.py <prefix>_norm.h5ad    -> <prefix>_norm_cc.h5ad
#   3. reduce_data_tum.py      <prefix>_norm_cc.h5ad -> <prefix>_reduced.h5ad
#                                                    +  <prefix>_hvg_2k_list.csv
#                                                    +  <prefix>_hvg_2k.h5ad
#   4. clustering_tum.py       <prefix>_reduced.h5ad -> <prefix>.h5ad
#
# Everything runs locally in sequence: ~36k cells, half of 04 and a fifth of 03, so scran
# does not need the cluster. Each step is a separate process, which also means the R session
# rpy2 opens for scran is torn down before the reduction step starts.
#
# ## The cell set
#
# CELL_SET picks what the phase is about and, with it, the file prefix - so the two sets can
# be run one after the other in the same directory without ever overwriting each other:
#
#   CELL_SET=tum  (default)  the malignant cells        -> shiao_tum_*
#   CELL_SET=epi             all epithelium, post-CNV   -> shiao_epicnv_*
#
# `tum` is the primary line of the phase; `epi` is a control that answers "how much did the
# wrong labels cost phase 04?". See cell_set.py for the full argument. Note that resuming is
# per cell set, because the outputs are named differently: running `epi` after `tum` re-runs
# all four steps, as it should.
#
# Resuming is the default: a step whose output already exists is reported as [have] and
# skipped, so re-running after a crash picks up where it stopped. --force re-runs everything.
# The check is existence only, so delete an .h5ad truncated by a crash mid-write before
# resuming.
#
# Usage:
#   export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
#   ./subsetting_all.sh                 # every step, resuming, on the malignant subset
#   CELL_SET=epi ./subsetting_all.sh    # the same chain on all epithelium (the control)
#   ./subsetting_all.sh --force         # re-run everything, overwriting
#   ./subsetting_all.sh --dry-run       # print what would run, do nothing
#   ./subsetting_all.sh reduce cluster  # only the named step(s), in file order
#
# Step names: norm, cc, reduce, cluster.
# Logs go to 05_2_subsetting/logs/subsetting_all_<set>_<timestamp>.log and to the terminal.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="$SCRIPT_DIR/logs"
PYTHON="${PYTHON:-python3}"

export CELL_SET="${CELL_SET:-tum}"
case "$CELL_SET" in
  tum) PREFIX="shiao_tum" ;;
  epi) PREFIX="shiao_epicnv" ;;
  *) echo "CELL_SET must be 'tum' or 'epi', got '$CELL_SET'" >&2; exit 1 ;;
esac

STEP_NAMES=(norm cc reduce cluster)
STEP_SCRIPTS=(
  scran_norm_tum.py
  cell_cycle_score_tum.py
  reduce_data_tum.py
  clustering_tum.py
)
STEP_SUFFIXES=(
  _norm.h5ad
  _norm_cc.h5ad
  _reduced.h5ad
  .h5ad
)

FORCE=0; DRY_RUN=0
STEP_FILTER=()
while [ $# -gt 0 ]; do
  case "$1" in
    --force|-f) FORCE=1; shift ;;
    --dry-run|-n) DRY_RUN=1; shift ;;
    -h|--help) sed -n '2,45p' "${BASH_SOURCE[0]}"; exit 0 ;;
    -*) echo "unknown option: $1" >&2; exit 1 ;;
    *) STEP_FILTER+=("$1"); shift ;;
  esac
done

# Checked after the options, so --help works without an environment.
: "${DATA_DIR:?set DATA_DIR to the directory holding the datasets (outside the repo)}"
TUM_DIR="$DATA_DIR/05_tum"

# The Tirosh/Regev list needed by the cc step. Same default as 01_4 and 04_1.
export CC_GENES="${CC_GENES:-$DATA_DIR/regev_lab_cell_cycle_genes.txt}"

in_list() { local x="$1"; shift; local e; for e in "$@"; do [ "$e" = "$x" ] && return 0; done; return 1; }

for requested in "${STEP_FILTER[@]:-}"; do
  [ -n "$requested" ] || continue
  in_list "$requested" "${STEP_NAMES[@]}" || {
    echo "unknown step '$requested'; valid steps: ${STEP_NAMES[*]}" >&2; exit 1; }
done

# The input of the whole chain comes from the notebook, not from this script.
INPUT_H5AD="$TUM_DIR/${PREFIX}_raw.h5ad"
if [ ! -f "$INPUT_H5AD" ] && [ "$DRY_RUN" -eq 0 ]; then
  echo "missing $INPUT_H5AD" >&2
  echo "run subset_and_qc.ipynb first (with CELL_SET=$CELL_SET): it is the step that joins" >&2
  echo "the 05_1 call onto shiao.h5ad, takes the subset, restores the raw counts and" >&2
  echo "re-applies the gene/cell/cohort filters." >&2
  exit 1
fi

mkdir -p "$TUM_DIR" "$LOG_DIR"
LOG_FILE="$LOG_DIR/subsetting_all_${CELL_SET}_$(date +%Y%m%d_%H%M%S).log"
[ "$DRY_RUN" -eq 0 ] && exec > >(tee -a "$LOG_FILE") 2>&1

echo "DATA_DIR  : $DATA_DIR"
echo "CELL_SET  : $CELL_SET  (prefix $PREFIX)"
echo "output dir: $TUM_DIR"
echo "CC_GENES  : $CC_GENES"
echo "python    : $($PYTHON --version 2>&1)"
[ "$DRY_RUN" -eq 0 ] && echo "log       : $LOG_FILE"
echo

n_run=0; n_have=0; n_skip=0
overall_start=$SECONDS

for i in "${!STEP_NAMES[@]}"; do
  name="${STEP_NAMES[$i]}"
  script="$SCRIPT_DIR/${STEP_SCRIPTS[$i]}"
  output="$TUM_DIR/${PREFIX}${STEP_SUFFIXES[$i]}"

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
  ( cd "$SCRIPT_DIR" && "$PYTHON" "$script" )   # cd: the scripts import cell_set.py

  # A step that returns 0 without writing its output means the chain is broken somewhere it
  # did not raise; stop before the next step reads a stale file.
  [ -f "$output" ] || { echo "[error] $name finished but $output is missing" >&2; exit 1; }

  echo "[ok] $name in $(( (SECONDS - step_start) / 60 ))m$(( (SECONDS - step_start) % 60 ))s -> $output"
  n_run=$((n_run + 1))
done

echo "==================================================================="
echo "done in $(( (SECONDS - overall_start) / 60 ))m$(( (SECONDS - overall_start) % 60 ))s: \
$n_run step(s) run, $n_have already done, $n_skip filtered out"
[ "$n_have" -gt 0 ] && [ "$FORCE" -eq 0 ] && \
  echo "      (--force to re-run and overwrite the $n_have existing output(s))"
echo "next: visualization_tum.ipynb on $TUM_DIR/${PREFIX}.h5ad"
exit 0
