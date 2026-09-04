#!/usr/bin/env bash
#
# 05_1 inferCNV: the headless half of the step, in one command.
#
# Two stages, two conda environments:
#
#   1. prepare_infercnv_input.py   shiao.h5ad -> per-cohort .mtx + annotations   [benchmark-py-r]
#   2. run_infercnv.R (per cohort) those inputs -> summary/<cohort>_cnv.csv      [infercnv-r]
#                                                + figures/05_1_infercnv/*.png
#
# Stage 2 is a loop: one process per cohort, sequential, because inferCNV holds the whole
# residual matrix of the cohort in memory and already parallelises internally over
# --threads. Running several cohorts at once buys little and multiplies both the memory and
# the working directory on disk.
#
# Resuming is the default: a cohort whose summary .csv exists is reported as [have] and
# skipped, so a run interrupted at cohort 20 of 33 picks up at 20. --force re-runs
# everything. Each cohort deletes its own working directory when it succeeds - see the disk
# note in the README; do not disable that with --keep-work for a full run.
#
# The malignant call is NOT here. It needs thresholds chosen by looking at the two-axis
# distribution, so it lives in call_malignant.ipynb, which reads what this script writes.
#
# Usage:
#   export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
#   ./infercnv_all.sh                          # prepare + every cohort, resuming
#   ./infercnv_all.sh --force                  # re-run everything, overwriting
#   ./infercnv_all.sh --dry-run                # print what would run, do nothing
#   ./infercnv_all.sh --threads 16             # threads for infercnv::run() [8]
#   ./infercnv_all.sh --cohorts Patient52 Patient16   # only these cohorts
#   ./infercnv_all.sh prepare                  # only stage 1
#   ./infercnv_all.sh infercnv                 # only stage 2
#
# Environments are overridable the way the SLURM wrappers of 01/02 do it:
#   PREP_ENV=benchmark-py-r  INFERCNV_ENV=infercnv-r
# Logs go to 05_1_infercnv/logs/infercnv_all_<timestamp>.log as well as to the terminal.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PHASE_DIR="$(dirname "$SCRIPT_DIR")"
LOG_DIR="$SCRIPT_DIR/logs"

PREP_ENV="${PREP_ENV:-benchmark-py-r}"
INFERCNV_ENV="${INFERCNV_ENV:-infercnv-r}"
THREADS="${THREADS:-8}"

FORCE=0; DRY_RUN=0; HMM=0
COHORTS=()
STAGE_FILTER=()
while [ $# -gt 0 ]; do
  case "$1" in
    --force|-f) FORCE=1; shift ;;
    --dry-run|-n) DRY_RUN=1; shift ;;
    --hmm) HMM=1; shift ;;
    --threads) THREADS="$2"; shift 2 ;;
    --cohorts) shift; while [ $# -gt 0 ] && [[ "$1" != -* ]]; do COHORTS+=("$1"); shift; done ;;
    -h|--help) sed -n '2,40p' "${BASH_SOURCE[0]}"; exit 0 ;;
    -*) echo "unknown option: $1" >&2; exit 1 ;;
    prepare|infercnv) STAGE_FILTER+=("$1"); shift ;;
    *) echo "unknown stage '$1'; valid stages: prepare infercnv" >&2; exit 1 ;;
  esac
done

: "${DATA_DIR:?set DATA_DIR to the directory holding the datasets (outside the repo)}"
CNV_DIR="$DATA_DIR/05_tum"
export FIG_DIR="${FIG_DIR:-$PHASE_DIR/figures/05_1_infercnv}"

run_stage() { [ ${#STAGE_FILTER[@]} -eq 0 ] || printf '%s\n' "${STAGE_FILTER[@]}" | grep -qx "$1"; }

# conda run rather than `conda activate`, so the two environments cannot leak into each
# other across the loop. Neither script uses rpy2, so full activation is not needed here
# (see the note in environments/README.md about why the metrics env is different).
conda_run() { local env="$1"; shift; conda run --no-capture-output -n "$env" "$@"; }

mkdir -p "$CNV_DIR" "$LOG_DIR" "$FIG_DIR"
LOG_FILE="$LOG_DIR/infercnv_all_$(date +%Y%m%d_%H%M%S).log"
[ "$DRY_RUN" -eq 0 ] && exec > >(tee -a "$LOG_FILE") 2>&1

echo "DATA_DIR    : $DATA_DIR"
echo "output dir  : $CNV_DIR"
echo "FIG_DIR     : $FIG_DIR"
echo "environments: $PREP_ENV (prepare) | $INFERCNV_ENV (inferCNV)"
echo "threads     : $THREADS | HMM: $HMM"
[ "$DRY_RUN" -eq 0 ] && echo "log         : $LOG_FILE"
echo
df -h "$DATA_DIR" | tail -1
echo

overall_start=$SECONDS

# ---------------------------------------------------------------------------- stage 1
if run_stage prepare; then
  echo "==================================================================="
  echo ">>> prepare  (prepare_infercnv_input.py)  $(date '+%Y-%m-%d %H:%M:%S')"
  echo "==================================================================="
  prep_args=()
  [ "$FORCE" -eq 1 ] && prep_args+=(--force)
  [ ${#COHORTS[@]} -gt 0 ] && prep_args+=(--cohorts "${COHORTS[@]}")
  if [ "$DRY_RUN" -eq 1 ]; then
    echo "[dry] conda run -n $PREP_ENV python3 prepare_infercnv_input.py ${prep_args[*]:-}"
  else
    # ${a[@]+"${a[@]}"} and not "${a[@]:-}": the latter expands an EMPTY array to one
    # empty argument, which argparse then rejects as an unrecognised argument.
    conda_run "$PREP_ENV" python3 "$SCRIPT_DIR/prepare_infercnv_input.py" ${prep_args[@]+"${prep_args[@]}"}
  fi
  echo
fi

# ---------------------------------------------------------------------------- stage 2
if run_stage infercnv; then
  CENSUS="$CNV_DIR/cohort_census.csv"
  [ -f "$CENSUS" ] || { echo "missing $CENSUS; run the prepare stage first" >&2; exit 1; }

  if [ ${#COHORTS[@]} -gt 0 ]; then
    TO_RUN=("${COHORTS[@]}")
  else
    # Column 1 is `cohort`, column 2 is `status`; only 'prepared' rows have inputs on disk.
    mapfile -t TO_RUN < <(awk -F, 'NR>1 && $2=="prepared" {print $1}' "$CENSUS")
  fi
  echo "cohorts to run: ${#TO_RUN[@]}"
  echo

  n_run=0; n_have=0; n_fail=0
  FAILED=()
  for cohort in "${TO_RUN[@]}"; do
    summary="$CNV_DIR/summary/${cohort}_cnv.csv"
    if [ "$FORCE" -eq 0 ] && [ -f "$summary" ]; then
      echo "[have] $cohort: $(basename "$summary") already exists, skipping"
      n_have=$((n_have + 1)); continue
    fi
    r_args=(--cohort "$cohort" --threads "$THREADS")
    [ "$FORCE" -eq 1 ] && r_args+=(--force)
    [ "$HMM" -eq 1 ] && r_args+=(--hmm)
    if [ "$DRY_RUN" -eq 1 ]; then
      echo "[dry] conda run -n $INFERCNV_ENV Rscript run_infercnv.R ${r_args[*]}"
      n_run=$((n_run + 1)); continue
    fi

    echo "==================================================================="
    echo ">>> $cohort  $(date '+%Y-%m-%d %H:%M:%S')"
    echo "==================================================================="
    step_start=$SECONDS
    # A cohort that fails must not take the other 32 down with it: inferCNV can fail on a
    # cohort whose reference is degenerate, and that is a fact about that cohort, not a
    # reason to stop. Failures are counted and listed at the end.
    if conda_run "$INFERCNV_ENV" Rscript "$SCRIPT_DIR/run_infercnv.R" "${r_args[@]}"; then
      echo "[ok] $cohort in $(( (SECONDS - step_start) / 60 ))m$(( (SECONDS - step_start) % 60 ))s"
      n_run=$((n_run + 1))
    else
      echo "[FAIL] $cohort (exit $?), continuing" >&2
      n_fail=$((n_fail + 1))
      FAILED+=("$cohort")
    fi
    # The working directory is removed by run_infercnv.R itself on success; on failure it
    # is left behind on purpose, so remove it here before the next cohort fills the disk.
    [ -d "$CNV_DIR/work/$cohort" ] && rm -rf "$CNV_DIR/work/$cohort"
    df -h "$DATA_DIR" | tail -1
  done

  echo "==================================================================="
  echo "inferCNV: $n_run run, $n_have already done, $n_fail failed"
  [ "${#FAILED[@]}" -gt 0 ] && echo "failed cohorts: ${FAILED[*]}"
fi

echo "done in $(( (SECONDS - overall_start) / 60 ))m$(( (SECONDS - overall_start) % 60 ))s"
echo "next: call_malignant.ipynb on $CNV_DIR/summary/"
exit 0
