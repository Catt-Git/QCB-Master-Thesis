#!/usr/bin/env bash
#
# 06_2 SoupX: the headless half of the phase, in one command.
#
# Three stages, two conda environments:
#
#   1. prepare_soupx_input.py   shiao.h5ad -> per-sample .mtx.gz + soup + clusters  [benchmark-py-r]
#   2. run_soupx.R (per sample) those inputs -> adjusted/<sample>/removed.mtx.gz    [soupx-r]
#   3. assemble_soupx.py        the removals -> shiao_soupx_all_cells.h5ad          [benchmark-py-r]
#
# Stage 2 is a loop: one process per channel, sequential. SoupX is fast (seconds to a
# minute per sample, against inferCNV's tens of minutes per cohort), so the 148 channels
# are a coffee break rather than an overnight run, and running them concurrently would only
# multiply the memory for no useful gain.
#
# Resuming is the default: a sample whose rho.csv exists is reported as [have] and skipped,
# so a run interrupted at channel 90 of 148 picks up at 90. --force re-runs everything.
#
# The ambient CALL is not here. Which cells have lost so much of their signal that they
# should not go into phase 05 is a threshold read off a distribution, so it lives in
# 06_3_ambient_qc/ambient_qc.ipynb, which reads what stage 3 writes. Same split as 05_1.
#
# Usage:
#   export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
#   ./soupx_all.sh                             # prepare + every channel + assemble, resuming
#   ./soupx_all.sh --force                     # re-run everything, overwriting
#   ./soupx_all.sh --dry-run                   # print what would run, do nothing
#   ./soupx_all.sh --samples P01_A_P P02_B_P   # only these channels
#   ./soupx_all.sh --plot                      # save autoEstCont's diagnostic per channel
#   ./soupx_all.sh --method soupOnly           # adjustCounts method [subtraction]
#   ./soupx_all.sh --clean                     # remove 06_amb/input/ once assembly succeeds
#   ./soupx_all.sh prepare                     # only stage 1 (also: soupx, assemble)
#
# Environments are overridable the way the wrappers of 01/02/05 do it:
#   PREP_ENV=benchmark-py-r  SOUPX_ENV=soupx-r
# Logs go to 06_2_soupx/logs/soupx_all_<timestamp>.log as well as to the terminal.
#
# Disk. The .mtx.gz handoff is ~4 GB for all 148 channels and the assembled object another
# ~4 GB; --clean drops the first once the second exists. See the phase README before
# starting if `df -h` is close.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PHASE_DIR="$(dirname "$SCRIPT_DIR")"
LOG_DIR="$SCRIPT_DIR/logs"

PREP_ENV="${PREP_ENV:-benchmark-py-r}"
SOUPX_ENV="${SOUPX_ENV:-soupx-r}"
METHOD="${METHOD:-subtraction}"

FORCE=0; DRY_RUN=0; PLOT=0; CLEAN=0
SAMPLES=()
STAGE_FILTER=()
while [ $# -gt 0 ]; do
  case "$1" in
    --force|-f) FORCE=1; shift ;;
    --dry-run|-n) DRY_RUN=1; shift ;;
    --plot) PLOT=1; shift ;;
    --clean) CLEAN=1; shift ;;
    --method) METHOD="$2"; shift 2 ;;
    --samples) shift; while [ $# -gt 0 ] && [[ "$1" != -* ]]; do SAMPLES+=("$1"); shift; done ;;
    -h|--help) sed -n '2,40p' "${BASH_SOURCE[0]}"; exit 0 ;;
    -*) echo "unknown option: $1" >&2; exit 1 ;;
    prepare|soupx|assemble) STAGE_FILTER+=("$1"); shift ;;
    *) echo "unknown stage '$1'; valid stages: prepare soupx assemble" >&2; exit 1 ;;
  esac
done

: "${DATA_DIR:?set DATA_DIR to the directory holding the datasets (outside the repo)}"
AMB_DIR="$DATA_DIR/06_amb"
export FIG_DIR="${FIG_DIR:-$PHASE_DIR/figures/06_2_soupx}"

run_stage() { [ ${#STAGE_FILTER[@]} -eq 0 ] || printf '%s\n' "${STAGE_FILTER[@]}" | grep -qx "$1"; }

# conda run rather than `conda activate`, so the two environments cannot leak into each
# other across the loop. Neither script uses rpy2, so full activation is not needed.
conda_run() { local env="$1"; shift; conda run --no-capture-output -n "$env" "$@"; }

mkdir -p "$AMB_DIR" "$LOG_DIR"
[ "$PLOT" -eq 1 ] && mkdir -p "$FIG_DIR"
LOG_FILE="$LOG_DIR/soupx_all_$(date +%Y%m%d_%H%M%S).log"
[ "$DRY_RUN" -eq 0 ] && exec > >(tee -a "$LOG_FILE") 2>&1

echo "DATA_DIR    : $DATA_DIR"
echo "output dir  : $AMB_DIR"
echo "FIG_DIR     : $FIG_DIR (plots: $PLOT)"
echo "environments: $PREP_ENV (python) | $SOUPX_ENV (SoupX)"
echo "method      : $METHOD"
[ "$DRY_RUN" -eq 0 ] && echo "log         : $LOG_FILE"
echo
df -h "$DATA_DIR" | tail -1
echo

overall_start=$SECONDS

# ---------------------------------------------------------------------------- stage 1
if run_stage prepare; then
  echo "==================================================================="
  echo ">>> prepare  (prepare_soupx_input.py)  $(date '+%Y-%m-%d %H:%M:%S')"
  echo "==================================================================="
  prep_args=()
  [ "$FORCE" -eq 1 ] && prep_args+=(--force)
  [ ${#SAMPLES[@]} -gt 0 ] && prep_args+=(--samples "${SAMPLES[@]}")
  if [ "$DRY_RUN" -eq 1 ]; then
    echo "[dry] conda run -n $PREP_ENV python3 prepare_soupx_input.py ${prep_args[*]:-}"
  else
    # ${a[@]+"${a[@]}"} and not "${a[@]:-}": the latter expands an EMPTY array to one
    # empty argument, which argparse then rejects as unrecognised.
    conda_run "$PREP_ENV" python3 "$SCRIPT_DIR/prepare_soupx_input.py" ${prep_args[@]+"${prep_args[@]}"}
  fi
  echo
fi

# ---------------------------------------------------------------------------- stage 2
if run_stage soupx; then
  CENSUS="$AMB_DIR/sample_census.csv"
  [ -f "$CENSUS" ] || { echo "missing $CENSUS; run the prepare stage first" >&2; exit 1; }

  if [ ${#SAMPLES[@]} -gt 0 ]; then
    TO_RUN=("${SAMPLES[@]}")
  else
    # Column 1 is `sample`, column 2 is `status`; only 'prepared' rows have inputs on disk.
    mapfile -t TO_RUN < <(awk -F, 'NR>1 && $2=="prepared" {print $1}' "$CENSUS")
  fi
  echo "channels to run: ${#TO_RUN[@]}"
  echo

  n_run=0; n_have=0; n_fail=0
  FAILED=()
  for sample in "${TO_RUN[@]}"; do
    rho_csv="$AMB_DIR/adjusted/$sample/rho.csv"
    if [ "$FORCE" -eq 0 ] && [ -f "$rho_csv" ]; then
      echo "[have] $sample: rho.csv already exists, skipping"
      n_have=$((n_have + 1)); continue
    fi
    r_args=(--sample "$sample" --method "$METHOD")
    [ "$FORCE" -eq 1 ] && r_args+=(--force)
    [ "$PLOT" -eq 1 ] && r_args+=(--plot)
    if [ "$DRY_RUN" -eq 1 ]; then
      echo "[dry] conda run -n $SOUPX_ENV Rscript run_soupx.R ${r_args[*]}"
      n_run=$((n_run + 1)); continue
    fi

    echo "--- $sample  $(date '+%H:%M:%S')"
    # A channel that fails must not take the other 147 down with it: autoEstCont can fail
    # on a degenerate sample, and that is a fact about the sample. Failures are counted and
    # listed at the end, and assemble refuses to run without them unless told otherwise.
    if conda_run "$SOUPX_ENV" Rscript "$SCRIPT_DIR/run_soupx.R" "${r_args[@]}"; then
      n_run=$((n_run + 1))
    else
      echo "[FAIL] $sample (exit $?), continuing" >&2
      n_fail=$((n_fail + 1))
      FAILED+=("$sample")
    fi
  done

  echo
  echo "==================================================================="
  echo "SoupX: $n_run run, $n_have already done, $n_fail failed"
  [ "${#FAILED[@]}" -gt 0 ] && echo "failed channels: ${FAILED[*]}"
  echo
fi

# ---------------------------------------------------------------------------- stage 3
if run_stage assemble; then
  echo "==================================================================="
  echo ">>> assemble  (assemble_soupx.py)  $(date '+%Y-%m-%d %H:%M:%S')"
  echo "==================================================================="
  if [ "$DRY_RUN" -eq 1 ]; then
    echo "[dry] conda run -n $PREP_ENV python3 assemble_soupx.py"
  else
    conda_run "$PREP_ENV" python3 "$SCRIPT_DIR/assemble_soupx.py"
    if [ "$CLEAN" -eq 1 ]; then
      # Only after the object exists: the input directory is reproducible from shiao.h5ad
      # in a few minutes, the assembled object is not.
      OUT="$AMB_DIR/shiao_soupx_all_cells.h5ad"
      if [ -s "$OUT" ]; then
        echo "--clean: removing $AMB_DIR/input"
        rm -rf "$AMB_DIR/input"
        df -h "$DATA_DIR" | tail -1
      else
        echo "--clean: $OUT is missing or empty, keeping the input directory" >&2
      fi
    fi
  fi
  echo
fi

echo "done in $(( (SECONDS - overall_start) / 60 ))m$(( (SECONDS - overall_start) % 60 ))s"
echo "next: 06_3_ambient_qc/ambient_qc.ipynb"
exit 0
