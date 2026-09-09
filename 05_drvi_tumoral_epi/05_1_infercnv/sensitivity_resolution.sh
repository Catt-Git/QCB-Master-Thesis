#!/usr/bin/env bash
#
# 05_1 sensitivity 1 of 2: does the per-patient verdict depend on the subcluster granularity?
#
# The objection this answers is "the call depends on the resolution of your clustering".
# It cannot be answered by argument, only by sweeping the parameter and showing what the
# verdict does. The subclustering happens inside infercnv::run(), so every resolution costs
# a full run - hence four representative cohorts rather than all 33:
#
#   Patient16   8,028 epithelial cells, tumour-rich        (95% of epithelium CNV-positive)
#   Patient43   5,394                    intermediate      (~54%)
#   Patient64     502                    thin stromal null (~90%)
#   Patient52   9,243                    no detectable tumour, the cohort gate's test case
#
# Runs into a SEPARATE directory tree so the production summaries are never touched:
#   $DATA_DIR/05_tum/sensitivity/resolution/<cohort>_res<r>.csv
# analyse with sensitivity_resolution.py, which writes the tables the README quotes.
#
# Result as of the run in the README: over a 10x range of resolution the number of
# subclusters grows ~7x and the fraction of epithelium called CNV-positive moves by at most
# 0.035 - the verdict is identical for all four cohorts at all four resolutions. What does
# degrade is the aggregation itself: at 0.05 only 23-57% of epithelial cells sit in
# subclusters of >= 20 cells, so the median stops being an aggregate. That is why the
# operating value is the coarse end of the range, 0.005, and not simply "any of them".
#
# Usage:
#   export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
#   ./sensitivity_resolution.sh                       # ~75 min, 16 runs, resumes
#   ./sensitivity_resolution.sh --threads 20
#   python3 sensitivity_resolution.py                 # then the tables

set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
: "${DATA_DIR:?set DATA_DIR}"

THREADS="${THREADS:-20}"
COHORTS=(Patient64 Patient43 Patient16 Patient52)
RESOLUTIONS=(0.005 0.01 0.02 0.05)
[ "${1:-}" = "--threads" ] && THREADS="$2"

OUT="$DATA_DIR/05_tum/sensitivity/resolution"
SUMMARY="$DATA_DIR/05_tum/summary"
mkdir -p "$OUT"

for res in "${RESOLUTIONS[@]}"; do
  for coh in "${COHORTS[@]}"; do
    dest="$OUT/${coh}_res${res}.csv"
    [ -f "$dest" ] && { echo "[have] $coh res=$res"; continue; }
    # run_infercnv.R always writes to summary/<cohort>_cnv.csv; stash anything already
    # there, take the sweep's output out of the way, and put the original back.
    keep=""
    if [ -f "$SUMMARY/${coh}_cnv.csv" ]; then keep="$SUMMARY/${coh}_cnv.csv.sensitivity_bak"
      mv "$SUMMARY/${coh}_cnv.csv" "$keep"; fi
    echo ">>> $coh res=$res  $(date '+%H:%M:%S')"
    if conda run --no-capture-output -n "${INFERCNV_ENV:-infercnv-r}" Rscript \
         "$SCRIPT_DIR/run_infercnv.R" --cohort "$coh" --threads "$THREADS" --force \
         --leiden-resolution "$res" > "$OUT/${coh}_res${res}.log" 2>&1 \
       && [ -f "$SUMMARY/${coh}_cnv.csv" ]; then
      mv "$SUMMARY/${coh}_cnv.csv" "$dest"; echo "[ok] $coh res=$res"
    else
      echo "[FAIL] $coh res=$res - see $OUT/${coh}_res${res}.log" >&2
      tail -5 "$OUT/${coh}_res${res}.log" >&2
    fi
    [ -n "$keep" ] && mv "$keep" "$SUMMARY/${coh}_cnv.csv"
  done
done
echo "done - now: python3 sensitivity_resolution.py"
