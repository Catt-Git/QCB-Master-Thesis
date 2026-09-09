#!/usr/bin/env bash
#
# 05_1 sensitivity 2 of 2, the driver: prepare four reference configurations for two
# cohorts, run each, and write the comparison table. ~8 runs, ~60 min.
#
# Nothing under $DATA_DIR/05_tum/input or summary is read or written; everything lives in
# $DATA_DIR/05_tum/sensitivity/reference/.
#
#   export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
#   ./sensitivity_reference.sh

set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
: "${DATA_DIR:?set DATA_DIR}"
THREADS="${THREADS:-20}"
ROOT="$DATA_DIR/05_tum/sensitivity/reference"

conda run --no-capture-output -n "${PREP_ENV:-benchmark-py-r}" \
  python3 "$SCRIPT_DIR/sensitivity_reference.py" --mode prepare || exit 1

for cfg in "$ROOT"/*/*/; do
  [ -f "$cfg/annotations.tsv" ] || continue
  echo ">>> $(basename "$(dirname "$cfg")")/$(basename "$cfg")  $(date '+%H:%M:%S')"
  conda run --no-capture-output -n "${INFERCNV_ENV:-infercnv-r}" Rscript \
    "$SCRIPT_DIR/sensitivity_reference.R" --config-dir "${cfg%/}" --threads "$THREADS" \
    || echo "[FAIL] ${cfg%/}" >&2
done

conda run --no-capture-output -n "${PREP_ENV:-benchmark-py-r}" \
  python3 "$SCRIPT_DIR/sensitivity_reference.py" --mode analyze
