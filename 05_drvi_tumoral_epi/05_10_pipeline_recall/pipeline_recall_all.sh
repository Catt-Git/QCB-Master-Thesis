#!/usr/bin/env bash
#
# 05_10: the standard pipeline against the factor pipeline, on all three collections.
#
#   1. pipeline_recall_tum.py --collection scie
#   2. pipeline_recall_tum.py --collection emt
#   3. pipeline_recall_tum.py --collection gavish
#
# One step, three collections, its own driver - as 05_9. The chain 05_4 - 05_8 is Route A and
# Route B on DRVI and cannot answer for a space with no decoder; this step does not ask it to,
# it runs the clustering pipeline end to end and compares the two OUTPUTS.
#
# WHAT IT READS AND NEVER WRITES: 05_4's .gmt, 05_6's per-cell scores, 05_3's embedding and
# 05_3b's, all read-only. It owns its own tables and figures plus one cache per arm in
# $DATA_DIR/05_tum/ (`pipeline_lists_<arm>_<run>.json` and the partitions beside it). The
# cache is what makes the three collections cheap: the Leiden partitions and their DE do not
# depend on the collection, so they are computed once and reused - and dropped whole if k,
# the list depth or either DE filter moves.
#
# Usage:
#   export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
#   ./pipeline_recall_all.sh                             # all three collections, resuming
#   ./pipeline_recall_all.sh --collection gavish         # one of them
#   ./pipeline_recall_all.sh --arms "harmony_leiden drvi_axes"
#   ./pipeline_recall_all.sh --force                     # redo the ORA and the tables
#   ./pipeline_recall_all.sh --recompute                 # drop the partition cache as well
#   ./pipeline_recall_all.sh --dry-run
#
# `gavish` is the one to read first: 22 metaprograms is the only collection of the three wide
# enough for a recall number to have room to move. `emt` has four lists inside the HVG
# background and is reported for completeness rather than as evidence.
#
# Environment: benchmark-py-r, the environment 05_3b and 05_9 run in.

set -euo pipefail

STEP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PHASE_DIR="$(cd "$STEP_DIR/.." && pwd)"
LOG_DIR="$PHASE_DIR/logs"
PYTHON="${PYTHON:-python3}"

COLLECTIONS="${COLLECTIONS:-gavish scie emt}"
ARMS="${ARMS:-harmony_leiden drvi_leiden pca_leiden drvi_axes}"
RESOLUTIONS="${RESOLUTIONS:-0.2 0.4 0.6 0.8 1.0 1.5 2.0 3.0}"
REF_RESOLUTION="${REF_RESOLUTION:-1.0}"
N_LATENT="${N_LATENT:-64}"
FORCE=0; RECOMPUTE=0; DRY_RUN=0

while [ $# -gt 0 ]; do
  case "$1" in
    --collection) COLLECTIONS="${2:?--collection needs a value}"; shift 2 ;;
    --collection=*) COLLECTIONS="${1#*=}"; shift ;;
    --arms) ARMS="${2:?--arms needs a value}"; shift 2 ;;
    --arms=*) ARMS="${1#*=}"; shift ;;
    --resolutions) RESOLUTIONS="${2:?--resolutions needs a value}"; shift 2 ;;
    --resolutions=*) RESOLUTIONS="${1#*=}"; shift ;;
    --force|-f) FORCE=1; shift ;;
    --recompute) FORCE=1; RECOMPUTE=1; shift ;;
    --dry-run|-n) DRY_RUN=1; shift ;;
    -h|--help) sed -n '2,32p' "${BASH_SOURCE[0]}"; exit 0 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done

: "${DATA_DIR:?set DATA_DIR to the datasets directory}"
mkdir -p "$LOG_DIR"
export N_LATENT

RUN_ID="$("$PYTHON" -c "
import sys; sys.path.insert(0, '$PHASE_DIR/05_2_subsetting')
import cell_set as CS; print(CS.run_id())")"
TABLE_DIR="$PHASE_DIR/tables${OUT_TAG:-}"

echo "05_10 pipeline recall"
echo "  run id      : $RUN_ID"
echo "  collections : $COLLECTIONS"
echo "  arms        : $ARMS"
echo "  resolutions : $RESOLUTIONS  (tables at $REF_RESOLUTION)"
echo "  DATA_DIR    : $DATA_DIR"
echo

for coll in $COLLECTIONS; do
  slug="$coll"; [ "$coll" = "gavish" ] && slug="gavish_tnbc"
  out="$TABLE_DIR/$slug/$RUN_ID/pipeline_recall_${slug}_${RUN_ID}.csv"
  log="$LOG_DIR/05_10_pipeline_recall_${slug}_${RUN_ID}.log"

  if [ "$FORCE" -eq 0 ] && [ -f "$out" ]; then
    echo "[have] $out"
    continue
  fi

  cmd=("$PYTHON" "$STEP_DIR/pipeline_recall_tum.py" --collection "$coll"
       --arms $ARMS --resolutions $RESOLUTIONS --ref-resolution "$REF_RESOLUTION")
  # --force redoes the step; only --recompute throws away the partitions and their DE,
  # which are the expensive half and do not depend on the collection.
  [ "$RECOMPUTE" -eq 1 ] && cmd+=(--overwrite)

  if [ "$DRY_RUN" -eq 1 ]; then
    echo "[dry] ${cmd[*]}   > $log"
    continue
  fi

  echo "[run] $coll  -> $log"
  "${cmd[@]}" 2>&1 | tee "$log"
done

echo
echo "[done] tables in $TABLE_DIR/<collection>/$RUN_ID/, figures in "
echo "       $PHASE_DIR/figures${OUT_TAG:-}/05_10_pipeline_recall/<collection>/$RUN_ID/"
