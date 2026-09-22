#!/usr/bin/env bash
#
# 05_9 embedding control: are the signatures PRESENT in the Harmony space?
#
#   1. signature_presence_tum.py --collection scie
#   2. signature_presence_tum.py --collection emt
#   3. signature_presence_tum.py --collection gavish
#
# One step, three collections, and its own driver rather than a branch of the 05_4 - 05_8
# chain. That chain is Route A and Route B on DRVI, and Route B (05_7, 05_8) reads an
# additive decoder no other method has, so wiring an embedding flag through it would put two
# steps in every run that cannot answer for a Harmony space. The DRVI chain stays as it is.
#
# WHAT THIS STEP NEEDS AND NEVER WRITES. It reads 05_6's per-cell scores, 05_4's .gmt, 05_3's
# embedding and 05_3b's - all of them read-only. It writes only its own tables and figures,
# plus one cached null per collection in $DATA_DIR/05_tum/. The per-cell score files stay
# owned by the 05_6 run, which is the point: the cells being placed are the same cells in
# every space.
#
# Usage:
#   export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
#   ./embedding_control_all.sh                          # all three collections, resuming
#   ./embedding_control_all.sh --collection gavish      # one of them
#   ./embedding_control_all.sh --spaces "harmony drvi"  # without the uncorrected PCA arm
#   ./embedding_control_all.sh --force                  # redraw the random gene sets too
#   ./embedding_control_all.sh --dry-run                # print what would run
#
# Resuming is the default and the check is existence only, as in the other drivers of the
# phase: delete a table truncated by a crash before resuming. Logs go to ../logs/.
#
# Environment: benchmark-py-r, the one holding harmony and drvi both - not because this step
# needs either (it does not; it reads h5ad files and sklearn), but because it is the
# environment 05_3b runs in and a second one is a second set of versions to report.

set -euo pipefail

STEP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PHASE_DIR="$(cd "$STEP_DIR/.." && pwd)"
LOG_DIR="$PHASE_DIR/logs"
PYTHON="${PYTHON:-python3}"

COLLECTIONS="${COLLECTIONS:-scie emt gavish}"
SPACES="${SPACES:-harmony drvi pca}"
N_LATENT="${N_LATENT:-64}"
FORCE=0; DRY_RUN=0

while [ $# -gt 0 ]; do
  case "$1" in
    --collection) COLLECTIONS="${2:?--collection needs a value}"; shift 2 ;;
    --collection=*) COLLECTIONS="${1#*=}"; shift ;;
    --spaces) SPACES="${2:?--spaces needs a value}"; shift 2 ;;
    --spaces=*) SPACES="${1#*=}"; shift ;;
    --force|-f) FORCE=1; shift ;;
    --dry-run|-n) DRY_RUN=1; shift ;;
    -h|--help) sed -n '2,40p' "${BASH_SOURCE[0]}"; exit 0 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done

: "${DATA_DIR:?set DATA_DIR to the datasets directory}"
mkdir -p "$LOG_DIR"
export N_LATENT

# The run id every output is named after, rebuilt the only way it is ever allowed to be.
RUN_ID="$("$PYTHON" -c "
import sys; sys.path.insert(0, '$PHASE_DIR/05_2_subsetting')
import cell_set as CS; print(CS.run_id())")"
TABLE_DIR="$PHASE_DIR/tables"

echo "05_9 embedding control"
echo "  run id      : $RUN_ID"
echo "  collections : $COLLECTIONS"
echo "  spaces      : $SPACES"
echo "  DATA_DIR    : $DATA_DIR"
echo

for coll in $COLLECTIONS; do
  # `gavish` writes under the `gavish_tnbc` slug unless --all-metaprograms is given; the
  # script resolves that and this only has to know where to look for the resume check.
  slug="$coll"; [ "$coll" = "gavish" ] && slug="gavish_tnbc"
  out="$TABLE_DIR/$slug/$RUN_ID/signature_presence_${slug}_${RUN_ID}.csv"
  log="$LOG_DIR/05_9_signature_presence_${slug}_${RUN_ID}.log"

  if [ "$FORCE" -eq 0 ] && [ -f "$out" ]; then
    echo "[have] $out"
    continue
  fi

  cmd=("$PYTHON" "$STEP_DIR/signature_presence_tum.py" --collection "$coll" --spaces $SPACES)
  [ "$FORCE" -eq 1 ] && cmd+=(--overwrite)

  if [ "$DRY_RUN" -eq 1 ]; then
    echo "[dry] ${cmd[*]}   > $log"
    continue
  fi

  echo "[run] $coll  -> $log"
  "${cmd[@]}" 2>&1 | tee "$log"
done

echo
echo "[done] tables in $TABLE_DIR/<collection>/$RUN_ID/, figures in "
echo "       $PHASE_DIR/figures/05_9_embedding_control/<collection>/$RUN_ID/"
