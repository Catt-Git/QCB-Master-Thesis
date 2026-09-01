#!/usr/bin/env bash
#
# 04_9 embedding control: Route A on a coordinate system other than DRVI's, end to end.
#
#   1. run_harmony_epi.py          shiao_epi_hvg_2k.h5ad -> embed_harmony_epi_50.h5ad
#                                                        (+ embed_pca_epi_50.h5ad, the null arm)
#   2. cell_first_epi.py           the SAME Route A, --embedding <name>, per collection
#   3. compare_embeddings_epi.py   the two runs side by side -> embedding_comparison table
#
# Its own driver rather than a branch of signature_interpretation_all.sh, and deliberately:
# that chain is 04_3 -> 04_7 on DRVI and includes Route B and Route C, which are DRVI-only
# (they read the additive decoder off embed.varm and no method here has one). Wiring an
# embedding flag through it would put two steps in every run that cannot answer for a
# Harmony space. The DRVI chain therefore stays exactly as it is.
#
# WHAT STEP 2 DOES AND DOES NOT WRITE. A1 - A5 of Route A - the scoring, the within-stratum
# standardisation, the target region, the consensus vote and the confounder checks - are
# computed from shiao_epi.h5ad and never see an embedding. A control run recomputes them
# identically and writes only the three that depend on the space:
# dim_signature_spearman, dim_target_effect_size, dimension_row_order, plus its own heatmap.
# The cell-level tables stay owned by the drvi_epi_64 run, which is the point: the cells
# being placed are the same cells in both spaces.
#
# Usage:
#   export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
#   ./embedding_control_all.sh                        # harmony, both collections, resuming
#   ./embedding_control_all.sh --collection scie      # one collection
#   ./embedding_control_all.sh --embeddings "harmony pca"   # with the uncorrected null arm
#   ./embedding_control_all.sh --force                # re-run everything
#   ./embedding_control_all.sh --dry-run              # print what would run
#   ./embedding_control_all.sh cellfirst compare      # only the named steps
#
# Step names: harmony, cellfirst, compare. Resuming is the default: a step whose output
# exists is reported [have] and skipped, as in the other drivers of the phase; the check is
# existence only, so delete a file truncated by a crash before resuming.
# Environment benchmark-py-r (harmony-pytorch 0.1.8 is already in it, from phase 02); the
# driver does not activate it. Logs go to 04_drvi_epithelial/logs/.

set -euo pipefail

STEP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PHASE_DIR="$(cd "$STEP_DIR/.." && pwd)"
LOG_DIR="$PHASE_DIR/logs"
PYTHON="${PYTHON:-python3}"

COLLECTIONS="${COLLECTIONS:-scie emt}"
EMBEDDINGS="${EMBEDDINGS:-harmony}"
FORCE=0; DRY_RUN=0
STEP_FILTER=()

while [ $# -gt 0 ]; do
  case "$1" in
    --collection) COLLECTIONS="${2:?--collection needs a value}"; shift 2 ;;
    --collection=*) COLLECTIONS="${1#*=}"; shift ;;
    --embeddings) EMBEDDINGS="${2:?--embeddings needs a value}"; shift 2 ;;
    --embeddings=*) EMBEDDINGS="${1#*=}"; shift ;;
    --force|-f) FORCE=1; shift ;;
    --dry-run|-n) DRY_RUN=1; shift ;;
    -h|--help) sed -n '2,40p' "${BASH_SOURCE[0]}"; exit 0 ;;
    -*) echo "unknown option: $1" >&2; exit 1 ;;
    *) STEP_FILTER+=("$1"); shift ;;
  esac
done

: "${DATA_DIR:?set DATA_DIR to the directory holding the datasets (outside the repo)}"
EPI_DIR="$DATA_DIR/04_epi"

in_list() { local x="$1"; shift; local e; for e in "$@"; do [ "$e" = "$x" ] && return 0; done; return 1; }
want() { [ "${#STEP_FILTER[@]}" -eq 0 ] || in_list "$1" "${STEP_FILTER[@]}"; }

for s in "${STEP_FILTER[@]:-}"; do
  [ -z "$s" ] && continue
  in_list "$s" harmony cellfirst compare || { echo "unknown step: $s (have: harmony, cellfirst, compare)" >&2; exit 1; }
done
for c in $COLLECTIONS; do
  case "$c" in scie|emt) ;; *) echo "unknown collection: $c (have: scie, emt)" >&2; exit 1 ;; esac
done
for e in $EMBEDDINGS; do
  case "$e" in harmony|pca) ;; *) echo "unknown embedding: $e (have: harmony, pca; drvi is the reference and is not re-run here)" >&2; exit 1 ;; esac
done

# The run id of an embedding, which is what every output is named with. Kept in step with
# the EMBEDDINGS registry of utils/signature_common.py - if LINEAR_N_DIMS moves there, the
# two numbers below move with it.
run_id_of() {
  case "$1" in
    harmony) echo "harmony_epi_50" ;;
    pca) echo "pca_epi_50" ;;
    drvi) echo "drvi_epi_64" ;;
  esac
}

mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/embedding_control_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "04_9 embedding control"
echo "collections $COLLECTIONS"
echo "embeddings  $EMBEDDINGS  (compared against drvi_epi_64, the phase's own run)"
echo "DATA_DIR    $DATA_DIR"
echo "log         $LOG_FILE"
echo

run() {  # run <name> <dir> <script> <output-that-marks-it-done> [args...]
  local name="$1" dir="$2" script="$3" out="$4"; shift 4
  if [ "$FORCE" -eq 0 ] && [ -e "$out" ]; then
    echo "[have] $name  -> $out"
    return 0
  fi
  if [ "$DRY_RUN" -eq 1 ]; then
    echo "[dry ] $name  -> $PYTHON $script $*"
    return 0
  fi
  echo "[run ] $name  -> $script $*"
  ( cd "$dir" && "$PYTHON" "$script" "$@" )
  echo "[done] $name"
  echo
}

# 1. the embeddings. One call: run_harmony_epi.py writes Harmony and, unless told not to,
#    the uncorrected PCA in the same pass - they share the PCA, so splitting them would
#    compute it twice and risk two different ones.
if want harmony; then
  pca_flag=()
  in_list pca $EMBEDDINGS || pca_flag=(--no-pca-arm)
  force_flag=()
  [ "$FORCE" -eq 1 ] && force_flag=(--overwrite)
  marker="$EPI_DIR/embed_$(run_id_of harmony).h5ad"
  # If the PCA arm was asked for and is not on disk, the step is not done even when Harmony
  # is: point the marker at the missing file. run_harmony_epi.py then skips whichever of the
  # two is already there, so nothing is recomputed twice.
  if in_list pca $EMBEDDINGS && [ ! -e "$EPI_DIR/embed_$(run_id_of pca).h5ad" ]; then
    marker="$EPI_DIR/embed_$(run_id_of pca).h5ad"
  fi
  run harmony "$STEP_DIR" run_harmony_epi.py "$marker" "${pca_flag[@]}" "${force_flag[@]}"
fi

# 2. Route A on each of them, per collection.
if want cellfirst; then
  for e in $EMBEDDINGS; do
    for c in $COLLECTIONS; do
      run "cellfirst $e/$c" "$PHASE_DIR/04_5_cell_first" cell_first_epi.py \
        "$PHASE_DIR/tables/$c/dim_signature_spearman_${c}_$(run_id_of "$e").csv" \
        --collection "$c" --embedding "$e"
    done
  done
fi

# 3. the comparison, one per collection, always against the DRVI run.
if want compare; then
  for c in $COLLECTIONS; do
    run "compare $c" "$STEP_DIR" compare_embeddings_epi.py \
      "$PHASE_DIR/tables/$c/embedding_comparison_${c}_epi_embeddings.csv" \
      --collection "$c" --embeddings drvi $EMBEDDINGS
  done
fi

echo "all requested steps finished."
