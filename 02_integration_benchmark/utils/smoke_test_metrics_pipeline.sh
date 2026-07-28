#!/usr/bin/env bash
#
# 02 utils - smoke test of the METRICS PIPELINE (local)
#
# Companion to smoke_test_metrics.py. That one validates the *stack* (scib, rpy2,
# R kBET) on an identity integration; this one validates the *pipeline glue* of
# 02_4_metrics on real integration outputs, at smoke scale: it runs the actual
# check_integrations.py / metrics.py / metrics_kbet.py / merge_metrics.py over the
# nine tiny objects in $DATA_DIR/smoke_out/ (produced by the 02_2 dispatchers),
# scoring them against $DATA_DIR/smoke_hvg.h5ad.
#
# It walks the same benchmark_grid.tsv as the real run, but LOCAL and on 5,252
# cells. The real grid runs on the cluster via 02_4_metrics/submit_metrics.slurm + 
# submit_kbet.slurm, never here.
#
# It checks the CODE PATH only: every expected metric per type must come back
# finite (13 for full, 12 for embed, 7 for the knn graph output). DRVI is a
# notebook, so it has no smoke output and is skipped; only the unscaled half is
# scored (the smoke outputs were integrated unscaled).
#
# Usage:
#   export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
#   conda activate benchmark-py-r        # so rpy2 can find R (kBET)
#   02_integration_benchmark/utils/smoke_test_metrics_pipeline.sh --dry-run
#   02_integration_benchmark/utils/smoke_test_metrics_pipeline.sh
#   ...pipeline.sh --method harmony      # one method only
#   ...pipeline.sh bbknn_unscaled        # one run_id only
#   ...pipeline.sh --no-kbet             # skip the kBET job
#
# Outputs (all disposable): $DATA_DIR/smoke_metrics/ + smoke_metrics_merged.csv.

set -euo pipefail

: "${DATA_DIR:?set DATA_DIR to the directory holding the prepared inputs}"

UTILS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
METRICS_DIR="$UTILS_DIR/../02_4_metrics"
GRID="$UTILS_DIR/../benchmark_grid.tsv"
ENV="${METRICS_ENV:-benchmark-py-r}"
TASK="shiao"
BATCH_KEY="cohort"; LABEL_KEY="cell_type"; ORGANISM="human"

ROOT="$DATA_DIR/smoke_metrics"
MERGED="$DATA_DIR/smoke_metrics_merged.csv"

[ -f "$GRID" ] || { echo "grid not found: $GRID" >&2; exit 1; }

conda_guarded() { set +u; conda "$@"; set -u; }

# Parse the command line: --method, run_id filters, --no-kbet, --no-check, --dry-run.
NO_KBET=0; NO_CHECK=0; DRY_RUN=0
METHOD_FILTER=(); FILTER=()
while [ $# -gt 0 ]; do
  case "$1" in
    --no-kbet) NO_KBET=1; shift ;;
    --no-check) NO_CHECK=1; shift ;;
    --dry-run|-n) DRY_RUN=1; shift ;;
    --method) METHOD_FILTER+=("${2:-}"); shift 2 ;;
    --method=*) METHOD_FILTER+=("${1#*=}"); shift ;;
    *) FILTER+=("$1"); shift ;;
  esac
done

# A row is wanted if it matches every active filter (method, run_id). The smoke
# only covers the unscaled half, so scaled rows are always skipped below.
want() {
  local run_id="$1" method="$2" f matched
  if [ ${#METHOD_FILTER[@]} -gt 0 ]; then
    matched=0
    for f in "${METHOD_FILTER[@]}"; do [ "$f" = "$method" ] && matched=1; done
    [ "$matched" -eq 1 ] || return 1
  fi
  [ ${#FILTER[@]} -eq 0 ] && return 0
  for f in "${FILTER[@]}"; do [ "$f" = "$run_id" ] && return 0; done
  return 1
}

if [ "$DRY_RUN" -eq 0 ]; then
  set +u; eval "$(conda shell.bash hook)"; set -u
  conda_guarded activate "$ENV"
fi

n_done=0; n_skip=0

{
  read -r _header
  while IFS=$'\t' read -r run_id method language env scaling input output types reference hvgs; do
    [ -n "$run_id" ] || continue
    [ "$scaling" = "unscaled" ] || { n_skip=$((n_skip + 1)); continue; }
    if ! want "$run_id" "$method"; then n_skip=$((n_skip + 1)); continue; fi
    if [ "$language" = "notebook" ]; then
      echo "[skip] $run_id: notebook method, no smoke output"; n_skip=$((n_skip + 1)); continue
    fi

    int_path="$DATA_DIR/smoke_out/${method}.h5ad"
    ref_path="$DATA_DIR/smoke_hvg.h5ad"
    if [ ! -f "$int_path" ]; then
      echo "[skip] $run_id: no smoke output ($int_path)"; n_skip=$((n_skip + 1)); continue
    fi

    out_dir="$ROOT/$TASK/metrics/$scaling/hvg"

    if [ "$DRY_RUN" -eq 1 ]; then
      echo "[dry] $run_id  types=$types"
      IFS=',' read -ra tps <<< "$types"
      for t in "${tps[@]}"; do echo "      $method $t -> $out_dir/${method}_${t}.csv"; done
      n_done=$((n_done + 1)); continue
    fi

    mkdir -p "$out_dir"
    echo "==================================================================="
    echo ">>> $run_id  ($method, types=$types)"
    echo "==================================================================="

    if [ "$NO_CHECK" -eq 0 ]; then
      python "$METRICS_DIR/check_integrations.py" \
        -i "$int_path" -u "$ref_path" --types "$types" -b "$BATCH_KEY" -l "$LABEL_KEY"
    fi

    IFS=',' read -ra tps <<< "$types"
    for t in "${tps[@]}"; do
      t="$(echo "$t" | tr -d '[:space:]')"
      csv="$out_dir/${method}_${t}.csv"
      echo "--- metrics: $method $t ---"
      python "$METRICS_DIR/metrics.py" \
        -u "$ref_path" -i "$int_path" -o "$csv" \
        -m "$method" --type "$t" -b "$BATCH_KEY" -l "$LABEL_KEY" \
        --organism "$ORGANISM" --hvgs "$hvgs" -v
      if [ "$NO_KBET" -eq 0 ]; then
        echo "--- kBET: $method $t ---"
        python "$METRICS_DIR/metrics_kbet.py" \
          -u "$ref_path" -i "$int_path" -o "$csv" \
          -m "$method" --type "$t" -b "$BATCH_KEY" -l "$LABEL_KEY"
      fi
    done

    n_done=$((n_done + 1))
    echo "[ok] $run_id"
  done
} < "$GRID"

if [ "$DRY_RUN" -eq 0 ]; then
  echo "==================================================================="
  echo ">>> merging per-run CSVs"
  python "$METRICS_DIR/merge_metrics.py" \
    -o "$MERGED" -r "$ROOT" --glob "$ROOT/$TASK/metrics/*/hvg/*.csv" || \
    echo "[warn] merge produced nothing (no CSVs yet?)"
  conda_guarded deactivate
fi

echo "==================================================================="
echo "done: $n_done run(s) scored, $n_skip skipped   ->  $MERGED"
