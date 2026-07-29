#!/usr/bin/env bash
#
# 02_4 metrics: the metrics step, runnable locally OR on SLURM.
#
# Walks benchmark_grid.tsv and, for every run whose integrated object exists,
# scores each of its output type(s) with the 13 scib metrics (12 + kBET). One row
# with types "full,embed" becomes two jobs, so 17 runs expand into 21.
#
#   default (local): runs check_integrations.py + metrics.py + metrics_kbet.py for
#                    each (run, type) here, in sequence, then merges the CSVs.
#   --slurm        : submits submit_metrics.slurm + submit_kbet.slurm restricted to
#                    the matching array indices, throttled to 3 concurrent tasks
#                    (at most 3 nodes in `normal`), kBET afterok-dependent on the
#                    metrics array. Merge + summary are run by hand once both
#                    arrays finish (the command is printed at the end).
#
# CSVs go to $DATA_DIR/02_metrics/<task>/metrics/<scaling>/hvg/<method>_<type>.csv.
#
# Usage:
#   export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
#   ./run_all_metrics.sh                       # local, whole grid
#   ./run_all_metrics.sh --slurm               # submit the whole grid (3 at a time)
#   ./run_all_metrics.sh --method harmony      # one method, both scalings
#   ./run_all_metrics.sh --slurm --method harmony
#   ./run_all_metrics.sh --scaling unscaled    # the unscaled half
#   ./run_all_metrics.sh bbknn_unscaled ...    # named run_id(s)
#   ./run_all_metrics.sh --no-kbet             # skip the kBET job/array
#   ./run_all_metrics.sh --dry-run             # print what would run, do nothing
#   ./run_all_metrics.sh --force               # re-score tasks whose CSV is complete
#
# Resuming is the default, as in run_all.sh and plot_all.sh: a (run, type) task
# whose CSV already holds a kBET value is reported [have] and skipped, so
# re-running the command after a crash picks up where it stopped. The kBET row is
# what makes the CSV complete - metrics.py writes the file, metrics_kbet.py fills
# that row afterwards - so a pair killed between the two is correctly re-run.
# --force re-scores.
#
# --scaling, --method and run_id filters combine (a task must match all active).

set -euo pipefail

: "${DATA_DIR:?set DATA_DIR to the directory holding the prepared inputs and integration outputs}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GRID="$SCRIPT_DIR/../benchmark_grid.tsv"
TASK="shiao"
ENV="${METRICS_ENV:-benchmark-py-r}"    # local env; the .slurm files use catalano_env
BATCH_KEY="cohort"; LABEL_KEY="cell_type"; ORGANISM="human"
ROOT="$DATA_DIR/02_metrics"
MERGED="$DATA_DIR/02_metrics_merged.csv"
MAX_CONCURRENT=3                        # SLURM: at most 3 nodes in `normal`

[ -f "$GRID" ] || { echo "grid not found: $GRID" >&2; exit 1; }

conda_guarded() { set +u; conda "$@"; set -u; }
in_list() { local x="$1"; shift; local e; for e in "$@"; do [ "$e" = "$x" ] && return 0; done; return 1; }

# Parse the command line.
USE_SLURM=0; NO_KBET=0; NO_CHECK=0; DRY_RUN=0; FORCE=0
SCALING_FILTER=""
METHOD_FILTER=(); RUNID_FILTER=()
while [ $# -gt 0 ]; do
  case "$1" in
    --slurm) USE_SLURM=1; shift ;;
    --no-kbet) NO_KBET=1; shift ;;
    --no-check) NO_CHECK=1; shift ;;
    --force|-f) FORCE=1; shift ;;
    --dry-run|-n) DRY_RUN=1; shift ;;
    --scaling) SCALING_FILTER="${2:-}"; shift 2 ;;
    --scaling=*) SCALING_FILTER="${1#*=}"; shift ;;
    --method) METHOD_FILTER+=("${2:-}"); shift 2 ;;
    --method=*) METHOD_FILTER+=("${1#*=}"); shift ;;
    *) RUNID_FILTER+=("$1"); shift ;;
  esac
done
if [ -n "$SCALING_FILTER" ] && [ "$SCALING_FILTER" != "scaled" ] && [ "$SCALING_FILTER" != "unscaled" ]; then
  echo "--scaling must be 'scaled' or 'unscaled', got '$SCALING_FILTER'" >&2; exit 1
fi

# A grid row (by run_id/scaling/method) is wanted if it matches every filter.
want_row() {
  local run_id="$1" scaling="$2" method="$3"
  [ -n "$SCALING_FILTER" ] && [ "$scaling" != "$SCALING_FILTER" ] && return 1
  [ ${#METHOD_FILTER[@]} -gt 0 ] && ! in_list "$method" "${METHOD_FILTER[@]}" && return 1
  [ ${#RUNID_FILTER[@]} -gt 0 ] && ! in_list "$run_id" "${RUNID_FILTER[@]}" && return 1
  return 0
}

# Where a (run, type) pair's CSV lands. One place, so the local loop and the
# SLURM index filter can never disagree about what "already scored" means.
csv_path() { echo "$ROOT/$TASK/metrics/$2/hvg/${1}_${3}.csv"; }

# Already scored? The CSV existing is not enough: metrics.py writes it, then
# metrics_kbet.py patches the kBET row into that same file, so a pair killed
# between the two leaves a complete-looking CSV with kBET still empty. Require a
# non-empty kBET value unless --no-kbet, which is the only case where a CSV
# without it is finished work.
#
# (metrics_kbet.py also writes NaN deliberately when --max-cells is exceeded;
# this driver never passes it and the default is 0 = no cap, so that cannot
# happen here. It would show up as a pair that never stops being re-run.)
already_done() {
  local csv="$1" v
  [ "$FORCE" -eq 0 ] || return 1
  [ -f "$csv" ] || return 1
  [ "$NO_KBET" -eq 1 ] && return 0
  v="$(awk -F, 'NR>1 && $1=="kBET" { gsub(/[ \t"]/, "", $2); print $2; exit }' "$csv")"
  case "$v" in ""|nan|NaN|NA|None) return 1 ;; *) return 0 ;; esac
}

# SLURM mode: resolve the matching array indices and submit the arrays.

run_slurm() {
  # Expand the grid exactly as the .slurm scripts do (index = SLURM_ARRAY_TASK_ID).
  # Scored pairs are dropped from the array spec rather than submitted and skipped
  # inside the job: a queued task that exits immediately still costs a slot
  # against the %3 throttle. Both arrays share the spec, so a pair is resubmitted
  # whole when its kBET is missing - metrics.py is the cheap half of the two.
  local indices=() total=0 n_have=0
  while IFS=$'\t' read -r idx run_id method scaling type; do
    total=$((total + 1))
    want_row "$run_id" "$scaling" "$method" || continue
    if already_done "$(csv_path "$method" "$scaling" "$type")"; then
      echo "[have] $run_id $type: already scored, not submitting"
      n_have=$((n_have + 1)); continue
    fi
    indices+=("$idx")
  done < <(awk -F'\t' 'NR>1 && $1!="" { n=split($8, ts, ","); for (i=1;i<=n;i++){ idx++; gsub(/[ \t]/,"",ts[i]); print idx"\t"$1"\t"$2"\t"$5"\t"ts[i] } }' "$GRID")

  if [ ${#indices[@]} -eq 0 ]; then
    if [ "$n_have" -gt 0 ]; then
      echo "nothing to submit: all $n_have matching task(s) are already scored (--force to re-run)"
      return 0
    fi
    echo "no array task matches the given filters" >&2; exit 1
  fi

  local spec; spec="$(IFS=,; echo "${indices[*]}")%${MAX_CONCURRENT}"
  echo "selected ${#indices[@]}/${total} array task(s), ${n_have} already scored, throttled to ${MAX_CONCURRENT}: $spec"

  if [ "$DRY_RUN" -eq 1 ]; then
    echo "[dry] sbatch --export=ALL,DATA_DIR=$DATA_DIR,METRICS_DIR=$SCRIPT_DIR --array=$spec $SCRIPT_DIR/submit_metrics.slurm"
    [ "$NO_KBET" -eq 0 ] && echo "[dry] sbatch --dependency=afterok:<jobid> --export=ALL,DATA_DIR=$DATA_DIR,METRICS_DIR=$SCRIPT_DIR --array=$spec $SCRIPT_DIR/submit_kbet.slurm"
    return 0
  fi

  local m; m="$(sbatch --parsable --export="ALL,DATA_DIR=$DATA_DIR,METRICS_DIR=$SCRIPT_DIR" --array="$spec" "$SCRIPT_DIR/submit_metrics.slurm")"
  echo "submitted metrics array: job $m"
  if [ "$NO_KBET" -eq 0 ]; then
    local k; k="$(sbatch --parsable --dependency=afterok:"$m" --export="ALL,DATA_DIR=$DATA_DIR,METRICS_DIR=$SCRIPT_DIR" --array="$spec" "$SCRIPT_DIR/submit_kbet.slurm")"
    echo "submitted kBET array:    job $k  (afterok:$m)"
  else
    echo "kBET array skipped (--no-kbet)"
  fi
  echo "when both arrays finish, merge + plot locally:"
  echo "  python $SCRIPT_DIR/merge_metrics.py -o $MERGED -r $ROOT --glob '$ROOT/$TASK/metrics/*/hvg/*.csv'"
  echo "  Rscript $SCRIPT_DIR/make_summary_table.R -i $MERGED -o figures"
}

# Local mode: run every matching (run, type) here, then merge.

run_local() {
  if [ "$DRY_RUN" -eq 0 ]; then
    set +u; eval "$(conda shell.bash hook)"; set -u
    conda_guarded activate "$ENV"
  fi

  local n_done=0 n_skip=0 n_have=0
  {
    read -r _header
    while IFS=$'\t' read -r run_id method language env scaling input output types reference hvgs; do
      [ -n "$run_id" ] || continue
      if ! want_row "$run_id" "$scaling" "$method"; then n_skip=$((n_skip + 1)); continue; fi

      local int_path="$DATA_DIR/$output" ref_path="$DATA_DIR/$reference"
      local out_dir="$ROOT/$TASK/metrics/$scaling/hvg"

      # Scoring is per (run, type): a "full,embed" row is two independent jobs, so
      # a row half-scored resumes on the half that is missing.
      IFS=',' read -ra tps <<< "$types"
      local todo=() t_clean
      for t in "${tps[@]}"; do
        t_clean="$(echo "$t" | tr -d '[:space:]')"
        if already_done "$(csv_path "$method" "$scaling" "$t_clean")"; then
          echo "[have] $run_id $t_clean: already scored ($out_dir/${method}_${t_clean}.csv)"
          n_have=$((n_have + 1))
        else
          todo+=("$t_clean")
        fi
      done
      [ ${#todo[@]} -eq 0 ] && continue

      # Missing-input check after the done check: whether a task's CSV is complete
      # does not depend on its input still being readable, and reporting a scored
      # task as "not integrated yet" would be a false alarm.
      if [ ! -f "$int_path" ]; then
        echo "[skip] $run_id: not integrated yet ($int_path)"; n_skip=$((n_skip + 1)); continue
      fi

      if [ "$DRY_RUN" -eq 1 ]; then
        echo "[dry] $run_id  types=${todo[*]}  ref=$ref_path"
        for t in "${todo[@]}"; do echo "      $method $t -> $out_dir/${method}_${t}.csv"; done
        n_done=$((n_done + 1)); continue
      fi

      mkdir -p "$out_dir"
      echo "==================================================================="
      echo ">>> $run_id  ($method, types=$types)"
      echo "==================================================================="
      if [ "$NO_CHECK" -eq 0 ]; then
        python "$SCRIPT_DIR/check_integrations.py" -i "$int_path" -u "$ref_path" \
          --types "$types" -b "$BATCH_KEY" -l "$LABEL_KEY"
      fi
      for t in "${todo[@]}"; do
        local csv="$out_dir/${method}_${t}.csv"
        echo "--- metrics: $method $t ---"
        python "$SCRIPT_DIR/metrics.py" -u "$ref_path" -i "$int_path" -o "$csv" \
          -m "$method" --type "$t" -b "$BATCH_KEY" -l "$LABEL_KEY" \
          --organism "$ORGANISM" --hvgs "$hvgs" -v
        if [ "$NO_KBET" -eq 0 ]; then
          echo "--- kBET: $method $t ---"
          python "$SCRIPT_DIR/metrics_kbet.py" -u "$ref_path" -i "$int_path" -o "$csv" \
            -m "$method" --type "$t" -b "$BATCH_KEY" -l "$LABEL_KEY"
        fi
      done
      n_done=$((n_done + 1)); echo "[ok] $run_id"
    done
  } < "$GRID"

  if [ "$DRY_RUN" -eq 0 ]; then
    echo ">>> merging per-run CSVs"
    python "$SCRIPT_DIR/merge_metrics.py" -o "$MERGED" -r "$ROOT" \
      --glob "$ROOT/$TASK/metrics/*/hvg/*.csv" || echo "[warn] merge produced nothing"
    conda_guarded deactivate
  fi
  echo "==================================================================="
  echo "done: $n_done run(s) scored, $n_have task(s) already scored, $n_skip skipped   ->  $MERGED"
  [ "$n_have" -gt 0 ] && [ "$FORCE" -eq 0 ] && \
    echo "      (--force to re-score the $n_have existing task(s))"
  return 0
}


if [ "$USE_SLURM" -eq 1 ]; then run_slurm; else run_local; fi
