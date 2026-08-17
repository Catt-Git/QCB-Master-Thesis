#!/usr/bin/env bash
#
# 02_4 metrics: the metrics step, runnable locally OR on SLURM.
#
# Walks benchmark_grid.tsv and, for every run whose integrated object exists,
# scores each of its output type(s) with the 12 scib metrics. One row with types
# "full,embed" becomes two jobs, so the 19 rows expand into 23.
#
# The two `unintegrated_*` rows are the baseline and are scored like any other:
# their `output` and `reference` are the same prepared object, so metrics.py ends
# up comparing it with itself, which is exactly what the Unintegrated row of the
# scIB tables is. Unlike 02_2 and 02_3, this step needs no special case for them.
#
#   default (local): runs check_integrations.py + metrics.py for each (run, type)
#                    here, in sequence, then merges the CSVs.
#   --slurm        : submits submit_metrics.slurm restricted to the matching array
#                    indices, 5 concurrent tasks. The indices are split into two
#                    groups by output type - full/embed at MEM_STD, knn at MEM_KNN
#                    - so each is submitted with the memory it actually needs; that
#                    is two arrays, two job ids. Merge + summary are run by hand
#                    once they finish (the command is printed at the end).
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
#   ./run_all_metrics.sh --dry-run             # print what would run, do nothing
#   ./run_all_metrics.sh --force               # re-score tasks whose CSV exists
#
# Resuming is the default, as in run_all.sh and plot_all.sh: a (run, type) task
# whose CSV exists is reported [have] and skipped, so re-running the command after
# a crash picks up where it stopped. metrics.py writes that CSV in one go, at the
# end, so its existence means the task finished. --force re-scores.
#
# --scaling, --method and run_id filters combine (a task must match all active).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GRID="$SCRIPT_DIR/../benchmark_grid.tsv"

# DATA_DIR: an explicit one wins, otherwise the cluster data root when running
# there. That root is the same one sync_to_cluster.sh uploads into, defined once in
# utils/cluster_env.sh - the two used to disagree, silently.
. "$SCRIPT_DIR/../utils/cluster_env.sh"
resolve_data_dir || exit 1
TASK="shiao"
ENV="${METRICS_ENV:-benchmark-py-r}"    # local env; the .slurm files use catalano_env
BATCH_KEY="cohort"; LABEL_KEY="cell_type"; ORGANISM="human"
ROOT="$DATA_DIR/02_metrics"
MERGED="$DATA_DIR/02_metrics_merged.csv"

# SLURM sizing. `normal` has six nodes and a task takes 4 CPUs, not a whole node
# (two tasks shared node08 on the first run), so five at a time is comfortable.
MAX_CONCURRENT=5

# Memory is per output type, not per run. full/embed load a dense corrected matrix
# (5-10 GB), rebuild a 15-neighbour graph and peak somewhere under 64G - 96G is that
# with margin. knn (BBKNN) is a different animal: its graph carries 8.7e8 nonzeros,
# scanpy's leiden materialises it as ~1e9 Python tuples, and the job was measured at
# 197 GiB before finishing the clustering. 300G, which also confines it to node04.
MEM_STD="${METRICS_MEM:-96G}"
MEM_KNN="${METRICS_MEM_KNN:-300G}"

# Nodes to keep away from. node02 was accepting jobs and failing them at launch on
# 2026-07-30 (ExitCode 0:53, no log file written at all) while `sinfo -R` reported
# nothing wrong; six tasks died there in fourteen seconds. Clear METRICS_EXCLUDE
# once it is fixed.
EXCLUDE_NODES="${METRICS_EXCLUDE-node02}"

# leiden implementation for the NMI/ARI clustering, passed through to metrics.py.
# Benchmark-wide: mixing flavors makes the NMI/ARI rows incomparable across methods.
CLUSTER_FLAVOR="${CLUSTER_FLAVOR:-igraph}"

[ -f "$GRID" ] || { echo "grid not found: $GRID" >&2; exit 1; }

conda_guarded() { set +u; conda "$@"; set -u; }
in_list() { local x="$1"; shift; local e; for e in "$@"; do [ "$e" = "$x" ] && return 0; done; return 1; }

# Parse the command line.
USE_SLURM=0; NO_CHECK=0; DRY_RUN=0; FORCE=0
SCALING_FILTER=""
METHOD_FILTER=(); RUNID_FILTER=()
while [ $# -gt 0 ]; do
  case "$1" in
    --slurm) USE_SLURM=1; shift ;;
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

# Already scored? metrics.py writes the CSV once, at the end, from the DataFrame
# scib returns, so the file existing means that (run, type) ran to completion; a
# task killed halfway leaves no file at all.
already_done() {
  local csv="$1"
  [ "$FORCE" -eq 0 ] || return 1
  [ -f "$csv" ] || return 1
  return 0
}

# SLURM mode: resolve the matching array indices and submit the arrays.

# Submit one group of array indices. A group exists because `--mem` is a
# submission-level option and the knn tasks need four times what the others do;
# everything else about the two submissions is identical.
submit_group() {
  local label="$1" mem="$2"; shift 2
  [ $# -gt 0 ] || return 0

  local spec; spec="$(IFS=,; echo "$*")%${MAX_CONCURRENT}"
  local opts=(--export="ALL,DATA_DIR=$DATA_DIR,METRICS_DIR=$SCRIPT_DIR,CLUSTER_FLAVOR=$CLUSTER_FLAVOR"
              --array="$spec" --mem="$mem")
  [ -n "$EXCLUDE_NODES" ] && opts+=(--exclude="$EXCLUDE_NODES")

  if [ "$DRY_RUN" -eq 1 ]; then
    echo "[dry] $label: sbatch ${opts[*]} $SCRIPT_DIR/submit_metrics.slurm"
    return 0
  fi

  local m; m="$(sbatch --parsable "${opts[@]}" "$SCRIPT_DIR/submit_metrics.slurm")"
  echo "submitted $label metrics array: job $m   (--mem=$mem, tasks $spec)"
}

run_slurm() {
  # Expand the grid exactly as the .slurm scripts do (index = SLURM_ARRAY_TASK_ID).
  # Scored pairs are dropped from the array spec rather than submitted and skipped
  # inside the job: a queued task that exits immediately still costs a slot
  # against the throttle.
  #
  # The indices are split by output type because the two groups are submitted with
  # different --mem (see MEM_STD / MEM_KNN).
  local indices=() indices_knn=() total=0 n_have=0 n_skip=0
  while IFS=$'\t' read -r idx run_id method scaling type output reference; do
    total=$((total + 1))
    want_row "$run_id" "$scaling" "$method" || continue
    if already_done "$(csv_path "$method" "$scaling" "$type")"; then
      echo "[have] $run_id $type: already scored, not submitting"
      n_have=$((n_have + 1)); continue
    fi
    # Same check run_local makes, on both files the job reads. Submitting a task
    # whose input is missing costs a queue slot to fail two seconds later, and
    # when integration runs locally and metrics on the cluster (sync_to_cluster.sh
    # in between) a whole half of the grid is routinely still in transit. The
    # reference is checked too: it travels separately from the integrated object,
    # and a scaled run scored against a missing scaled reference fails just as
    # late and far less obviously.
    local missing=""
    [ -f "$DATA_DIR/$output" ]    || missing="$DATA_DIR/$output"
    [ -f "$DATA_DIR/$reference" ] || missing="${missing:+$missing, }$DATA_DIR/$reference"
    if [ -n "$missing" ]; then
      echo "[skip] $run_id $type: input not here yet ($missing)"
      n_skip=$((n_skip + 1)); continue
    fi
    if [ "$type" = "knn" ]; then indices_knn+=("$idx"); else indices+=("$idx"); fi
  done < <(awk -F'\t' 'NR>1 && $1!="" { n=split($8, ts, ","); for (i=1;i<=n;i++){ idx++; gsub(/[ \t]/,"",ts[i]); print idx"\t"$1"\t"$2"\t"$5"\t"ts[i]"\t"$7"\t"$9 } }' "$GRID")

  local n_sel=$(( ${#indices[@]} + ${#indices_knn[@]} ))
  if [ "$n_sel" -eq 0 ]; then
    if [ "$n_have" -gt 0 ]; then
      echo "nothing to submit: all $n_have matching task(s) are already scored (--force to re-run)"
      return 0
    fi
    echo "no array task matches the given filters" >&2; exit 1
  fi

  echo "selected ${n_sel}/${total} array task(s), ${n_have} already scored, ${n_skip} not integrated yet, ${MAX_CONCURRENT} at a time"
  echo "  clustering: leiden ($CLUSTER_FLAVOR)${EXCLUDE_NODES:+   excluding: $EXCLUDE_NODES}"
  submit_group "full/embed" "$MEM_STD" ${indices[@]+"${indices[@]}"}
  submit_group "knn       " "$MEM_KNN" ${indices_knn[@]+"${indices_knn[@]}"}

  echo "when the arrays finish, merge + plot locally:"
  echo "  python $SCRIPT_DIR/merge_metrics.py -o $MERGED -r $ROOT --glob '$ROOT/$TASK/metrics/*/hvg/*.csv'"
  echo "  Rscript $SCRIPT_DIR/make_summary_table.R -i $MERGED"
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
          --organism "$ORGANISM" --hvgs "$hvgs" --cluster-flavor "$CLUSTER_FLAVOR" -v
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
