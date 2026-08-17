#!/usr/bin/env bash
#
# 02_3 method plots: driver.
#
# Similarly to 02_2's run_all.sh, but for the QC panels. Reads benchmark_grid.tsv 
# and runs plot_methods_umaps.py for every run whose integrated .h5ad already exists,
# writing figures/<run_id>/ (i.e. figures/<method>_<scaling>/). Rows not yet
# integrated are skipped, so it is safe to run after a partial run_all.sh.
#
# The plotting --type is derived from the grid's `types` column: the corrected
# low-dimensional embedding is the best QC view, so `embed` wins when a run
# produces one (scanorama, fastmnn are full,embed); otherwise `full`, else `knn`.
#
# Reads the integrated object from the grid's `output` column and the matching
# unintegrated reference from its `reference` column (which already tracks the
# scaling variant), so a scaled run is compared against the scaled reference: the
# "before" picture is the input the method actually received, z-scoring included.
#
# Both references must therefore carry obsm['X_umap']: the unscaled one gets it
# from 01_6, the scaled one from 02_1_prepare/scale_batch.py, which rebuilds it on
# the scaled matrix after dropping the stale 01_6 layout. A scaled object produced
# before that step existed has no UMAP, and every scaled row aborts on the
# plotter's assertion.
#
# Runnable locally OR on SLURM, like 02_4's run_all_metrics.sh:
#
#   default (local): plots each matching run here, in sequence. A run that fails
#                    is reported and the batch continues; the failures are listed
#                    at the end and the exit code is non-zero.
#   --slurm        : submits submit_plots.slurm restricted to the matching array
#                    indices, MAX_CONCURRENT at a time. The indices are split by
#                    output type - knn at MEM_KNN, the rest at MEM_STD - so each
#                    group is submitted with the memory it actually needs.
#
# Usage:
#   export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
#   ./plot_all.sh                       # every integrated run
#   ./plot_all.sh harmony_unscaled ...  # only the named run_id(s)
#   ./plot_all.sh --scaling unscaled    # only the unscaled half
#   ./plot_all.sh --dry-run             # print what would run, do nothing
#   ./plot_all.sh --force               # redraw runs whose panels already exist
#   ./plot_all.sh --force --recompute-umap   # redraw *and* lay the UMAPs out again
#   ./plot_all.sh --slurm --force       # the whole grid as an array job
#
# On the cluster the panels land in the repo's figures/, NOT under DATA_DIR (see
# FIG_DIR below), so they have to be rsynced back; the command is printed at the
# end of a --slurm submission.
#
# The integrated layout is cached into each integrated .h5ad the first time it is
# computed (obsm['X_umap'] + uns['integrated_umap']), so a redraw costs seconds
# rather than the ~25 min a 620k-cell UMAP takes. --recompute-umap discards that
# cache; without it a --force redraw reuses the stored layout, which is what you
# want when only colours, sizes or panel layout changed.
#
# Resuming is the default, as in run_all.sh: a run whose five panels are all in
# figures/<run_id>/ is reported [have] and skipped, so re-running the command
# after a crash picks up where it stopped. --force redraws.
#
# --scaling and run_id filters combine (a row must match both), exactly as in
# run_all.sh.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PHASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
GRID="$PHASE_DIR/benchmark_grid.tsv"
FIG_DIR="$PHASE_DIR/figures"          # the phase's one figures/, not a local one
PLOTTER="$SCRIPT_DIR/plot_methods_umaps.py"

# DATA_DIR: an explicit one wins, otherwise the cluster data root when running
# there - the same resolution 02_4's wrappers use, defined once in cluster_env.sh.
. "$SCRIPT_DIR/../utils/cluster_env.sh"
resolve_data_dir || exit 1

ENV="${PLOT_ENV:-benchmark-py-r}"       # local env; the .slurm file uses catalano_env

# SLURM sizing, per output type, as in 02_4. A `full` run rebuilds a PCA and a
# 15-neighbour graph on 620k x 2k; `embed` only the graph. The two `knn` rows are
# the outlier: BBKNN's own graph carries 8.7e8 nonzeros and sc.tl.umap lays that
# out directly, so they get their own group with the memory the 02_4 knn tasks
# needed.
MAX_CONCURRENT="${PLOT_MAX_CONCURRENT:-5}"
MEM_STD="${PLOT_MEM:-96G}"
MEM_KNN="${PLOT_MEM_KNN:-300G}"

# node02 accepted and killed jobs at launch on 2026-07-30 (ExitCode 0:53); same
# exclusion as 02_4. Clear PLOT_EXCLUDE once it is fixed.
EXCLUDE_NODES="${PLOT_EXCLUDE-node02}"

[ -f "$GRID" ] || { echo "grid not found: $GRID" >&2; exit 1; }

# conda's env activation scripts (e.g. MKL's) read unbound variables, which
# aborts under `set -u`; run conda with nounset temporarily disabled.
conda_guarded() { set +u; conda "$@"; set -u; }

# Parse the command line: --scaling, --dry-run, --force, plus optional run_id filters.
SCALING_FILTER=""
DRY_RUN=0
FORCE=0
RECOMPUTE_UMAP=0
USE_SLURM=0
FILTER=()
while [ $# -gt 0 ]; do
  case "$1" in
    --scaling) SCALING_FILTER="${2:-}"; shift 2 ;;
    --scaling=*) SCALING_FILTER="${1#*=}"; shift ;;
    --force|-f) FORCE=1; shift ;;
    --dry-run|-n) DRY_RUN=1; shift ;;
    --recompute-umap) RECOMPUTE_UMAP=1; shift ;;
    --slurm) USE_SLURM=1; shift ;;
    # A mistyped option must not become a run_id filter: it would match no row,
    # plot nothing and still exit 0 - a typo indistinguishable from success.
    # Same guard as run_all.sh.
    -*) echo "unknown option '$1'" >&2; exit 1 ;;
    *) FILTER+=("$1"); shift ;;
  esac
done
if [ -n "$SCALING_FILTER" ] && [ "$SCALING_FILTER" != "scaled" ] && [ "$SCALING_FILTER" != "unscaled" ]; then
  echo "--scaling must be 'scaled' or 'unscaled', got '$SCALING_FILTER'" >&2
  exit 1
fi

# A row is wanted if it matches every active filter (run_id and scaling).
#
# `baseline` rows (Unintegrated) are never wanted: their `output` and `reference`
# are the same file, so all five panels would compare the object with itself. The
# unintegrated space already has its own figures, from 01_6.
want() {
  local run_id="$1" scaling="$2" language="$3" f
  [ "$language" = "baseline" ] && return 1
  if [ -n "$SCALING_FILTER" ] && [ "$scaling" != "$SCALING_FILTER" ]; then
    return 1
  fi
  [ ${#FILTER[@]} -eq 0 ] && return 0
  for f in "${FILTER[@]}"; do [ "$f" = "$run_id" ] && return 0; done
  return 1
}

# The five panels plot_methods_umaps.py writes, in the order it writes them.
PANELS=(cohort_integrated cohort_int_vs_unint celltype_integrated
        celltype_int_vs_unint cohort_celltype_integrated)

# Already plotted? All five panels must be there: the plotter saves them one by
# one, so a crash halfway leaves a partial directory that must be redone. Each
# run costs a UMAP on 620k cells, hence the default is to leave a finished one
# alone; --force redraws.
already_done() {
  local dir="$1" m="$2" p
  [ "$FORCE" -eq 0 ] || return 1
  for p in "${PANELS[@]}"; do [ -f "$dir/${m}_${p}.png" ] || return 1; done
  return 0
}

# Best QC view for a run's output type(s): prefer the corrected embedding.
plot_type() {
  case ",$1," in
    *,embed,*) echo embed ;;
    *,full,*)  echo full ;;
    *,knn,*)   echo knn ;;
    *)         echo "$1" ;;
  esac
}

# Both files the plotter reads. The reference is checked too, not just the
# integrated object: it travels separately (sync_to_cluster.sh sends it as its own
# file) and a scaled run whose scaled reference has not arrived fails on the
# plotter's assertion, deep into the run and far less obviously.
missing_inputs() {
  local int_path="$1" ref_path="$2" missing=""
  [ -f "$int_path" ] || missing="$int_path"
  [ -f "$ref_path" ] || missing="${missing:+$missing, }$ref_path"
  echo "$missing"
}

# SLURM mode: submit one array task per matching run.

# Submit one group of array indices. A group exists only because `--mem` is a
# submission-level option and the two knn rows need three times what the others
# do; everything else about the two submissions is identical.
submit_group() {
  local label="$1" mem="$2"; shift 2
  [ $# -gt 0 ] || return 0

  local spec; spec="$(IFS=,; echo "$*")%${MAX_CONCURRENT}"
  local export_vars="ALL,DATA_DIR=$DATA_DIR,PLOT_DIR=$SCRIPT_DIR"
  [ "$RECOMPUTE_UMAP" -eq 1 ] && export_vars="$export_vars,PLOT_RECOMPUTE_UMAP=1"
  # --chdir, and the logs/ directory created HERE: the #SBATCH --output path is
  # relative to the job's working directory and SLURM opens it before the script
  # runs, so the `mkdir -p logs` inside submit_plots.slurm is too late to help on
  # a first submission. --chdir also makes the log location independent of where
  # this driver was invoked from.
  local opts=(--chdir="$SCRIPT_DIR" --export="$export_vars" --array="$spec" --mem="$mem")
  [ -n "$EXCLUDE_NODES" ] && opts+=(--exclude="$EXCLUDE_NODES")

  if [ "$DRY_RUN" -eq 1 ]; then
    echo "[dry] $label: sbatch ${opts[*]} $SCRIPT_DIR/submit_plots.slurm"
    return 0
  fi
  mkdir -p "$SCRIPT_DIR/logs"
  local j; j="$(sbatch --parsable "${opts[@]}" "$SCRIPT_DIR/submit_plots.slurm")"
  echo "submitted $label array: job $j   (--mem=$mem, tasks $spec)"
}

run_slurm() {
  # The array index is the grid's data row number, which is how submit_plots.slurm
  # finds its row. Rows already plotted or missing an input are dropped from the
  # spec rather than submitted and skipped inside the job: a queued task that
  # exits immediately still costs a slot against the throttle.
  local indices=() indices_knn=() total=0 n_have=0 n_skip=0 idx=0
  local run_id method language env scaling input output types reference hvgs
  {
    read -r _header
    while IFS=$'\t' read -r run_id method language env scaling input output types reference hvgs; do
      [ -n "$run_id" ] || continue
      idx=$((idx + 1)); total=$((total + 1))
      want "$run_id" "$scaling" "$language" || continue

      if already_done "$FIG_DIR/$run_id" "$method"; then
        echo "[have] $run_id: already plotted, not submitting"
        n_have=$((n_have + 1)); continue
      fi
      local missing; missing="$(missing_inputs "$DATA_DIR/$output" "$DATA_DIR/$reference")"
      if [ -n "$missing" ]; then
        echo "[skip] $run_id: input not here yet ($missing)"
        n_skip=$((n_skip + 1)); continue
      fi
      if [ "$(plot_type "$types")" = "knn" ]; then indices_knn+=("$idx"); else indices+=("$idx"); fi
    done
  } < "$GRID"

  local n_sel=$(( ${#indices[@]} + ${#indices_knn[@]} ))
  if [ "$n_sel" -eq 0 ]; then
    if [ "$n_have" -gt 0 ]; then
      echo "nothing to submit: all $n_have matching run(s) are already plotted (--force to redraw)"
      return 0
    fi
    echo "no grid row matches the given filters" >&2; exit 1
  fi

  echo "selected ${n_sel}/${total} run(s), ${n_have} already plotted, ${n_skip} not integrated yet, ${MAX_CONCURRENT} at a time"
  echo "  figures -> $FIG_DIR/<run_id>${EXCLUDE_NODES:+   excluding: $EXCLUDE_NODES}"
  submit_group "std" "$MEM_STD" ${indices[@]+"${indices[@]}"}
  submit_group "knn" "$MEM_KNN" ${indices_knn[@]+"${indices_knn[@]}"}
  echo "the panels land in the repo, not in DATA_DIR; bring them back with:"
  echo "  rsync -av <cluster>:$FIG_DIR/ ./figures/"
}

# Local mode: plot every matching run here, in sequence.

run_local() {
  if [ "$DRY_RUN" -eq 0 ]; then
    set +u; eval "$(conda shell.bash hook)"; set -u
    conda_guarded activate "$ENV"
  fi

  local n_done=0 n_skip=0 n_have=0 n_match=0
  local failed=()
  local run_id method language env scaling input output types reference hvgs
  {
    read -r _header
    while IFS=$'\t' read -r run_id method language env scaling input output types reference hvgs; do
      [ -n "$run_id" ] || continue
      if ! want "$run_id" "$scaling" "$language"; then continue; fi
      n_match=$((n_match + 1))

      local int_path="$DATA_DIR/$output" ref_path="$DATA_DIR/$reference"
      local out_dir="$FIG_DIR/$run_id" ptype; ptype="$(plot_type "$types")"

      # Done first, missing-input second: what a run needed to produce its panels is
      # not part of whether those panels exist. Reporting a finished run as
      # not-integrated because its input is unreadable would be a false alarm.
      if already_done "$out_dir" "$method"; then
        echo "[have] $run_id: already plotted, skipping ($out_dir)"
        n_have=$((n_have + 1)); continue
      fi

      local missing; missing="$(missing_inputs "$int_path" "$ref_path")"
      if [ -n "$missing" ]; then
        echo "[skip] $run_id: input missing ($missing)"
        n_skip=$((n_skip + 1)); continue
      fi

      local umap_args=()
      [ "$RECOMPUTE_UMAP" -eq 1 ] && umap_args+=(--recompute-umap)

      if [ "$DRY_RUN" -eq 1 ]; then
        echo "[dry] $run_id  (--type $ptype)  -> $out_dir"
        echo "      python $PLOTTER -m $method --type $ptype -i $int_path -u $ref_path -o $out_dir ${umap_args[*]}"
        n_done=$((n_done + 1)); continue
      fi

      echo "==================================================================="
      echo ">>> $run_id  ($method, --type $ptype)"
      echo "==================================================================="
      # One failure must not take the rest of the batch with it: the runs are
      # independent, each costs ~25 min of UMAP, and a whole-batch abort is how
      # fifteen of them were lost on 2026-07-31. Failures are collected and
      # reported at the end, and the script exits non-zero so a caller still
      # notices. </dev/null so the plotter cannot consume the grid this loop is
      # reading from stdin.
      if python "$PLOTTER" -m "$method" --type "$ptype" \
           -i "$int_path" -u "$ref_path" -o "$out_dir" "${umap_args[@]}" </dev/null; then
        n_done=$((n_done + 1)); echo "[ok] $run_id -> $out_dir"
      else
        failed+=("$run_id"); echo "[fail] $run_id: see the traceback above, continuing" >&2
      fi
    done
  } < "$GRID"

  [ "$DRY_RUN" -eq 0 ] && conda_guarded deactivate

  if [ "$n_match" -eq 0 ]; then
    echo "no grid row matches the given filters" >&2; exit 1
  fi

  echo "==================================================================="
  echo "done: $n_done run(s) plotted, $n_have already plotted, $n_skip skipped, ${#failed[@]} failed"
  if [ "$n_have" -gt 0 ] && [ "$FORCE" -eq 0 ]; then
    echo "      (--force to redraw the $n_have existing figure set(s))"
  fi
  if [ ${#failed[@]} -gt 0 ]; then
    echo "      failed: ${failed[*]}"
    return 1
  fi
  return 0
}


if [ "$USE_SLURM" -eq 1 ]; then run_slurm; else run_local; fi
