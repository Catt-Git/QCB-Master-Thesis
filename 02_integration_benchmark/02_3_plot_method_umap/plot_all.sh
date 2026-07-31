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
# Usage:
#   export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
#   ./plot_all.sh                       # every integrated run
#   ./plot_all.sh harmony_unscaled ...  # only the named run_id(s)
#   ./plot_all.sh --scaling unscaled    # only the unscaled half
#   ./plot_all.sh --dry-run             # print what would run, do nothing
#   ./plot_all.sh --force               # redraw runs whose panels already exist
#   ./plot_all.sh --force --recompute-umap   # redraw *and* lay the UMAPs out again
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

: "${DATA_DIR:?set DATA_DIR to the directory holding the prepared inputs}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GRID="$SCRIPT_DIR/../benchmark_grid.tsv"
FIG_DIR="$SCRIPT_DIR/../figures"
PLOTTER="$SCRIPT_DIR/plot_methods_umaps.py"

[ -f "$GRID" ] || { echo "grid not found: $GRID" >&2; exit 1; }

# conda's env activation scripts (e.g. MKL's) read unbound variables, which
# aborts under `set -u`; run conda with nounset temporarily disabled.
conda_guarded() { set +u; conda "$@"; set -u; }

# Parse the command line: --scaling, --dry-run, --force, plus optional run_id filters.
SCALING_FILTER=""
DRY_RUN=0
FORCE=0
RECOMPUTE_UMAP=0
FILTER=()
while [ $# -gt 0 ]; do
  case "$1" in
    --scaling) SCALING_FILTER="${2:-}"; shift 2 ;;
    --scaling=*) SCALING_FILTER="${1#*=}"; shift ;;
    --force|-f) FORCE=1; shift ;;
    --dry-run|-n) DRY_RUN=1; shift ;;
    --recompute-umap) RECOMPUTE_UMAP=1; shift ;;
    *) FILTER+=("$1"); shift ;;
  esac
done
if [ -n "$SCALING_FILTER" ] && [ "$SCALING_FILTER" != "scaled" ] && [ "$SCALING_FILTER" != "unscaled" ]; then
  echo "--scaling must be 'scaled' or 'unscaled', got '$SCALING_FILTER'" >&2
  exit 1
fi

# A row is wanted if it matches every active filter (run_id and scaling).
want() {
  local run_id="$1" scaling="$2" f
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

if [ "$DRY_RUN" -eq 0 ]; then
  set +u; eval "$(conda shell.bash hook)"; set -u
  conda_guarded activate benchmark-py-r
fi

n_done=0
n_skip=0
n_have=0

{
  read -r _header
  while IFS=$'\t' read -r run_id method language env scaling input output types reference hvgs; do
    [ -n "$run_id" ] || continue
    if ! want "$run_id" "$scaling"; then n_skip=$((n_skip + 1)); continue; fi

    int_path="$DATA_DIR/$output"
    ref_path="$DATA_DIR/$reference"
    out_dir="$FIG_DIR/$run_id"
    ptype="$(plot_type "$types")"

    # Done first, missing-input second: what a run needed to produce its panels is
    # not part of whether those panels exist. Reporting a finished run as
    # not-integrated because its input is unreadable would be a false alarm.
    if already_done "$out_dir" "$method"; then
      echo "[have] $run_id: already plotted, skipping ($out_dir)"
      n_have=$((n_have + 1))
      continue
    fi

    if [ ! -f "$int_path" ]; then
      echo "[skip] $run_id: not integrated yet ($int_path)"
      n_skip=$((n_skip + 1))
      continue
    fi

    UMAP_ARGS=()
    [ "$RECOMPUTE_UMAP" -eq 1 ] && UMAP_ARGS+=(--recompute-umap)

    if [ "$DRY_RUN" -eq 1 ]; then
      echo "[dry] $run_id  (--type $ptype)  -> $out_dir"
      echo "      python $PLOTTER -m $method --type $ptype -i $int_path -u $ref_path -o $out_dir ${UMAP_ARGS[*]}"
      n_done=$((n_done + 1))
      continue
    fi

    echo "==================================================================="
    echo ">>> $run_id  ($method, --type $ptype)"
    echo "==================================================================="
    python "$PLOTTER" -m "$method" --type "$ptype" \
      -i "$int_path" -u "$ref_path" -o "$out_dir" "${UMAP_ARGS[@]}"
    n_done=$((n_done + 1))
    echo "[ok] $run_id -> $out_dir"
  done
} < "$GRID"

[ "$DRY_RUN" -eq 0 ] && conda_guarded deactivate

echo "==================================================================="
echo "done: $n_done run(s) plotted, $n_have already plotted, $n_skip skipped"
if [ "$n_have" -gt 0 ] && [ "$FORCE" -eq 0 ]; then
  echo "      (--force to redraw the $n_have existing figure set(s))"
fi
