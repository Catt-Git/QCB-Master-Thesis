#!/usr/bin/env bash
#
# 02_2 integration: the integration step, runnable locally OR on SLURM.
#
# Reads benchmark_grid.tsv and runs each integration, picking the dispatcher and
# conda environment from the row. This is the single place that turns the grid
# into commands, so the benchmark matrix stays a piece of data rather than logic
# scattered across scripts.
#
#   default (local): runs each integration here, in sequence, using the local
#                    environments (benchmark-py-r, and scgen-py for scGen).
#   --slurm        : submits submit_integration.slurm restricted to the matching
#                    grid rows, throttled to 3 concurrent tasks (at most 3 nodes in
#                    `normal`, CPU). On the cluster every method uses catalano_env
#                    except scGen, which uses $SCGEN_ENV (see that file's header).
#
# Per row, by language:
#   python (benchmark-py-r)  run_integration.py         -> <run_id>.h5ad
#   python (scgen-py)        run_scgen.py               -> <run_id>.h5ad
#   R      (benchmark-py-r)  run_integration.R          -> <run_id>.rds
#                            then rds_to_h5ad.py        -> <run_id>.h5ad
#                            the .rds is then deleted (--keep-rds to keep it):
#                            nothing downstream reads it and it is the largest
#                            artefact of the phase (15 GB for fastMNN alone).
#   notebook                 skipped (DRVI, run by hand to check and adjust n_latent)
#
# Embeddings (embed / full,embed rows) are also written to
# 02_embeddings/<run_id>.npy: small, durable artifacts for the figures step.
#
# The three trained methods checkpoint their fitted model to
# 02_<method>/<run_id>_model.pt (02_scvi/, 02_scanvi/, 02_scgen/, alongside the
# 02_drvi/ the DRVI notebook writes) and reload it instead of training again on a
# re-run. Training is the long half of those runs and what follows it is the
# fragile half, so a crash after training no longer costs the training. They are
# kept out of 02_integration/ so that directory holds only the integrated objects
# the benchmark scores. Delete a .pt to force that run to train again.
#
# Usage:
#   export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
#   ./run_all.sh                        # local, every row
#   ./run_all.sh --slurm                # submit every row (3 at a time)
#   ./run_all.sh harmony_unscaled ...   # only the named run_id(s)
#   ./run_all.sh --method scanorama     # one method, both scalings
#   ./run_all.sh --scaling unscaled     # every unscaled run (half the grid)
#   ./run_all.sh --dry-run              # print what would run, do nothing
#   ./run_all.sh --force                # re-run and overwrite finished runs
#   ./run_all.sh --keep-rds             # keep the R intermediates
#
# Resuming is the default: a row whose `output` already exists is reported as
# [have] and skipped, so re-running the same command after a crash picks up where
# it stopped. --force overwrites instead. The check is existence only, so delete
# an .h5ad truncated by a crash mid-write before resuming.
# --scaling, --method and run_id filters combine (a row must match all active).
# Integrations are long; scVI/scANVI/scGen run on CPU on the cluster (slower).
#
# Usage examples:
# Run all scaled locally
# ./run_all.sh --scaling scaled
# Run all scaled on slurm server
# ./run_all.sh --scaling scaled --slurm

set -euo pipefail

: "${DATA_DIR:?set DATA_DIR to the directory holding the prepared inputs}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GRID="$SCRIPT_DIR/../benchmark_grid.tsv"
HVG_RDS="$DATA_DIR/shiao_hvg_2k_unintegrated_list.rds"
MAX_CONCURRENT=3                        # SLURM: at most 3 nodes in `normal`

# Reference *batches* for the anchor-based Seurat methods (CCA/RPCA). This is a
# method parameter, unrelated to the grid's `reference` column, which names the
# unintegrated .h5ad the metrics score against. Left empty, run_integration.R's
# own default (the three largest patients) applies; override with
#   SEURAT_REFERENCE=Patient53,Patient16 ./run_all.sh --method seurat_cca
SEURAT_REF_ARG=()
if [ -n "${SEURAT_REFERENCE:-}" ]; then SEURAT_REF_ARG=(--reference "$SEURAT_REFERENCE"); fi

[ -f "$GRID" ] || { echo "grid not found: $GRID" >&2; exit 1; }

# conda's env activation scripts (e.g. MKL's) read unbound variables, which
# aborts under `set -u`; run conda with nounset temporarily disabled.
conda_guarded() { set +u; conda "$@"; set -u; }
in_list() { local x="$1"; shift; local e; for e in "$@"; do [ "$e" = "$x" ] && return 0; done; return 1; }

# Called as `need_value "$@"` from an option branch: aborts if the option has no
# value after it. Without it the `shift 2` fails and the script dies silently.
need_value() {
  [ $# -ge 2 ] || { echo "$1 needs a value" >&2; exit 1; }
}

# Abort unless every given value appears in column `col` of the grid.
check_in_grid() {
  local what="$1" col="$2"; shift 2
  [ $# -gt 0 ] || return 0
  local valid; valid="$(awk -F'\t' -v c="$col" 'NR>1 && $1!="" { print $c }' "$GRID" | sort -u)"
  local v
  for v in "$@"; do
    if ! grep -qxF -- "$v" <<< "$valid"; then
      { echo "unknown $what '$v'; the grid has:"; sed 's/^/  /' <<< "$valid"; } >&2
      exit 1
    fi
  done
}

# Parse the command line.
USE_SLURM=0; DRY_RUN=0; KEEP_RDS=0; FORCE=0
SCALING_FILTER=""
METHOD_FILTER=(); RUNID_FILTER=()
while [ $# -gt 0 ]; do
  case "$1" in
    --slurm) USE_SLURM=1; shift ;;
    --keep-rds) KEEP_RDS=1; shift ;;
    --force|-f) FORCE=1; shift ;;
    --dry-run|-n) DRY_RUN=1; shift ;;
    --scaling) need_value "$@"; SCALING_FILTER="$2"; shift 2 ;;
    --scaling=*) SCALING_FILTER="${1#*=}"; shift ;;
    --method) need_value "$@"; METHOD_FILTER+=("$2"); shift 2 ;;
    --method=*) METHOD_FILTER+=("${1#*=}"); shift ;;
    # Anything else that looks like an option is a typo, not a run id. Without
    # this it would land in RUNID_FILTER, match no grid row, and the script would
    # report "0 run(s) executed" and exit 0 -- a typo indistinguishable from a
    # successful no-op.
    -*) echo "unknown option '$1'" >&2; exit 1 ;;
    *) RUNID_FILTER+=("$1"); shift ;;
  esac
done
if [ -n "$SCALING_FILTER" ] && [ "$SCALING_FILTER" != "scaled" ] && [ "$SCALING_FILTER" != "unscaled" ]; then
  echo "--scaling must be 'scaled' or 'unscaled', got '$SCALING_FILTER'" >&2; exit 1
fi

# Same reasoning for the values: a run id or method that is not in the grid is a
# typo, and silently matching nothing would look like success. Checked here, up
# front, so local and --slurm mode reject it identically.
check_in_grid "run id" 1 "${RUNID_FILTER[@]+"${RUNID_FILTER[@]}"}"
check_in_grid "method" 2 "${METHOD_FILTER[@]+"${METHOD_FILTER[@]}"}"

# A row is wanted if it matches every active filter (scaling, method, run_id).
want_row() {
  local run_id="$1" scaling="$2" method="$3"
  [ -n "$SCALING_FILTER" ] && [ "$scaling" != "$SCALING_FILTER" ] && return 1
  [ ${#METHOD_FILTER[@]} -gt 0 ] && ! in_list "$method" "${METHOD_FILTER[@]}" && return 1
  [ ${#RUNID_FILTER[@]} -gt 0 ] && ! in_list "$run_id" "${RUNID_FILTER[@]}" && return 1
  return 0
}

# Already integrated? The grid's `output` column is the single source of truth for
# where a run lands, so "does that file exist" is the whole check. Integrations
# cost hours, so the default is to leave a finished one alone and say so; --force
# re-runs and overwrites. Note this tests existence only: an .h5ad truncated by a
# crash mid-write looks finished here. Delete it (or check it with
# 02_4_metrics/check_integrations.py) before resuming.
already_done() { [ "$FORCE" -eq 0 ] && [ -f "$1" ]; }

# Whether the run has an embedding worth exporting (embed appears in types).
has_embed() { case ",$1," in *,embed,*) return 0 ;; *) return 1 ;; esac }

# Where a trained method checkpoints its model: its own flat per-method directory,
# rather than beside the output, so 02_integration/ holds only the integrated
# objects the metrics read. The run id goes in the file name, not in a
# subdirectory (see model_paths.py). Sets model_arg, empty for untrained methods.
set_model_arg() {
  model_arg=()
  case "$1" in
    scvi|scanvi|scgen) model_arg=(--model-dir "$DATA_DIR/02_${1}") ;;
  esac
}

# Remove an integrated .rds once it has been converted, unless --keep-rds.
drop_rds() {
  local p="$1"
  if [ "$KEEP_RDS" -eq 1 ]; then echo "[keep] $p"; return 0; fi
  [ -f "$p" ] || return 0
  echo "[clean] removing intermediate $(du -h "$p" | cut -f1)  $p"
  rm -f "$p"
}

# SLURM mode: submit submit_integration.slurm for the matching grid rows.

run_slurm() {
  # One array task = one grid row (integration is per run, not per type). Notebook
  # rows (DRVI) cannot be submitted, so they are left out of the array spec.
  # Finished runs are dropped from the array spec rather than submitted and
  # skipped inside the job: a queued task that exits immediately still costs a
  # slot against the %3 throttle.
  local indices=() total=0 n_have=0
  while IFS=$'\t' read -r idx run_id method scaling language output; do
    total=$((total + 1))
    [ "$language" = "notebook" ] && continue
    want_row "$run_id" "$scaling" "$method" || continue
    if already_done "$DATA_DIR/$output"; then
      echo "[have] $run_id: already integrated, not submitting ($DATA_DIR/$output)"
      n_have=$((n_have + 1)); continue
    fi
    indices+=("$idx")
  done < <(awk -F'\t' 'NR>1 && $1!="" { idx++; print idx"\t"$1"\t"$2"\t"$5"\t"$3"\t"$7 }' "$GRID")

  if [ ${#indices[@]} -eq 0 ]; then
    if [ "$n_have" -gt 0 ]; then
      echo "nothing to submit: all $n_have matching row(s) are already integrated (--force to re-run)"
      return 0
    fi
    echo "no grid row matches the given filters (or all are notebooks)" >&2; exit 1
  fi

  local spec; spec="$(IFS=,; echo "${indices[*]}")%${MAX_CONCURRENT}"
  echo "selected ${#indices[@]}/${total} row(s), ${n_have} already integrated, throttled to ${MAX_CONCURRENT}: $spec"

  local exports="ALL,DATA_DIR=$DATA_DIR,KEEP_RDS=$KEEP_RDS"
  if [ "$DRY_RUN" -eq 1 ]; then
    echo "[dry] sbatch --export=$exports --array=$spec $SCRIPT_DIR/submit_integration.slurm"
    return 0
  fi
  local j; j="$(sbatch --parsable --export="$exports" --array="$spec" "$SCRIPT_DIR/submit_integration.slurm")"
  echo "submitted integration array: job $j"
  echo "DRVI is a notebook: run it by hand (shiao_drvi_128.ipynb). Then score with 02_4_metrics/run_all_metrics.sh."
}

# Local mode: run every matching integration here, in sequence.

run_local() {
  # conda activate inside a non-interactive shell needs the hook sourced first.
  [ "$DRY_RUN" -eq 0 ] && { set +u; eval "$(conda shell.bash hook)"; set -u; }

  local n_done=0 n_skip=0 n_have=0
  {
    read -r _header
    while IFS=$'\t' read -r run_id method language env scaling input output types ref_h5ad hvgs; do
      [ -n "$run_id" ] || continue
      if ! want_row "$run_id" "$scaling" "$method"; then n_skip=$((n_skip + 1)); continue; fi

      local in_path="$DATA_DIR/$input" out_path="$DATA_DIR/$output"
      local emb_path="$DATA_DIR/02_embeddings/${run_id}.npy"

      if already_done "$out_path"; then
        echo "[have] $run_id: already integrated, skipping ($out_path)"
        n_have=$((n_have + 1)); continue
      fi

      if [ "$DRY_RUN" -eq 1 ]; then
        echo "[dry] $run_id  ($method, $language, $env)  -> $out_path"
        n_done=$((n_done + 1)); continue
      fi
      mkdir -p "$(dirname "$out_path")"

      echo "==================================================================="
      echo ">>> $run_id  ($method, $language, $env)"
      echo "==================================================================="

      case "$language" in
        python)
          conda_guarded activate "$env"
          local model_arg; set_model_arg "$method"
          if [ "$method" = "scgen" ]; then
            python "$SCRIPT_DIR/run_scgen.py" -i "$in_path" -o "$out_path" "${model_arg[@]}"
          else
            local emb_arg=(); has_embed "$types" && emb_arg=(--emb-out "$emb_path")
            python "$SCRIPT_DIR/run_integration.py" -m "$method" -i "$in_path" -o "$out_path" \
              "${emb_arg[@]}" "${model_arg[@]}"
          fi
          conda_guarded deactivate
          ;;
        R)
          conda_guarded activate "$env"
          local rds_out="$DATA_DIR/02_integration/${run_id}.rds"
          Rscript "$SCRIPT_DIR/run_integration.R" -m "$method" -i "$in_path" -o "$rds_out" \
            -v "$HVG_RDS" "${SEURAT_REF_ARG[@]}"
          local emb_arg=(); has_embed "$types" && emb_arg=(--emb-out "$emb_path")
          python "$SCRIPT_DIR/rds_to_h5ad.py" -i "$rds_out" -o "$out_path" --types "$types" "${emb_arg[@]}"
          # The .rds is a pure intermediate: nothing downstream reads it, and it is
          # the largest artefact of the phase (fastMNN alone writes 15 GB). Drop it
          # once rds_to_h5ad.py has returned 0, which means the .h5ad exists and
          # passed its type / finiteness checks. --keep-rds opts out.
          drop_rds "$rds_out"
          conda_guarded deactivate
          ;;
        notebook)
          echo "[skip] $run_id is a notebook (DRVI); run it by hand"; n_skip=$((n_skip + 1)); continue ;;
        *)
          echo "[error] unknown language '$language' for $run_id" >&2; exit 1 ;;
      esac
      n_done=$((n_done + 1)); echo "[ok] $run_id -> $out_path"
    done
  } < "$GRID"

  echo "==================================================================="
  echo "done: $n_done run(s) executed, $n_have already integrated, $n_skip filtered out"
  [ "$n_have" -gt 0 ] && [ "$FORCE" -eq 0 ] && \
    echo "      (--force to re-run and overwrite the $n_have existing output(s))"
  return 0
}


if [ "$USE_SLURM" -eq 1 ]; then run_slurm; else run_local; fi
