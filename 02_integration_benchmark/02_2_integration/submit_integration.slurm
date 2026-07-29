#!/usr/bin/env bash
#SBATCH --job-name=integration
#SBATCH --output=logs/integration_%A_%a.out
#SBATCH --error=logs/integration_%A_%a.err
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --array=1-17%3
#
# 02_2 integration: SLURM array (CPU).
#
# One array task = one grid row (integration is per run, not per output type), so
# --array over the 17 rows. The `%3` throttle keeps at most three tasks running at
# once (at most three nodes in `normal`); the rest queue and start as those finish.
# No `--time`: the partition's own default wall-clock limit applies.
#
# Environments on the cluster (override with the vars below):
#   INTEGRATION_ENV  every method except scGen         [default: catalano_env]
#   SCGEN_ENV        scGen only (its pinned stack)      [default: scgen-py]

# catalano_env (benchmark-hpc.yml) carries the full integration stack
# scGen needs its own env: create one on the cluster equivalent to the local 
# scgen-py and point SCGEN_ENV at it, or run scGen locally
# (./run_all.sh --method scgen) and submit the rest here.
#
# scVI/scANVI/scGen run on CPU here (no GPU on `normal`), so they are slower than
# locally; every other method is CPU-native. DRVI is a notebook and is not in the
# array (run_all.sh leaves it out; run it by hand).
#
# Usually launched through run_all.sh --slurm, which submits a filtered subset
# with the same %3 throttle. Direct use:
#   sbatch --export=ALL,DATA_DIR=/path/datasets submit_integration.slurm

set -euo pipefail
mkdir -p logs

: "${DATA_DIR:?set DATA_DIR (sbatch --export=ALL,DATA_DIR=... or export it first)}"
INTEGRATION_ENV="${INTEGRATION_ENV:-catalano_env}"
SCGEN_ENV="${SCGEN_ENV:-scgen-py}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GRID="$SCRIPT_DIR/../benchmark_grid.tsv"
HVG_RDS="$DATA_DIR/shiao_hvg_2k_unintegrated_list.rds"

has_embed() { case ",$1," in *,embed,*) return 0 ;; *) return 1 ;; esac }

# Where a trained method checkpoints its model: its own flat per-method directory,
# so 02_integration/ holds only the integrated objects the metrics read. Same
# layout as run_all.sh's local mode. Empty for the untrained methods.
set_model_arg() {
  model_arg=()
  case "$1" in
    scvi|scanvi|scgen) model_arg=(--model-dir "$DATA_DIR/02_${1}") ;;
  esac
}

# Pick this task's grid row (index = row number, matching run_all.sh's expansion).
JOB=$(awk -F'\t' 'NR>1 && $1!="" { idx++; if (idx == '"$SLURM_ARRAY_TASK_ID"') print }' "$GRID")
[ -n "$JOB" ] || { echo "no grid row for array index $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
IFS=$'\t' read -r run_id method language env scaling input output types ref_h5ad hvgs <<< "$JOB"

# Reference *batches* for the anchor-based Seurat methods (CCA/RPCA): a method
# parameter, not the grid's `reference` column (which names the unintegrated
# .h5ad the metrics score against). Empty means run_integration.R's own default.
SEURAT_REF_ARG=()
if [ -n "${SEURAT_REFERENCE:-}" ]; then SEURAT_REF_ARG=(--reference "$SEURAT_REFERENCE"); fi

if [ "$language" = "notebook" ]; then
  echo "[skip] $run_id is a notebook (DRVI); run it by hand"; exit 0
fi

# Local env names in the grid map to cluster envs: scGen keeps its own, the rest
# collapse onto catalano_env.
if [ "$method" = "scgen" ]; then ENVUSE="$SCGEN_ENV"; else ENVUSE="$INTEGRATION_ENV"; fi

in_path="$DATA_DIR/$input"
out_path="$DATA_DIR/$output"
emb_path="$DATA_DIR/02_embeddings/${run_id}.npy"
mkdir -p "$(dirname "$out_path")" "$(dirname "$emb_path")"

echo ">>> [$SLURM_ARRAY_TASK_ID] $run_id  ($method, $language, env=$ENVUSE)"

module load Miniconda3/202411
eval "$(/soft/system/software/Miniconda3/202411/bin/conda shell.bash hook)"
conda activate "$ENVUSE"

case "$language" in
  python)
    set_model_arg "$method"
    if [ "$method" = "scgen" ]; then
      python "$SCRIPT_DIR/run_scgen.py" -i "$in_path" -o "$out_path" "${model_arg[@]}"
    else
      emb_arg=(); has_embed "$types" && emb_arg=(--emb-out "$emb_path")
      python "$SCRIPT_DIR/run_integration.py" -m "$method" -i "$in_path" -o "$out_path" \
        "${emb_arg[@]}" "${model_arg[@]}"
    fi
    ;;
  R)
    rds_out="$DATA_DIR/02_integration/${run_id}.rds"
    Rscript "$SCRIPT_DIR/run_integration.R" -m "$method" -i "$in_path" -o "$rds_out" \
      -v "$HVG_RDS" "${SEURAT_REF_ARG[@]}"
    emb_arg=(); has_embed "$types" && emb_arg=(--emb-out "$emb_path")
    python "$SCRIPT_DIR/rds_to_h5ad.py" -i "$rds_out" -o "$out_path" --types "$types" "${emb_arg[@]}"
    # Pure intermediate, and the largest artefact of the phase: drop it once the
    # converter has returned 0 (the .h5ad exists and passed its checks). Set
    # KEEP_RDS=1 in the sbatch --export to keep it.
    if [ "${KEEP_RDS:-0}" = "1" ]; then echo "[keep] $rds_out"; else
      echo "[clean] removing intermediate $rds_out"; rm -f "$rds_out"; fi
    ;;
  *)
    echo "[error] unknown language '$language' for $run_id" >&2; exit 1 ;;
esac

echo "[ok] $run_id -> $out_path"
