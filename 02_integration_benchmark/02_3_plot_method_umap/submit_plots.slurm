#!/usr/bin/env bash
#SBATCH --job-name=scib_plots
#SBATCH --output=logs/plots_%A_%a.out
#SBATCH --error=logs/plots_%A_%a.err
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --array=1-17%5
#
# 02_3 method plots: SLURM array, one task per integration run.
#
# One array task = one grid row = one figures/<run_id>/ with its five panels. The
# grid has 17 data rows, hence --array=1-17; the index IS the row number, unlike
# 02_4 where a "full,embed" row expands into two tasks (here the two output types
# of such a row collapse to one plot, see plot_type() in plot_all.sh). The `%5`
# throttle keeps five running at once. No `--time`: the partition's own limit
# applies.
#
#   sbatch submit_plots.slurm
#
# 96G covers `full` (a 50-component PCA plus a 15-neighbour graph on 620k x 2k)
# and `embed` (the graph alone). It does NOT cover the two `knn` rows: BBKNN's
# corrected graph carries 8.7e8 nonzeros and sc.tl.umap lays that out directly.
# Those two get their own --mem, which also confines them to the big nodes:
#
#   sbatch --mem=300G --array=1,2 ... submit_plots.slurm      # the bbknn rows
#
# Usually launched through 02_3_plot_method_umap/plot_all.sh --slurm, which
# submits a filtered subset already split into those two memory groups, skips the
# runs whose panels exist and refuses to queue a task whose input has not arrived.
#
# Requires DATA_DIR in the environment (or the cluster data root below) and the
# conda environment on the node.

set -euo pipefail
mkdir -p logs

ENV="${PLOT_ENV:-catalano_env}"

# SLURM copies this file to the node's spool directory before running it, so
# "$0"/BASH_SOURCE point there and NOT at the repo - they cannot locate
# benchmark_grid.tsv or plot_methods_umaps.py. plot_all.sh exports PLOT_DIR;
# SLURM_SUBMIT_DIR (where sbatch was called from) covers a direct `sbatch` issued
# from 02_3_plot_method_umap/.
SCRIPT_DIR="${PLOT_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
PHASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
GRID="$PHASE_DIR/benchmark_grid.tsv"
FIG_DIR="$PHASE_DIR/figures"
PLOTTER="$SCRIPT_DIR/plot_methods_umaps.py"
[ -f "$GRID" ] || { echo "grid not found: $GRID (set PLOT_DIR, or sbatch from 02_3_plot_method_umap/)" >&2; exit 1; }

# DATA_DIR normally arrives through --export from plot_all.sh; a bare sbatch on
# the cluster falls back to the same data root sync_to_cluster.sh uploads into.
# One definition for both, in utils/cluster_env.sh.
. "$SCRIPT_DIR/../utils/cluster_env.sh"
resolve_data_dir || exit 1

# The SLURM_ARRAY_TASK_ID-th data row of the grid is this task.
# Columns: run_id method scaling output types reference.
JOB=$(awk -F'\t' 'NR>1 && $1!="" { print $1"\t"$2"\t"$5"\t"$7"\t"$8"\t"$9 }' "$GRID" \
      | sed -n "${SLURM_ARRAY_TASK_ID}p")
[ -n "$JOB" ] || { echo "no job for array index $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
IFS=$'\t' read -r run_id method scaling output types reference <<< "$JOB"

# Best QC view for this run's output type(s): prefer the corrected embedding.
# Kept identical to plot_type() in plot_all.sh - the two must not drift, or the
# driver's memory grouping would be computed for a different type than the job
# actually plots.
case ",$types," in
  *,embed,*) ptype=embed ;;
  *,full,*)  ptype=full ;;
  *,knn,*)   ptype=knn ;;
  *)         ptype="$types" ;;
esac

int_path="$DATA_DIR/$output"
ref_path="$DATA_DIR/$reference"
out_dir="$FIG_DIR/$run_id"

echo ">>> [$SLURM_ARRAY_TASK_ID] $run_id  ($method, --type $ptype)"
echo "[node] $(hostname)  mem=${SLURM_MEM_PER_NODE:-?}M  cpus=${SLURM_CPUS_PER_TASK:-?}"
echo "[out ] $out_dir"

# Unbuffered stdout: a task killed by the OOM killer loses whatever is still in
# Python's buffer, which is exactly the progress trail saying whether it died
# loading the object or 20 minutes later in sc.tl.umap.
export PYTHONUNBUFFERED=1
# Headless node: no display, and matplotlib must not try to find one.
export MPLBACKEND=Agg

# `set -u` off across the whole conda block, the shell hook included: the hook
# auto-activates `base`, which sources this Miniconda's activate.d/openjdk_activate.sh,
# and that script expands ${JAVA_HOME} before defining it - fatal under -u, before
# `conda activate` is even reached. Same guard as conda_guarded() in plot_all.sh.
set +u
module load Miniconda3/202411
eval "$(/soft/system/software/Miniconda3/202411/bin/conda shell.bash hook)"
conda activate "$ENV"
set -u

# plot_all.sh --slurm forwards --recompute-umap through this variable, since
# --export is the only channel into the job.
UMAP_ARGS=()
[ "${PLOT_RECOMPUTE_UMAP:-0}" = "1" ] && UMAP_ARGS+=(--recompute-umap)

python "$PLOTTER" -m "$method" --type "$ptype" \
  -i "$int_path" -u "$ref_path" -o "$out_dir" ${UMAP_ARGS[@]+"${UMAP_ARGS[@]}"}

echo "[ok] $run_id -> $out_dir"
