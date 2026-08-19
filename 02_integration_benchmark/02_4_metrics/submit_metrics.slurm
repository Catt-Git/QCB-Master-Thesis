#!/usr/bin/env bash
#SBATCH --job-name=scib_metrics
#SBATCH --output=logs/metrics_%A_%a.out
#SBATCH --error=logs/metrics_%A_%a.err
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --array=1-22%5
#
# 02_4 metrics: SLURM array, the 12 metrics.
#
# One array task = one (run, output type). The task list is the grid expanded on
# its `types` column (18 rows -> 22 jobs), so --array=1-22. The `%5` throttle keeps
# at most five tasks running at once; the rest queue and start as those finish. A
# task takes 4 CPUs, not a whole node, so five fit across `normal` comfortably. No
# `--time` is set: the partition's own default wall-clock limit applies.
#
#   sbatch submit_metrics.slurm
#
# 96G covers every `full`/`embed` task: those score a 15-neighbour graph rebuilt by
# reduce_data. It does NOT cover the two `knn` (BBKNN) tasks. scib.integration.bbknn
# forces neighbors_within_batch=25 above 1e5 cells, which on 34 cohorts and 619,693
# cells gives 850 neighbours per cell in obsp['distances'] and 8.7e8 nonzeros in
# obsp['connectivities']. scanpy's leiden turns that adjacency into a Python
# list(zip(sources, targets)) - ~1e9 tuples, >100G - before leidenalg ever runs, so
# the first clustering resolution is OOM-killed (SIGKILL, ExitCode 9) at 64G; the
# retry was measured at 197 GiB. run_all_metrics.sh --slurm submits those two
# indices as their own group at MEM_KNN; by hand it is one option, and --mem also
# restricts the job to nodes that can hold it:
#
#   sbatch --mem=300G --array=1 ... submit_metrics.slurm     # bbknn_unscaled knn
#
# Usually launched through 02_4_metrics/run_all_metrics.sh --slurm, which submits
# a filtered subset with the same throttle, split by output type so each half gets
# the memory it needs.
# Requires DATA_DIR in the environment (or edit the default below) and the
# catalano_env conda environment on the node.

set -euo pipefail
mkdir -p logs

ENV="${METRICS_ENV:-catalano_env}"
TASK="shiao"
BATCH_KEY="cohort"; LABEL_KEY="cell_type"; ORGANISM="human"
# leiden implementation for the NMI/ARI clustering; run_all_metrics.sh exports it.
# Benchmark-wide setting: see metrics_shared.use_igraph_leiden.
CLUSTER_FLAVOR="${CLUSTER_FLAVOR:-igraph}"

# SLURM copies this file to the node's spool directory before running it, so
# "$0"/BASH_SOURCE point there and NOT at the repo - they cannot locate
# benchmark_grid.tsv or the sibling scripts. run_all_metrics.sh exports
# METRICS_DIR; SLURM_SUBMIT_DIR (where sbatch was called from) covers a direct
# `sbatch` issued from 02_4_metrics/, as phase 00 does.
SCRIPT_DIR="${METRICS_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
GRID="$SCRIPT_DIR/../benchmark_grid.tsv"
[ -f "$GRID" ] || { echo "grid not found: $GRID (set METRICS_DIR, or sbatch from 02_4_metrics/)" >&2; exit 1; }

# DATA_DIR normally arrives through --export from run_all_metrics.sh; a bare
# sbatch on the cluster falls back to the same data root sync_to_cluster.sh
# uploads into. One definition for both, in utils/cluster_env.sh.
. "$SCRIPT_DIR/../utils/cluster_env.sh"
resolve_data_dir || exit 1
ROOT="$DATA_DIR/02_metrics"

# Expand the grid into one line per (run, type): the SLURM_ARRAY_TASK_ID-th line
# is this task. Columns: run_id method scaling output types reference hvgs.
JOB=$(awk -F'\t' 'NR>1 && $1!="" {
        n=split($8, ts, ",");
        for (i=1;i<=n;i++) print $1"\t"$2"\t"$5"\t"$7"\t"ts[i]"\t"$9"\t"$10
      }' "$GRID" | sed -n "${SLURM_ARRAY_TASK_ID}p")
[ -n "$JOB" ] || { echo "no job for array index $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
IFS=$'\t' read -r run_id method scaling output type reference hvgs <<< "$JOB"

int_path="$DATA_DIR/$output"
ref_path="$DATA_DIR/$reference"
out_dir="$ROOT/$TASK/metrics/$scaling/hvg"
csv="$out_dir/${method}_${type}.csv"
mkdir -p "$out_dir"

echo ">>> [$SLURM_ARRAY_TASK_ID] $run_id  ($method, $type)"

# Unbuffered stdout. A task killed by the OOM killer loses whatever is still in
# Python's stdout buffer, which is precisely the progress trail that says whether it
# died in the leiden clustering or hours later in cLISI. Cheap insurance.
export PYTHONUNBUFFERED=1

# scib's LISI writes the whole connectivities graph to a text .mtx under TMPDIR
# before handing it to knn_graph.o. On a BBKNN graph (~9e8 nonzeros) that file is
# tens of GB, so TMPDIR must be disk-backed: on a tmpfs /tmp it counts against
# --mem and OOM-kills the job. Point METRICS_TMPDIR at node-local scratch if /tmp
# is small or memory-backed (check `df -h /tmp` on the node).
export TMPDIR="${METRICS_TMPDIR:-/tmp}/scib_${SLURM_JOB_ID:-$$}"
mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT

echo "[node] $(hostname)  mem=${SLURM_MEM_PER_NODE:-?}M  cpus=${SLURM_CPUS_PER_TASK:-?}  TMPDIR=$TMPDIR"

# `set -u` off across the whole conda block, the shell hook included: the hook
# auto-activates `base`, which sources this Miniconda's activate.d/openjdk_activate.sh,
# and that script expands ${JAVA_HOME} before defining it - fatal under -u, before
# `conda activate` is even reached. Same guard as conda_guarded() in run_all_metrics.sh.
set +u
module load Miniconda3/202411
eval "$(/soft/system/software/Miniconda3/202411/bin/conda shell.bash hook)"
conda activate "$ENV"
set -u

python "$SCRIPT_DIR/check_integrations.py" -i "$int_path" -u "$ref_path" \
  --types "$type" -b "$BATCH_KEY" -l "$LABEL_KEY"

python "$SCRIPT_DIR/metrics.py" -u "$ref_path" -i "$int_path" -o "$csv" \
  -m "$method" --type "$type" -b "$BATCH_KEY" -l "$LABEL_KEY" \
  --organism "$ORGANISM" --hvgs "$hvgs" --cluster-flavor "$CLUSTER_FLAVOR" -v

echo "[ok] $run_id $type -> $csv"
