#!/usr/bin/env bash
#
# 01 pre-processing: the headless half of the phase, in one command.
#
# Runs the six scripts that turn the QC'd/Scrublet-filtered object written by
# qc_scrublet_filter.ipynb (01_2) into the definitive unintegrated object, in order:
#
#   1. 01_3_normalization/scran_norm.py            _scrublet.h5ad         -> _norm.h5ad
#   2. 01_4_cc_and_annotation/cell_cycle_score.py  _norm.h5ad             -> _norm_cc.h5ad
#   3. 01_4_cc_and_annotation/celltypist_annotation.py                    -> _norm_cc_annotated.h5ad
#   4. 01_4_cc_and_annotation/fraction_reassignment.py  _norm_cc_annotated.h5ad (IN PLACE)
#   5. 01_5_scib_pp/scib_reduce_data.py            _norm_cc_annotated.h5ad -> _reduced.h5ad
#                                                                          + shiao_hvg_2k_unintegrated_list.csv
#   6. 01_5_scib_pp/scib_clustering.py             _reduced.h5ad          -> shiao.h5ad
#
# (all filenames are prefixed all_samples_combined_scrublet and live flat in $DATA_DIR)
#
# Two modes, same step list and same resume logic:
#
#   default  every step here, in sequence, one process each. 619k cells, so scran is
#            slow but it does run locally.
#   --slurm  submit the whole chain as ONE job on `long` (submit_preprocessing_all.slurm,
#            470G / 8 cpus): the job re-invokes this script in default mode,
#            so the [have]/[run] decisions are taken inside the job, not at submit time.
#
#
# Resuming is the default: a step whose output already exists is reported as [have]
# and skipped, so re-running after a crash picks up where it stopped. --force re-runs
# everything. The check is existence only, so delete an .h5ad truncated by a crash
# mid-write before resuming. The exception is step 4, which rewrites its input in
# place and so has no output of its own: it is [have] when .obs['fraction'] of the
# annotated file already holds the recoded imm/non_imm categories instead of the CD45
# sort labels.
#
# Usage:
#   export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
#   ./preprocessing_all.sh                  # every step, locally, resuming
#   ./preprocessing_all.sh --slurm          # same chain as one job on `long`
#   ./preprocessing_all.sh --force          # re-run everything, overwriting
#   ./preprocessing_all.sh --dry-run        # print what would run, do nothing
#   ./preprocessing_all.sh reduce           # only the named step(s), in file order
#
# Step names: norm, cc, annot, fraction, reduce, cluster.
# Logs go to 01_pre_processing/logs/preprocessing_all_<timestamp>.log as well as to
# the terminal (in --slurm mode, to the job's logs/preprocessing_all_<jobid>.out too).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="$SCRIPT_DIR/logs"
PYTHON="${PYTHON:-python3}"

# One row per step: name, script, input it consumes, output that marks it done.
# Step 4 is in place, hence input == output; its done-check is special (see below).
STEP_NAMES=(norm cc annot fraction reduce cluster)
STEP_SCRIPTS=(
  01_3_normalization/scran_norm.py
  01_4_cc_and_annotation/cell_cycle_score.py
  01_4_cc_and_annotation/celltypist_annotation.py
  01_4_cc_and_annotation/fraction_reassignment.py
  01_5_scib_pp/scib_reduce_data.py
  01_5_scib_pp/scib_clustering.py
)
STEP_INPUTS=(
  all_samples_combined_scrublet.h5ad
  all_samples_combined_scrublet_norm.h5ad
  all_samples_combined_scrublet_norm_cc.h5ad
  all_samples_combined_scrublet_norm_cc_annotated.h5ad
  all_samples_combined_scrublet_norm_cc_annotated.h5ad
  all_samples_combined_scrublet_norm_cc_annotated_reduced.h5ad
)
STEP_OUTPUTS=(
  all_samples_combined_scrublet_norm.h5ad
  all_samples_combined_scrublet_norm_cc.h5ad
  all_samples_combined_scrublet_norm_cc_annotated.h5ad
  all_samples_combined_scrublet_norm_cc_annotated.h5ad
  all_samples_combined_scrublet_norm_cc_annotated_reduced.h5ad
  shiao.h5ad
)

FORCE=0; DRY_RUN=0; USE_SLURM=0
STEP_FILTER=()
while [ $# -gt 0 ]; do
  case "$1" in
    --force|-f) FORCE=1; shift ;;
    --dry-run|-n) DRY_RUN=1; shift ;;
    --slurm) USE_SLURM=1; shift ;;
    -h|--help) sed -n '2,50p' "${BASH_SOURCE[0]}"; exit 0 ;;
    -*) echo "unknown option: $1" >&2; exit 1 ;;
    *) STEP_FILTER+=("$1"); shift ;;
  esac
done

# Checked after the options, so --help works without an environment.
: "${DATA_DIR:?set DATA_DIR to the directory holding the datasets (outside the repo)}"

# The Tirosh/Regev list needed by the cc step. Same default as 01_4; override by
# exporting CC_GENES yourself.
export CC_GENES="${CC_GENES:-$DATA_DIR/regev_lab_cell_cycle_genes.txt}"

in_list() { local x="$1"; shift; local e; for e in "$@"; do [ "$e" = "$x" ] && return 0; done; return 1; }

# Validate the step names before doing any work, so a typo fails immediately rather
# than after the first hours of scran.
for requested in "${STEP_FILTER[@]:-}"; do
  [ -n "$requested" ] || continue
  in_list "$requested" "${STEP_NAMES[@]}" || {
    echo "unknown step '$requested'; valid steps: ${STEP_NAMES[*]}" >&2; exit 1; }
done

scheduled() {  # is this step name wanted on this invocation?
  local name="$1"
  [ ${#STEP_FILTER[@]} -eq 0 ] && return 0
  in_list "$name" "${STEP_FILTER[@]}"
}

# Step 4 has no output of its own: it recodes .obs['fraction'] inside the annotated
# file. Read the categorical's categories straight out of the HDF5 (cheap, no AnnData
# load): {imm, non_imm} means the recode already happened, {CD45+, CD45-} means it
# did not. Any failure (missing key, no h5py) counts as "not done" - the step is
# idempotent, and its own asserts give a better message than a probe could.
fraction_recoded() {
  local f="$1"
  [ -f "$f" ] || return 1
  "$PYTHON" - "$f" <<'PY' 2>/dev/null
import sys
try:
    import h5py
    with h5py.File(sys.argv[1], "r") as h:
        node = h["obs/fraction"]
        raw = node["categories"][:] if isinstance(node, h5py.Group) else []
    cats = {c.decode() if isinstance(c, bytes) else str(c) for c in raw}
except Exception:
    sys.exit(1)
sys.exit(0 if cats == {"imm", "non_imm"} else 1)
PY
}

step_done() {  # step index -> 0 if its result is already on disk
  local i="$1"
  [ "$FORCE" -eq 1 ] && return 1
  if [ "${STEP_NAMES[$i]}" = "fraction" ]; then
    fraction_recoded "$DATA_DIR/${STEP_OUTPUTS[$i]}"
  else
    [ -f "$DATA_DIR/${STEP_OUTPUTS[$i]}" ]
  fi
}

# Pre-flight: every scheduled step must be able to read its input, either because the
# file is already there or because an earlier scheduled step writes it. This catches
# the common case of asking for a late step whose intermediate has been deleted (the
# .h5ad chain of this phase is ~10GB a file, so those get cleaned up) before hours of
# compute go to waste. In --dry-run it only warns, so the plan stays inspectable.
produced=()
missing_input=0
for i in "${!STEP_NAMES[@]}"; do
  scheduled "${STEP_NAMES[$i]}" || continue
  step_done "$i" && { produced+=("${STEP_OUTPUTS[$i]}"); continue; }
  if [ ! -f "$DATA_DIR/${STEP_INPUTS[$i]}" ] && ! in_list "${STEP_INPUTS[$i]}" "${produced[@]:-}"; then
    echo "[error] ${STEP_NAMES[$i]}: missing input $DATA_DIR/${STEP_INPUTS[$i]}" >&2
    [ "$i" -eq 0 ] && echo "        run 01_2_qc_scrublet_filtering/qc_scrublet_filter.ipynb first: it is the step that" >&2
    [ "$i" -eq 0 ] && echo "        computes QC + Scrublet and applies the cell/gene filters." >&2
    [ "$i" -gt 0 ] && echo "        add the earlier step(s) to the command, or drop the step filter to run the whole chain." >&2
    missing_input=1
  fi
  produced+=("${STEP_OUTPUTS[$i]}")
done
if [ "$missing_input" -eq 1 ]; then
  [ "$DRY_RUN" -eq 0 ] && exit 1
  echo "[warn] --dry-run: continuing anyway, the plan below would fail as it stands"
fi

# --slurm: hand the whole chain to one job on `long` and return. The job runs this
# same script without --slurm, so there is exactly one place where the step list,
# the resume logic and the pre-flight live.
if [ "$USE_SLURM" -eq 1 ]; then
  mkdir -p "$LOG_DIR"
  # The filters travel as trailing sbatch arguments (sbatch forwards everything after
  # the script path to the script, which passes them on to this runner). Keeping them
  # out of --export matters: its value is a comma-separated list, and a variable
  # holding spaces is not something to rely on there.
  runner_args=()
  [ "$FORCE" -eq 1 ] && runner_args+=(--force)
  for s in "${STEP_FILTER[@]:-}"; do [ -n "$s" ] && runner_args+=("$s"); done
  # PREPROC_DIR: SLURM runs the batch script from a spool copy, so the job cannot
  # locate the repo on its own (same reason as 02's INTEGRATION_DIR).
  exports="ALL,DATA_DIR=$DATA_DIR,CC_GENES=$CC_GENES,PREPROC_DIR=$SCRIPT_DIR"
  # Absolute --output/--error override the batch script's relative ones, which would
  # otherwise land in whatever directory sbatch happened to be called from.
  sbatch_log=(--output "$LOG_DIR/preprocessing_all_%j.out" --error "$LOG_DIR/preprocessing_all_%j.err")
  if [ "$DRY_RUN" -eq 1 ]; then
    echo "[dry] sbatch ${sbatch_log[*]} --export=$exports" \
         "$SCRIPT_DIR/submit_preprocessing_all.slurm ${runner_args[*]:-}"
    exit 0
  fi
  # "${a[@]+"${a[@]}"}" and not "${a[@]:-}": the latter passes one empty argument when
  # the array is empty, which the runner would then read as a (never matching) step name.
  j="$(sbatch --parsable "${sbatch_log[@]}" --export="$exports" \
       "$SCRIPT_DIR/submit_preprocessing_all.slurm" "${runner_args[@]+"${runner_args[@]}"}")"
  echo "submitted the whole 01_3->01_5 chain as job $j (partition long)"
  echo "logs : $LOG_DIR/preprocessing_all_$j.out (and .err)"
  echo "watch: squeue -j $j ; tail -f $LOG_DIR/preprocessing_all_$j.out"
  exit 0
fi

mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/preprocessing_all_$(date +%Y%m%d_%H%M%S).log"
[ "$DRY_RUN" -eq 0 ] && exec > >(tee -a "$LOG_FILE") 2>&1

echo "DATA_DIR  : $DATA_DIR"
echo "CC_GENES  : $CC_GENES"
echo "python    : $($PYTHON --version 2>&1)"
[ -n "${SLURM_JOB_ID:-}" ] && echo "slurm job : $SLURM_JOB_ID on $(hostname)"
[ "$DRY_RUN" -eq 0 ] && echo "log       : $LOG_FILE"
echo

n_run=0; n_have=0; n_skip=0
overall_start=$SECONDS

for i in "${!STEP_NAMES[@]}"; do
  name="${STEP_NAMES[$i]}"
  script="$SCRIPT_DIR/${STEP_SCRIPTS[$i]}"
  output="$DATA_DIR/${STEP_OUTPUTS[$i]}"

  if ! scheduled "$name"; then
    n_skip=$((n_skip + 1)); continue
  fi

  if step_done "$i"; then
    if [ "$name" = "fraction" ]; then
      echo "[have] $name: fraction already recoded to imm/non_imm in $output, skipping"
    else
      echo "[have] $name: $output already exists, skipping"
    fi
    n_have=$((n_have + 1)); continue
  fi

  if [ "$DRY_RUN" -eq 1 ]; then
    echo "[dry] $name: $PYTHON $script -> $output"
    n_run=$((n_run + 1)); continue
  fi

  echo "==================================================================="
  echo ">>> $name  (${STEP_SCRIPTS[$i]})  $(date '+%Y-%m-%d %H:%M:%S')"
  echo "==================================================================="
  step_start=$SECONDS
  "$PYTHON" "$script"

  # A step that returns 0 without writing its output means the chain is broken
  # somewhere it did not raise; stop before the next step reads a stale file.
  [ -f "$output" ] || { echo "[error] $name finished but $output is missing" >&2; exit 1; }

  echo "[ok] $name in $(( (SECONDS - step_start) / 60 ))m$(( (SECONDS - step_start) % 60 ))s -> $output"
  n_run=$((n_run + 1))
done

echo "==================================================================="
echo "done in $(( (SECONDS - overall_start) / 60 ))m$(( (SECONDS - overall_start) % 60 ))s: \
$n_run step(s) run, $n_have already done, $n_skip filtered out"
[ "$n_have" -gt 0 ] && [ "$FORCE" -eq 0 ] && \
  echo "      (--force to re-run and overwrite the $n_have existing output(s))"
echo "next: 01_6_visualization/visualization_unintegrated.ipynb on $DATA_DIR/shiao.h5ad"
exit 0
