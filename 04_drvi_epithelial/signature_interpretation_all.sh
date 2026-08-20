#!/usr/bin/env bash
#
# 04_3 - 04_7 signature interpretation: the whole chain in one command.
#
# Phase-level driver, as in 01_pre_processing/preprocessing_all.sh: the steps live in their
# own folders and this walks them in order. Runs the four scripts that turn the 04_2
# embedding and the lab's text files into the convergence table:
#
#   1. 04_3_signatures/build_signatures_epi.py  11 .txt files    -> lab_signatures.gmt
#                                                                +  coverage / jaccard tables
#   2. 04_5_cell_first/cell_first_epi.py        shiao_epi.h5ad   -> signature_scores_<run>.csv
#                                                                +  confounder / quadrant tables
#   3. 04_6_factor_first/factor_first_epi.py    embed_<run>.h5ad -> factor_first_top200_<run>.tsv
#   4. 04_7_convergence/convergence_epi.py      the tables above -> convergence_<run>.csv
#
# 04_4_cytotrace2 is NOT in this chain, on purpose. Its dependency (cytotrace2-py) is not
# part of benchmark-py-r, its runtime is minutes per patient across 29 patients rather than
# seconds, and it is the one part of the chain meant to run on the cluster - see
# 04_4_cytotrace2/submit_cytotrace2_epi.slurm. Run it explicitly by name:
#
#   ./signature_interpretation_all.sh cytotrace
#
# and re-run 04_5 afterwards, which picks the .csv up automatically and adds CytoTRACE2 as a
# seventh, non-circular stemness readout. Without it Route A runs on the lab's six stemness
# lists alone and says so in its output.
#
# Resuming is the default: a step whose output already exists is reported [have] and skipped.
# --force re-runs everything. The check is existence only, so delete a file truncated by a
# crash before resuming.
#
# Usage:
#   export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
#   ./signature_interpretation_all.sh                    # the four local steps, resuming
#   ./signature_interpretation_all.sh --force            # re-run everything
#   ./signature_interpretation_all.sh --dry-run          # print what would run
#   ./signature_interpretation_all.sh cellfirst convergence   # only the named steps
#   ./signature_interpretation_all.sh cytotrace          # the cluster step, locally
#
# Step names: signatures, cytotrace, cellfirst, factorfirst, convergence.
# Logs go to 04_drvi_epithelial/logs/ as well as to the terminal.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="$SCRIPT_DIR/logs"
PYTHON="${PYTHON:-python3}"
RUN_ID="${RUN_ID:-drvi_epi_64}"

# One row per step: name, script, the output that marks it done, and whether the default
# chain includes it.
STEP_NAMES=(signatures cytotrace cellfirst factorfirst convergence)
# Paths are relative to this file, i.e. to the phase folder.
STEP_SCRIPTS=(
  04_3_signatures/build_signatures_epi.py
  04_4_cytotrace2/cytotrace2_epi.py
  04_5_cell_first/cell_first_epi.py
  04_6_factor_first/factor_first_epi.py
  04_7_convergence/convergence_epi.py
)
STEP_OUTPUTS=(
  "$SCRIPT_DIR/tables/lab_signatures.gmt"
  "__EPI__/cytotrace2_${RUN_ID}.csv"
  "__EPI__/signature_scores_${RUN_ID}.csv"
  "__EPI__/factor_first_top200_${RUN_ID}.tsv"
  "$SCRIPT_DIR/tables/convergence_${RUN_ID}.csv"
)
STEP_IN_DEFAULT=(1 0 1 1 1)
# The flag each script spells --force with, if any. 01 and 05 always recompute from what is
# already on disk and take no such flag, so --force there just means "run it again".
STEP_FORCE_FLAG=("" "--force" "--overwrite" "--overwrite" "")

FORCE=0; DRY_RUN=0
STEP_FILTER=()
while [ $# -gt 0 ]; do
  case "$1" in
    --force|-f) FORCE=1; shift ;;
    --dry-run|-n) DRY_RUN=1; shift ;;
    -h|--help) sed -n '2,45p' "${BASH_SOURCE[0]}"; exit 0 ;;
    -*) echo "unknown option: $1" >&2; exit 1 ;;
    *) STEP_FILTER+=("$1"); shift ;;
  esac
done

# Checked after the options, so --help works without an environment.
: "${DATA_DIR:?set DATA_DIR to the directory holding the datasets (outside the repo)}"
EPI_DIR="$DATA_DIR/04_epi"

in_list() { local x="$1"; shift; local e; for e in "$@"; do [ "$e" = "$x" ] && return 0; done; return 1; }

# Validate the step names before doing any work, so a typo fails immediately.
for want in "${STEP_FILTER[@]:-}"; do
  [ -z "$want" ] && continue
  in_list "$want" "${STEP_NAMES[@]}" || { echo "unknown step: $want (have: ${STEP_NAMES[*]})" >&2; exit 1; }
done

mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/signature_interpretation_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "04_3 - 04_7 signature interpretation"
echo "run id     $RUN_ID"
echo "DATA_DIR   $DATA_DIR"
echo "log        $LOG_FILE"
echo

for i in "${!STEP_NAMES[@]}"; do
  name="${STEP_NAMES[$i]}"
  script="${STEP_SCRIPTS[$i]}"
  output="${STEP_OUTPUTS[$i]//__EPI__/$EPI_DIR}"

  if [ "${#STEP_FILTER[@]}" -gt 0 ]; then
    in_list "$name" "${STEP_FILTER[@]}" || continue
  elif [ "${STEP_IN_DEFAULT[$i]}" -eq 0 ]; then
    echo "[skip] $name  (not in the default chain; run it by name, see --help)"
    continue
  fi

  if [ "$FORCE" -eq 0 ] && [ -e "$output" ]; then
    echo "[have] $name  -> $output"
    continue
  fi

  extra=()
  if [ "$FORCE" -eq 1 ] && [ -n "${STEP_FORCE_FLAG[$i]}" ]; then
    extra=("${STEP_FORCE_FLAG[$i]}")
  fi

  if [ "$DRY_RUN" -eq 1 ]; then
    echo "[dry ] $name  -> $PYTHON $script ${extra[*]:-}"
    continue
  fi

  echo "[run ] $name  -> $script"
  # Spelled out rather than "${extra[@]:-}": an empty array expanded that way hands the
  # script one empty argument, which argparse rejects. Same guard as submit_drvi_epi.slurm.
  step_dir="$SCRIPT_DIR/$(dirname "$script")"
  step_file="$(basename "$script")"
  if [ "${#extra[@]}" -gt 0 ]; then
    ( cd "$step_dir" && "$PYTHON" "$step_file" "${extra[@]}" )
  else
    ( cd "$step_dir" && "$PYTHON" "$step_file" )
  fi
  echo "[done] $name"
  echo
done

echo "all requested steps finished."
