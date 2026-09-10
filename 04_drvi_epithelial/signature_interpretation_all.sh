#!/usr/bin/env bash
#
# 04_3 - 04_7 signature interpretation: the whole chain in one command.
#
# Phase-level driver, as in 01_pre_processing/preprocessing_all.sh: the steps live in their
# own folders and this walks them in order. Runs the four scripts that turn the 04_2
# embedding and the collaborator's text files into the convergence table:
#
#   1. 04_3_signatures/build_signatures_epi.py  the .txt files   -> signatures_<coll>.gmt
#                                                                +  coverage / jaccard tables
#   2. 04_5_cell_first/cell_first_epi.py        shiao_epi.h5ad   -> signature_scores_<coll>_<run>.csv
#                                                                +  confounder / target tables
#   3. 04_3_signatures/signature_composition_epi.py  the two above -> signature_concentration /
#                                                                signature_gene_contribution tables
#   4. 04_6_factor_first/factor_first_epi.py    embed_<run>.h5ad -> factor_first_top200_<coll>_<run>.tsv
#   5. 04_7_convergence/convergence_epi.py      the tables above -> convergence_<coll>_<run>.csv
#   6. 04_8_cycle_confound/cycle_confound_epi.py  the tables above -> cycle_confound_* tables
#                                                                +  the four-panel cycle figure
#
# Step 3 lives in the 04_3 folder because it characterises the signature collection, but it runs
# after 04_5 so it can correlate each gene against the score 04_5 computed. It feeds nothing:
# drop it from a run and every table below is unchanged.
#
# THE COLLECTION. The same four steps are run over two independent bodies of prior knowledge,
# declared in utils/sig_collections.py:
#
#   scie     stemness x immunogenicity, the ten lab lists plus CytoTRACE2   (the default)
#   emt      the EMT axis, nine lists on epithelial / hybrid / mesenchymal
#   gavish   the pan-cancer metaprograms of Gavish et al. 2023, as a VOCABULARY - the 22
#            relevant to a TNBC by default, all 41 with --all-metaprograms
#
# They share no output: every table and figure goes to <tables|figures>/<collection>/ and
# carries the collection in its filename, and 04_6 corrects its FDR inside one collection, so
# running one cannot move a single number of the other. Pick one with --collection:
#
#   ./signature_interpretation_all.sh                       # scie
#   ./signature_interpretation_all.sh --collection emt      # the same chain, EMT lists
#   ./signature_interpretation_all.sh --collection gavish   # the metaprograms, see below
#
# --all-metaprograms widens that last one from the 22 TNBC-relevant metaprograms to all 41.
# It is a different collection on disk - slug `gavish` instead of `gavish_tnbc` - so the two
# widths never overwrite each other, and the driver's [have] check follows the slug rather
# than the flag. It does nothing to scie or emt.
#
# GAVISH IS A SHORTER CHAIN, and the driver does not have to be told: it defines no target
# region (`Collection.has_target` is False), so 04_5 skips the quadrant and everything read
# off it, 04_7 skips the target-axis test, and the two steps that exist ONLY to interrogate a
# target region stop with a message instead of running:
#
#   cycle   04_8 - its cell-cycle question is answered for this collection by the
#           `confounders` table 04_5 writes for every readout;
#   (04_9 is not in this driver at all and stops the same way.)
#
# So the useful run is the five steps before it, which is what the default chain does anyway:
#
#   ./signature_interpretation_all.sh --collection gavish signatures cellfirst composition factorfirst convergence
#
# 04_4_cytotrace2 is NOT in this chain, on purpose, and the reason is the environment, not
# the runtime: cytotrace2-py pins numpy<2.0.0 and cannot be installed into benchmark-py-r
# without rolling that stack back (see environments/cytotrace2-py.yml). It is the one step
# that runs under a different interpreter, so it is asked for by name and pointed at that
# env through PYTHON, which every step in this driver reads:
#
#   PYTHON=~/miniconda3/envs/cytotrace2-py/bin/python ./signature_interpretation_all.sh cytotrace
#
# Runtime is ~55 min for the 29 cohorts on 12 cores, ~4 GB of temporary .txt matrices that
# the script deletes at the end, and a few GB of peak memory - no cluster needed.
#
# and re-run 04_5 afterwards, which picks the .csv up automatically and adds CytoTRACE2 as a
# sixth, non-circular stemness readout. Without it Route A runs on the lab's five stemness
# lists alone and says so in its output. It belongs to the scie collection only: the emt run
# ignores the .csv even when it is there.
#
# Resuming is the default: a step whose output already exists is reported [have] and skipped.
# --force re-runs everything. The check is existence only, so delete a file truncated by a
# crash before resuming.
#
# Usage:
#   export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
#   ./signature_interpretation_all.sh                    # the four local steps, resuming
#   ./signature_interpretation_all.sh --collection emt   # the same, on the EMT lists
#   ./signature_interpretation_all.sh --force            # re-run everything
#   ./signature_interpretation_all.sh --dry-run          # print what would run
#   ./signature_interpretation_all.sh --collection gavish --all-metaprograms   # all 41
#   ./signature_interpretation_all.sh cellfirst convergence   # only the named steps
#   PYTHON=.../envs/cytotrace2-py/bin/python ./signature_interpretation_all.sh cytotrace
#
# Step names: signatures, cytotrace, cellfirst, composition, factorfirst, convergence, cycle.
# Logs go to 04_drvi_epithelial/logs/ as well as to the terminal.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="$SCRIPT_DIR/logs"
PYTHON="${PYTHON:-python3}"
RUN_ID="${RUN_ID:-drvi_epi_64}"
COLLECTION="${COLLECTION:-scie}"
ALL_METAPROGRAMS="${ALL_METAPROGRAMS:-0}"

# One row per step: name, script, the output that marks it done, and whether the default
# chain includes it.
STEP_NAMES=(signatures cytotrace cellfirst composition factorfirst convergence cycle)
# Paths are relative to this file, i.e. to the phase folder.
STEP_SCRIPTS=(
  04_3_signatures/build_signatures_epi.py
  04_4_cytotrace2/cytotrace2_epi.py
  04_5_cell_first/cell_first_epi.py
  04_3_signatures/signature_composition_epi.py
  04_6_factor_first/factor_first_epi.py
  04_7_convergence/convergence_epi.py
  04_8_cycle_confound/cycle_confound_epi.py
)
# __COLL__ is substituted with the collection below, so resuming is tracked per collection:
# a finished scie run does not make the emt run look done. Inside tables/ the run is a folder
# of its own as well (tables/<coll>/<run_id>/), which is where the ${RUN_ID} in these paths
# appears twice - once as the folder, once in the filename, exactly as `table_dir` builds it.
STEP_OUTPUTS=(
  "$SCRIPT_DIR/tables/__COLL__/signatures___COLL__.gmt"
  "__EPI__/cytotrace2_${RUN_ID}.csv"
  "__EPI__/signature_scores___COLL___${RUN_ID}.csv"
  "$SCRIPT_DIR/tables/__COLL__/${RUN_ID}/signature_concentration___COLL___${RUN_ID}.csv"
  "__EPI__/factor_first_top200___COLL___${RUN_ID}.tsv"
  "$SCRIPT_DIR/tables/__COLL__/${RUN_ID}/convergence___COLL___${RUN_ID}.csv"
  "$SCRIPT_DIR/tables/__COLL__/${RUN_ID}/cycle_confound_by_dimension___COLL___${RUN_ID}.csv"
)
STEP_IN_DEFAULT=(1 0 1 1 1 1 1)
# The flag each script spells --force with, if any. 01 and 05 always recompute from what is
# already on disk and take no such flag, so --force there just means "run it again".
STEP_FORCE_FLAG=("" "--force" "--overwrite" "" "--overwrite" "" "")

FORCE=0; DRY_RUN=0
STEP_FILTER=()
while [ $# -gt 0 ]; do
  case "$1" in
    --collection) COLLECTION="${2:?--collection needs a value}"; shift 2 ;;
    --collection=*) COLLECTION="${1#*=}"; shift ;;
    --all-metaprograms) ALL_METAPROGRAMS=1; shift ;;
    --force|-f) FORCE=1; shift ;;
    --dry-run|-n) DRY_RUN=1; shift ;;
    -h|--help) sed -n '2,89p' "${BASH_SOURCE[0]}"; exit 0 ;;
    -*) echo "unknown option: $1" >&2; exit 1 ;;
    *) STEP_FILTER+=("$1"); shift ;;
  esac
done

# Checked after the options, so --help works without an environment.
: "${DATA_DIR:?set DATA_DIR to the directory holding the datasets (outside the repo)}"
EPI_DIR="$DATA_DIR/04_epi"

# Fail on a typo here rather than three steps in, where it would surface as a missing table.
case "$COLLECTION" in
  scie|emt|gavish) ;;
  *) echo "unknown collection: $COLLECTION (have: scie, emt, gavish)" >&2; exit 1 ;;
esac

# THE SLUG IS NOT ALWAYS THE COLLECTION. `--collection gavish` scores the TNBC-relevant
# subset unless --all-metaprograms is given, and the subset is its own collection on disk:
# every output path below is built from the slug, not from what was typed. This mirrors
# `sig_collections.resolve()`, which is where the decision actually lives - the driver only
# needs it to look for the right file when it decides what is already done.
SLUG="$COLLECTION"
if [ "$COLLECTION" = "gavish" ] && [ "$ALL_METAPROGRAMS" -eq 0 ]; then
  SLUG="gavish_tnbc"
fi

in_list() { local x="$1"; shift; local e; for e in "$@"; do [ "$e" = "$x" ] && return 0; done; return 1; }

# Validate the step names before doing any work, so a typo fails immediately.
for want in "${STEP_FILTER[@]:-}"; do
  [ -z "$want" ] && continue
  in_list "$want" "${STEP_NAMES[@]}" || { echo "unknown step: $want (have: ${STEP_NAMES[*]})" >&2; exit 1; }
done

mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/signature_interpretation_${SLUG}_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "04_3 - 04_7 signature interpretation"
# The slug is echoed only when it differs from what was asked for, i.e. on the gavish subset.
echo "collection $COLLECTION$([ "$SLUG" != "$COLLECTION" ] && echo " -> $SLUG")"
echo "run id     $RUN_ID"
echo "DATA_DIR   $DATA_DIR"
echo "log        $LOG_FILE"
echo

for i in "${!STEP_NAMES[@]}"; do
  name="${STEP_NAMES[$i]}"
  script="${STEP_SCRIPTS[$i]}"
  output="${STEP_OUTPUTS[$i]//__EPI__/$EPI_DIR}"
  output="${output//__COLL__/$SLUG}"

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

  # 04_8 interrogates the target region, and `gavish` defines none. The script itself stops
  # with the explanation; the driver skips it first so a full chain still ends on "finished".
  if [ "$name" = "cycle" ] && [ "$COLLECTION" = "gavish" ]; then
    echo "[skip] $name  (the $COLLECTION collection defines no target region for it to read)"
    continue
  fi

  extra=()
  # 04_4 is the one step that is not collection-scoped: it computes a measurement, and it is
  # the scie collection that decides to use it.
  if [ "$name" != "cytotrace" ]; then
    extra+=(--collection "$COLLECTION")
    # Spelled as an `if`, not as `cond && ...`: a false test as the last command of the
    # loop body would take `set -e` down with it.
    if [ "$ALL_METAPROGRAMS" -eq 1 ]; then extra+=(--all-metaprograms); fi
  fi
  if [ "$FORCE" -eq 1 ] && [ -n "${STEP_FORCE_FLAG[$i]}" ]; then
    extra+=("${STEP_FORCE_FLAG[$i]}")
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
