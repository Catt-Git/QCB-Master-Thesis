#!/usr/bin/env bash
#
# 05_4 - 05_8 signature interpretation: the whole chain in one command.
#
# Phase-level driver, the counterpart of 04's `signature_interpretation_all.sh` and built the
# same way: the steps live in their own folders and this walks them in order. It turns the
# 05_3 embedding and the collaborator's text files into the convergence table:
#
#   1. 05_4_signatures/build_signatures_tum.py   the .txt files  -> signatures_<coll>_<cmp>.gmt
#                                                                +  coverage / jaccard tables
#   2. 05_6_cell_first/cell_first_tum.py         shiao_tum.h5ad  -> signature_scores_<coll>_<run>.csv
#                                                                +  confounder / target / row-order tables
#   3. 05_4_signatures/signature_composition_tum.py  the two above -> signature_concentration /
#                                                                signature_gene_contribution tables
#   4. 05_7_factor_first/factor_first_tum.py     embed_<run>.h5ad -> factor_first_top200_<coll>_<run>.tsv
#                                                                +  dim_geneset_signed_significance
#   5. 05_8_convergence/convergence_tum.py       the tables above -> convergence_<coll>_<run>.csv
#   6. 05_4_signatures/emt_vs_gavish_tum.py      the two .gmt     -> emt_vs_gavish_* tables (emt only)
#
# THE ORDER IS NOT COSMETIC. 05_7 reads the row order off the `dimension_row_order` table that
# 05_6 writes, and 05_8 reads five tables from 05_6 and 05_7; step 3 runs after 05_6 so it can
# correlate each gene against the score 05_6 computed. Steps 1, 2, 4, 5 are a chain; 3 and 6
# feed nothing and can be dropped from a run without moving a number below them.
#
# THE COLLECTION. The same chain is run over three independent bodies of prior knowledge,
# declared in utils/sig_collections.py:
#
#   scie     stemness x immunogenicity, the ten lab lists plus CytoTRACE2   (the default)
#   emt      the EMT axis, nine lists on epithelial / hybrid / mesenchymal
#   gavish   the pan-cancer metaprograms of Gavish et al. 2023, as a VOCABULARY - the 22
#            relevant to a TNBC by default (slug `gavish_tnbc`), all 41 with --all-metaprograms
#
# They share no output: every table and figure goes to <tables|figures>/<collection>/<run_id>/
# and carries the collection in its filename, and 05_7 corrects its FDR inside one collection,
# so running one cannot move a single number of another.
#
# WHICH RUN. `CELL_SET`, `N_LATENT` and `HVG_SET` select the 05_3 run this reads, exactly as
# in the step READMEs - the driver only passes the environment through:
#
#   N_LATENT=64 ./signature_interpretation_all.sh                  # drvi_tum_64_nomt
#   CELL_SET=epi N_LATENT=64 ./signature_interpretation_all.sh     # drvi_epicnv_64_nomt
#
# PRUNE_VANISHED, the control this driver exists to make cheap. Unset, the stage keeps every
# dimension and every direction, which is what the phase reports and what the existing
# `tables/` and `figures/` hold. Set, 05_6 drops the dimensions DRVI flagged vanished, 05_7
# drops the vanished DIRECTIONS (a finer cut - a dimension can vanish on one side only), 05_8
# inherits the shorter dimension list from 05_6's table, and the BH denominator of 05_7
# follows the directions actually tested:
#
#   PRUNE_VANISHED=1 N_LATENT=64 ./signature_interpretation_all.sh --collection scie
#
# NOTHING OF THE REPORTED RUN IS OVERWRITTEN when it is set: `OUT_TAG` in signature_common
# moves TABLE_DIR to `tables_pruned/`, the figure root to `figures_pruned/`, and the two 05_7
# .tsv files to a `_pruned` name. The per-cell files that pruning cannot change - CytoTRACE2
# and `signature_scores_<coll>_<run>.csv` - are READ from their existing paths and left alone,
# so the control costs no re-scoring and 05_5 never runs twice.
#
# 05_5_cytotrace2 is NOT in the default chain, on purpose, and the reason is the environment,
# not the runtime: cytotrace2-py pins numpy<2.0.0 and cannot live in the same env as the rest
# (see environments/cytotrace2-py.yml). It is asked for by name and pointed at that env
# through PYTHON, which every step in this driver reads:
#
#   PYTHON=~/miniconda3/envs/cytotrace2-py/bin/python ./signature_interpretation_all.sh cytotrace
#
# and 05_6 picks the .csv up automatically afterwards, adding CytoTRACE2 as a sixth,
# non-circular stemness readout. It belongs to the scie collection only.
#
# Resuming is the default: a step whose output already exists is reported [have] and skipped.
# --force re-runs everything. The check is existence only, so delete a file truncated by a
# crash before resuming.
#
# Usage:
#   export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
#   N_LATENT=64 ./signature_interpretation_all.sh                     # scie, resuming
#   N_LATENT=64 ./signature_interpretation_all.sh --collection emt    # the EMT lists
#   N_LATENT=64 ./signature_interpretation_all.sh --collection gavish # the 22 TNBC metaprograms
#   N_LATENT=64 ./signature_interpretation_all.sh --force             # re-run everything
#   N_LATENT=64 ./signature_interpretation_all.sh --dry-run           # print what would run
#   N_LATENT=64 ./signature_interpretation_all.sh cellfirst convergence  # only the named steps
#
# Step names: signatures, cytotrace, cellfirst, composition, factorfirst, convergence, emtgavish.
# Logs go to each step's own logs/ folder, as the hand-run logs already there do, and to the
# terminal.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="$SCRIPT_DIR/logs"
PYTHON="${PYTHON:-python3}"
COLLECTION="${COLLECTION:-scie}"
ALL_METAPROGRAMS="${ALL_METAPROGRAMS:-0}"

FORCE=0; DRY_RUN=0
STEP_FILTER=()
while [ $# -gt 0 ]; do
  case "$1" in
    --collection) COLLECTION="${2:?--collection needs a value}"; shift 2 ;;
    --collection=*) COLLECTION="${1#*=}"; shift ;;
    --all-metaprograms) ALL_METAPROGRAMS=1; shift ;;
    --force|-f) FORCE=1; shift ;;
    --dry-run|-n) DRY_RUN=1; shift ;;
    -h|--help) sed -n '2,82p' "${BASH_SOURCE[0]}"; exit 0 ;;
    -*) echo "unknown option: $1" >&2; exit 1 ;;
    *) STEP_FILTER+=("$1"); shift ;;
  esac
done

# Checked after the options, so --help works without an environment.
: "${DATA_DIR:?set DATA_DIR to the directory holding the datasets (outside the repo)}"
TUM_DATA_DIR="$DATA_DIR/05_tum"

case "$COLLECTION" in
  scie|emt|gavish) ;;
  *) echo "unknown collection: $COLLECTION (have: scie, emt, gavish)" >&2; exit 1 ;;
esac

# THE RUN ID AND THE OUTPUT TAG ARE DERIVED, NEVER RETYPED. Both live in Python - `cell_set`
# owns the run id, `signature_common` owns the tag - and a second copy of either rule in bash
# is exactly how a driver ends up checking for files the steps do not write. So ask.
RESOLVE_PY="$SCRIPT_DIR/utils/resolve_run.py"
read -r RUN_ID COMPARTMENT OUT_TAG < <(cd "$SCRIPT_DIR" && "$PYTHON" "$RESOLVE_PY")
# The tag is printed as '-' when it is empty, so the line always has three fields and
# `read` cannot leave a stale value in the third variable.
if [ "$OUT_TAG" = "-" ]; then OUT_TAG=""; fi
TABLE_ROOT="$SCRIPT_DIR/tables${OUT_TAG}"

# THE SLUG IS NOT ALWAYS THE COLLECTION - `--collection gavish` writes under `gavish_tnbc`
# unless --all-metaprograms widens it. This mirrors `sig_collections.resolve()`, which is
# where the decision actually lives; the driver needs it only to look for the right file.
SLUG="$COLLECTION"
if [ "$COLLECTION" = "gavish" ] && [ "$ALL_METAPROGRAMS" -eq 0 ]; then
  SLUG="gavish_tnbc"
fi

# One row per step: name, script, its log basename, the output that marks it done, and
# whether the default chain includes it. __COLL__ is the slug, __TUM__ the heavy-data dir.
STEP_NAMES=(signatures cytotrace cellfirst composition factorfirst convergence emtgavish)
STEP_SCRIPTS=(
  05_4_signatures/build_signatures_tum.py
  05_5_cytotrace2/cytotrace2_tum.py
  05_6_cell_first/cell_first_tum.py
  05_4_signatures/signature_composition_tum.py
  05_7_factor_first/factor_first_tum.py
  05_8_convergence/convergence_tum.py
  05_4_signatures/emt_vs_gavish_tum.py
)
STEP_LOGS=(signatures cytotrace2 cell_first composition factor_first convergence emt_vs_gavish)
# The .gmt sits ABOVE the run folder (it depends on the object, not the run); everything else
# is tables/<coll>/<run_id>/, which is where ${RUN_ID} appears twice - once as the folder,
# once in the filename, exactly as `table_dir` builds it.
STEP_OUTPUTS=(
  "$TABLE_ROOT/__COLL__/signatures___COLL___${COMPARTMENT}.gmt"
  "__TUM__/cytotrace2_${RUN_ID}.csv"
  "$TABLE_ROOT/__COLL__/${RUN_ID}/dimension_row_order___COLL___${RUN_ID}.csv"
  "$TABLE_ROOT/__COLL__/${RUN_ID}/signature_concentration___COLL___${RUN_ID}.csv"
  "$TABLE_ROOT/__COLL__/${RUN_ID}/dim_geneset_signed_significance___COLL___${RUN_ID}.csv"
  "$TABLE_ROOT/__COLL__/${RUN_ID}/convergence___COLL___${RUN_ID}.csv"
  "$TABLE_ROOT/emt/${RUN_ID}/emt_vs_gavish_pairs_emt_${RUN_ID}.csv"
)
STEP_IN_DEFAULT=(1 0 1 1 1 1 0)
# The flag each script spells --force with, if any.
STEP_FORCE_FLAG=("" "--force" "--overwrite" "" "--overwrite" "" "")

in_list() { local x="$1"; shift; local e; for e in "$@"; do [ "$e" = "$x" ] && return 0; done; return 1; }

for want in "${STEP_FILTER[@]:-}"; do
  [ -z "$want" ] && continue
  in_list "$want" "${STEP_NAMES[@]}" || { echo "unknown step: $want (have: ${STEP_NAMES[*]})" >&2; exit 1; }
done

mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/signature_interpretation_${SLUG}_${RUN_ID}${OUT_TAG}_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "05_4 - 05_8 signature interpretation"
echo "collection $COLLECTION$([ "$SLUG" != "$COLLECTION" ] && echo " -> $SLUG")"
echo "run id     $RUN_ID  (compartment $COMPARTMENT)"
echo "pruning    $([ -n "$OUT_TAG" ] && echo "PRUNE_VANISHED on -> outputs under tables${OUT_TAG}/ and figures${OUT_TAG}/" || echo "off, the run this phase reports")"
echo "DATA_DIR   $DATA_DIR"
echo "log        $LOG_FILE"
echo

for i in "${!STEP_NAMES[@]}"; do
  name="${STEP_NAMES[$i]}"
  script="${STEP_SCRIPTS[$i]}"
  output="${STEP_OUTPUTS[$i]//__TUM__/$TUM_DATA_DIR}"
  output="${output//__COLL__/$SLUG}"

  if [ "${#STEP_FILTER[@]}" -gt 0 ]; then
    in_list "$name" "${STEP_FILTER[@]}" || continue
  elif [ "${STEP_IN_DEFAULT[$i]}" -eq 0 ]; then
    echo "[skip] $name  (not in the default chain; run it by name, see --help)"
    continue
  fi

  # The cross-collection check reads the EMT lists against Gavish's and is meaningless for
  # any other collection: it is asked for by name, on an emt run.
  if [ "$name" = "emtgavish" ] && [ "$COLLECTION" != "emt" ]; then
    echo "[skip] $name  (it compares the emt collection with gavish; run it on --collection emt)"
    continue
  fi

  if [ "$FORCE" -eq 0 ] && [ -e "$output" ]; then
    echo "[have] $name  -> $output"
    continue
  fi

  extra=()
  # 05_5 computes a measurement on raw counts and is not collection-scoped; 05_4's cross
  # check is pinned to the emt collection in the script itself and takes no --collection.
  if [ "$name" != "cytotrace" ] && [ "$name" != "emtgavish" ]; then
    extra+=(--collection "$COLLECTION")
  fi
  if [ "$ALL_METAPROGRAMS" -eq 1 ] && [ "$name" != "cytotrace" ]; then
    # Spelled as an `if`, not as `cond && ...`: a false test as the last command of the
    # loop body would take `set -e` down with it.
    extra+=(--all-metaprograms)
  fi
  if [ "$FORCE" -eq 1 ] && [ -n "${STEP_FORCE_FLAG[$i]}" ]; then
    extra+=("${STEP_FORCE_FLAG[$i]}")
  fi

  if [ "$DRY_RUN" -eq 1 ]; then
    echo "[dry ] $name  -> $PYTHON $script ${extra[*]:-}"
    continue
  fi

  echo "[run ] $name  -> $script"
  step_dir="$SCRIPT_DIR/$(dirname "$script")"
  step_file="$(basename "$script")"
  step_log="$step_dir/logs/${STEP_LOGS[$i]}_${SLUG}_${RUN_ID}${OUT_TAG}.log"
  mkdir -p "$step_dir/logs"
  # Spelled out rather than "${extra[@]:-}": an empty array expanded that way hands the
  # script one empty argument, which argparse rejects.
  if [ "${#extra[@]}" -gt 0 ]; then
    ( cd "$step_dir" && "$PYTHON" "$step_file" "${extra[@]}" ) | tee "$step_log"
  else
    ( cd "$step_dir" && "$PYTHON" "$step_file" ) | tee "$step_log"
  fi
  echo "[done] $name  (step log: $step_log)"
  echo
done

echo "all requested steps finished."
