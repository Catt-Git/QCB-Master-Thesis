#!/usr/bin/env bash
#
# 02 utils - upload the locally produced objects to the cluster.
#
# The usual split in this phase is integration LOCALLY (the local environment is
# verified and has a GPU) and metrics ON THE CLUSTER, which means the integrated
# objects have to travel. This is that step, as one command.
#
# It walks the same benchmark_grid.tsv as every other wrapper, so the file names
# are never typed by hand: `output`, `input` and `reference` are read from the
# grid and sent to the SAME relative path under the cluster's DATA_DIR. Add a run
# to the grid and it is picked up here too.
#
#   --what metrics (default) : the integrated objects + the references they are
#                              scored against - everything 02_4_metrics needs.
#   --what inputs            : the 02_1 prepared inputs (.h5ad and .rds) + the HVG
#                              lists - everything 02_2_integration needs, if you
#                              want to integrate on the cluster instead.
#   --what all               : both.
#
# Usage:
#   export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
#   export CLUSTER_HOST=albertoc@hpc.example.it        # user@host, required
#   utils/sync_to_cluster.sh --dry-run                 # list what would be sent
#   utils/sync_to_cluster.sh                           # the whole grid
#   utils/sync_to_cluster.sh --scaling unscaled        # the unscaled half
#   utils/sync_to_cluster.sh --method harmony          # one method
#   utils/sync_to_cluster.sh bbknn_unscaled drvi_unscaled_128   # named run_id(s)
#   utils/sync_to_cluster.sh --what inputs             # the 02_1 inputs instead
#   utils/sync_to_cluster.sh --scp                     # scp instead of rsync
#
# The remote DATA_DIR defaults to $CLUSTER_DATA_DIR, or to the thesis path below;
# it mirrors the local layout exactly, so nothing downstream needs reconfiguring:
#
#   $DATA_DIR/02_integration/harmony_unscaled.h5ad
#     -> $CLUSTER_HOST:$CLUSTER_DATA_DIR/02_integration/harmony_unscaled.h5ad
#
# Resuming is the default, as in run_all.sh and run_all_metrics.sh: a file already
# on the cluster with the same byte size is reported [have] and not sent again, so
# re-running after a dropped connection picks up where it stopped. rsync also
# resumes a half-transferred file in place (--partial); scp cannot, and restarts
# that one file. --force re-sends everything.
#
# Both modes open a single SSH connection and reuse it (ControlMaster), so the
# password / 2FA prompt happens once for the whole transfer, not once per file.
#
# --scaling, --method and run_id filters combine (a run must match all active).

set -euo pipefail

: "${DATA_DIR:?set DATA_DIR to the local directory holding the objects to upload}"
: "${CLUSTER_HOST:?set CLUSTER_HOST to user@hostname of the cluster}"

UTILS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GRID="$UTILS_DIR/../benchmark_grid.tsv"

# Where this writes on the cluster is the same directory 02_4_metrics reads from,
# so it is defined once, in cluster_env.sh, and not here.
. "$UTILS_DIR/cluster_env.sh"
REMOTE="$CLUSTER_DATA_DIR"

# Read by the R methods only, so they are not in any grid column; they belong to
# the `inputs` set all the same.
HVG_LISTS=(shiao_hvg_2k_unintegrated_list.csv shiao_hvg_2k_unintegrated_list.rds)

[ -f "$GRID" ] || { echo "grid not found: $GRID" >&2; exit 1; }

in_list() { local x="$1"; shift; local e; for e in "$@"; do [ "$e" = "$x" ] && return 0; done; return 1; }
human() { numfmt --to=iec --suffix=B "$1" 2>/dev/null || echo "$1 bytes"; }

# Parse the command line.
WHAT="metrics"; USE_SCP=0; DRY_RUN=0; FORCE=0
SCALING_FILTER=""
METHOD_FILTER=(); RUNID_FILTER=()
while [ $# -gt 0 ]; do
  case "$1" in
    --what) WHAT="${2:-}"; shift 2 ;;
    --what=*) WHAT="${1#*=}"; shift ;;
    --scp) USE_SCP=1; shift ;;
    --force|-f) FORCE=1; shift ;;
    --dry-run|-n) DRY_RUN=1; shift ;;
    --scaling) SCALING_FILTER="${2:-}"; shift 2 ;;
    --scaling=*) SCALING_FILTER="${1#*=}"; shift ;;
    --method) METHOD_FILTER+=("${2:-}"); shift 2 ;;
    --method=*) METHOD_FILTER+=("${1#*=}"); shift ;;
    -h|--help) sed -n '2,50p' "${BASH_SOURCE[0]}"; exit 0 ;;
    -*) echo "unknown option: $1" >&2; exit 1 ;;
    *) RUNID_FILTER+=("$1"); shift ;;
  esac
done
case "$WHAT" in
  metrics|inputs|all) ;;
  *) echo "--what must be 'metrics', 'inputs' or 'all', got '$WHAT'" >&2; exit 1 ;;
esac
if [ -n "$SCALING_FILTER" ] && [ "$SCALING_FILTER" != "scaled" ] && [ "$SCALING_FILTER" != "unscaled" ]; then
  echo "--scaling must be 'scaled' or 'unscaled', got '$SCALING_FILTER'" >&2; exit 1
fi

# Reject typos up front, so a misspelled run id cannot silently select nothing.
mapfile -t ALL_RUNIDS < <(awk -F'\t' 'NR>1 && $1!="" { print $1 }' "$GRID")
mapfile -t ALL_METHODS < <(awk -F'\t' 'NR>1 && $1!="" { print $2 }' "$GRID" | sort -u)
for r in ${RUNID_FILTER+"${RUNID_FILTER[@]}"}; do
  in_list "$r" "${ALL_RUNIDS[@]}" || { echo "unknown run id: $r" >&2; exit 1; }
done
for m in ${METHOD_FILTER+"${METHOD_FILTER[@]}"}; do
  in_list "$m" "${ALL_METHODS[@]}" || { echo "unknown method: $m" >&2; exit 1; }
done

want_row() {
  local run_id="$1" scaling="$2" method="$3"
  [ -n "$SCALING_FILTER" ] && [ "$scaling" != "$SCALING_FILTER" ] && return 1
  [ ${#METHOD_FILTER[@]} -gt 0 ] && ! in_list "$method" "${METHOD_FILTER[@]}" && return 1
  [ ${#RUNID_FILTER[@]} -gt 0 ] && ! in_list "$run_id" "${RUNID_FILTER[@]}" && return 1
  return 0
}

# Collect the paths to send, as paths RELATIVE to DATA_DIR: that is what the grid
# stores and what keeps the two DATA_DIRs identical in shape.

WANTED=()
add_path() { in_list "$1" ${WANTED+"${WANTED[@]}"} || WANTED+=("$1"); }

{
  read -r _header
  while IFS=$'\t' read -r run_id method language env scaling input output types reference hvgs; do
    [ -n "$run_id" ] || continue
    want_row "$run_id" "$scaling" "$method" || continue
    case "$WHAT" in
      metrics) add_path "$output"; add_path "$reference" ;;
      inputs)  add_path "$input";  add_path "$reference" ;;
      all)     add_path "$output"; add_path "$input"; add_path "$reference" ;;
    esac
  done
} < "$GRID"

if [ "$WHAT" != "metrics" ]; then
  for f in "${HVG_LISTS[@]}"; do add_path "$f"; done
fi

[ ${#WANTED[@]} -gt 0 ] || { echo "no run matches the given filters" >&2; exit 1; }

# Drop what is not here yet (an integration still running, the scaled .rds never
# built) and, unless --force, what the cluster already holds at the same size.

SSH_CTL="${TMPDIR:-/tmp}/sync_to_cluster_$$.sock"
SSH_OPTS=(-o ControlPath="$SSH_CTL")
cleanup() { ssh -O exit "${SSH_OPTS[@]}" "$CLUSTER_HOST" 2>/dev/null || true; rm -f "$SSH_CTL"; }
trap cleanup EXIT
echo ">>> opening the SSH connection to $CLUSTER_HOST (authenticate once)"
ssh -fN -o ControlMaster=yes -o ControlPersist=8h "${SSH_OPTS[@]}" "$CLUSTER_HOST"

declare -A REMOTE_SIZE=()
if [ "$FORCE" -eq 0 ]; then
  # One round trip for the whole list: stat prints "<size> <relative path>" for
  # what exists and nothing for what does not.
  while read -r size rel; do
    [ -n "${rel:-}" ] && REMOTE_SIZE["$rel"]="$size"
  done < <(printf '%s\n' "${WANTED[@]}" \
    | ssh "${SSH_OPTS[@]}" "$CLUSTER_HOST" \
        "cd '$REMOTE' 2>/dev/null && xargs -d '\n' -r stat -c '%s %n' 2>/dev/null" || true)
fi

SEND=(); total=0; n_have=0; n_missing=0
for rel in "${WANTED[@]}"; do
  local_path="$DATA_DIR/$rel"
  if [ ! -f "$local_path" ]; then
    echo "[skip] $rel: not produced locally yet"; n_missing=$((n_missing + 1)); continue
  fi
  size="$(stat -c %s "$local_path")"
  if [ "${REMOTE_SIZE[$rel]:-}" = "$size" ]; then
    echo "[have] $rel: already on the cluster ($(human "$size"))"; n_have=$((n_have + 1)); continue
  fi
  SEND+=("$rel"); total=$((total + size))
done

if [ ${#SEND[@]} -eq 0 ]; then
  echo "nothing to send: $n_have file(s) already there, $n_missing missing locally"
  exit 0
fi

echo "==================================================================="
echo "sending ${#SEND[@]} file(s), $(human "$total") -> $CLUSTER_HOST:$REMOTE"
printf '  %s\n' "${SEND[@]}"
echo "==================================================================="

if [ "$DRY_RUN" -eq 1 ]; then
  echo "[dry] nothing transferred"
  exit 0
fi

if [ "$USE_SCP" -eq 1 ]; then
  # scp cannot create the remote tree, so do it first, then send file by file.
  mapfile -t subdirs < <(for rel in "${SEND[@]}"; do d="$(dirname "$rel")"; [ "$d" != "." ] && echo "$d"; done | sort -u)
  ssh "${SSH_OPTS[@]}" "$CLUSTER_HOST" \
    "mkdir -p '$REMOTE' $(printf "'%s' " ${subdirs+"${subdirs[@]/#/$REMOTE/}"})"
  for rel in "${SEND[@]}"; do
    echo "--- $rel"
    scp "${SSH_OPTS[@]}" -p "$DATA_DIR/$rel" "$CLUSTER_HOST:$REMOTE/$rel"
  done
else
  # --files-from implies --relative, so the 02_integration/ prefix is recreated
  # on the far side and half-transferred files resume in place (--partial).
  ssh "${SSH_OPTS[@]}" "$CLUSTER_HOST" "mkdir -p '$REMOTE'"
  printf '%s\n' "${SEND[@]}" | rsync -av --partial --progress \
    -e "ssh ${SSH_OPTS[*]}" --files-from=- "$DATA_DIR/" "$CLUSTER_HOST:$REMOTE/"
fi

echo "==================================================================="
echo "done: ${#SEND[@]} sent, $n_have already there, $n_missing missing locally"
echo "      -> $CLUSTER_HOST:$REMOTE"
[ "$WHAT" = "metrics" ] && cat <<EOF
next, on the cluster:
  export DATA_DIR=$REMOTE
  cd <repo>/02_integration_benchmark/02_4_metrics && mkdir -p logs
  ./run_all_metrics.sh --slurm${SCALING_FILTER:+ --scaling $SCALING_FILTER} --dry-run
EOF
exit 0
