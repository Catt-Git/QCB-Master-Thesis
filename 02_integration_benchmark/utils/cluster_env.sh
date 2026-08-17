#!/usr/bin/env bash
#
# 02 utils - the one place that says where the data lives on the cluster.
#
# Two halves of this phase name the same directory: sync_to_cluster.sh writes into
# it from the laptop, and 02_4_metrics reads from it on the cluster through
# DATA_DIR. They used to hardcode it independently, and drifted: the upload landed
# in .../hopes_and_dreams/datasets while the metrics jobs read
# .../hopes_and_dreams, so a freshly uploaded half of the grid was reported "not
# integrated yet" while ten stale objects from an older manual copy scored fine.
# One definition, sourced by both, is the fix.
#
# Sourced, never executed:
#   . "$SCRIPT_DIR/../utils/cluster_env.sh"
#
# Both variables are read from the environment first, so a different machine or a
# second dataset root needs no edit here:
#   CLUSTER_HOST      user@hostname; only sync_to_cluster.sh needs it
#   CLUSTER_DATA_DIR  the data root ON the cluster. It is both where the upload
#                     goes and the DATA_DIR the metrics wrappers fall back to when
#                     they run there.

CLUSTER_DATA_DIR="${CLUSTER_DATA_DIR:-/users/genomics/albertoc/Tesi/hopes_and_dreams/datasets}"

# Resolve DATA_DIR for a script that *reads* the data (the metrics wrappers).
#
# An explicit DATA_DIR always wins: a local run points it at the laptop's own
# datasets/. With none set, use the cluster root if it is actually present, which
# is precisely the case of a job running on the cluster. Otherwise fail, as the
# wrappers did before.
#
# The warning covers the drift that motivated this file. On the cluster both paths
# can exist at once, and then a DATA_DIR aimed at the one sync_to_cluster.sh does
# not fill is a silently stale read - the worst kind, because the run succeeds.
resolve_data_dir() {
  if [ -n "${DATA_DIR:-}" ]; then
    if [ "${DATA_DIR%/}" != "${CLUSTER_DATA_DIR%/}" ] && [ -d "$CLUSTER_DATA_DIR" ]; then
      echo "[warn] DATA_DIR=$DATA_DIR, but the cluster data root is" >&2
      echo "[warn]   $CLUSTER_DATA_DIR" >&2
      echo "[warn] and it exists on this machine. sync_to_cluster.sh fills the latter," >&2
      echo "[warn] so one of the two is stale - check which before trusting the scores." >&2
    fi
    return 0
  fi
  if [ -d "$CLUSTER_DATA_DIR" ]; then
    export DATA_DIR="$CLUSTER_DATA_DIR"
    echo "[data] DATA_DIR unset, using the cluster data root: $DATA_DIR"
    return 0
  fi
  echo "set DATA_DIR: it is not in the environment and the cluster data root" >&2
  echo "($CLUSTER_DATA_DIR) does not exist here either." >&2
  return 1
}
