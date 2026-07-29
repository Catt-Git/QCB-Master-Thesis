#!/bin/bash
#SBATCH --job-name=sample_map
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# submit_sample_map.slurm
# SLURM wrapper for step 2: build_sample_map.sh, which turns ena_full_metadata.tsv into the SRR -> patient/timepoint/fraction map.
# Pure text processing on a 266-row table: 1 core and 2G are generous. As for step 1, the wrapper exists to keep the whole phase submittable as a chain of jobs.
# No DATA_DIR: input and outputs are small TSVs that stay IN the repo, next to build_sample_map.sh.
#
# Usage (submit FROM the 00_2_sample_mapping folder, after step 1):
#   sbatch submit_sample_map.slurm

set -euo pipefail
# Exit on error, undefined variables, and pipe failures, so a parsing failure is reported as FAILED rather than producing a silently truncated map.

hostname; date
echo "Job name: $SLURM_JOB_NAME"; echo "Job ID: $SLURM_JOB_ID"
# Log the node, timestamp, and job identifiers into the .out file for traceability.

awk --version 2>/dev/null | head -1 | grep -qi 'gnu awk' \
  || { echo "[error] awk on this node is not GNU awk; build_sample_map.sh needs the 3-argument match(str, re, arr), a gawk extension. Load a gawk module or activate an environment providing it." >&2; exit 1; }
# Guard: the library_name parser relies on gawk's 3-argument match(). Cluster nodes sometimes default to mawk/busybox awk, where that call is a hard error halfway through the job; checking here fails in one second with an actionable message.

cd "$SLURM_SUBMIT_DIR"
bash build_sample_map.sh
# Build the map. Submit with 'sbatch submit_sample_map.slurm' FROM the 00_2_sample_mapping folder: SLURM copies this file to a spool directory, so "$0" cannot be used to locate build_sample_map.sh; SLURM_SUBMIT_DIR (the directory sbatch was called from) is the reliable anchor.

date
