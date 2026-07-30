#!/bin/bash
#SBATCH --job-name=verify_integrity
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# submit_verify_integrity.slurm
# SLURM wrapper for step 3b: verify_integrity.sh, the content check that runs once the whole download array has finished.
# Sized for I/O, not CPU: 'gzip -t' decompresses all 298 files (~2.4 TB) serially, which is hours of streaming at a few hundred MB/s, hence the 'long' partition. 2 cores and 8G are plenty (gzip is single-threaded and streams).
# Run it AFTER the array is complete, not per task: it is a full-dataset audit, and the download tasks are checked individually only for file presence.
#
# Usage (submit FROM the 00_3_download folder):
#   sbatch submit_verify_integrity.slurm
# Chain it to the download array so it starts by itself when every task is done:
#   sbatch --dependency=afterok:<download_array_jobid> submit_verify_integrity.slurm

set -uo pipefail
# Exit on undefined variables and pipe failures; -e is intentionally omitted here as in verify_integrity.sh itself, so a single corrupted file does not abort the audit before the remaining files are checked.

hostname; date
echo "Job name: $SLURM_JOB_NAME"; echo "Job ID: $SLURM_JOB_ID"
# Log the node, timestamp, and job identifiers into the .out file for traceability.

export DATA_DIR="${DATA_DIR:-/users/genomics/albertoc/Tesi/hopes_and_dreams/datasets}"
# Data root: the value used for the thesis run, kept overridable with 'sbatch --export=ALL,DATA_DIR=...' rather than hardcoded inside verify_integrity.sh, which stays site-agnostic.

cd "$SLURM_SUBMIT_DIR"
# Submit with 'sbatch submit_verify_integrity.slurm' FROM the 00_3_download folder: SLURM copies this file to a spool directory, so "$0" cannot be used to locate verify_integrity.sh; SLURM_SUBMIT_DIR (the directory sbatch was called from) is the reliable anchor.

REPORT="$SLURM_SUBMIT_DIR/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.report"
bash verify_integrity.sh | tee "$REPORT"
# Run the three checks (gzip integrity, R1/R2 read lengths, 5' TSO on R2) and keep a copy of the output next to the log: the read-length and TSO sections are the evidence quoted in Materials & Methods, so they are worth a file of their own rather than only a job log.

corrupted=$(grep -c 'CORRUPTED' "$REPORT")
if [ "$corrupted" -gt 0 ]; then
    echo "[error] $corrupted corrupted file(s): re-download them by resubmitting their index of the download array (see README)." >&2
    exit 1
fi
echo "[ok] no corrupted file; report: $REPORT"
date
# verify_integrity.sh reports corruption on stdout but always exits 0, so the wrapper turns 'CORRUPTED' lines into a non-zero exit: that is what makes the job show up as FAILED in 'sacct' and what lets a --dependency=afterok chain (rename, Cell Ranger) refuse to start on a broken dataset.
