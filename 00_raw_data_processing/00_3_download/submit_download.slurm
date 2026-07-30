#!/bin/bash
#SBATCH --job-name=download_fastq
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --array=1-150%4
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

# submit_download.slurm
# SLURM wrapper for the FASTQ download (~2.4 TB), which runs on the cluster and not locally.
# One array task = one line of sample_map_gex_final.tsv = one SRR, so the 150 runs are downloaded as 150 independent jobs instead of one multi-day serial loop: a run that fails (bad transfer, wall clock, the empty SRR26541168 at line 44) takes down only its own task, and the failed indices can be resubmitted on their own with 'sbatch --array=<idx>[,<idx>...] submit_download.slurm'.
# The '%4' throttle keeps at most 4 tasks downloading at once: this step is network- and I/O-bound, not CPU-bound, and NCBI throttles aggressive parallel transfers, so more concurrency buys nothing and risks truncated files.
# Partition 'long' because a single run means 37-40 GB to transfer, extract and recompress: hours per task, more than 'normal' is meant for. A task killed on the partition limit costs nothing anyway - download_fastq.sh skips any SRR whose _2.fastq.gz is already there, so resubmitting that index resumes the array.
# This step needs OUTBOUND INTERNET from the compute nodes (prefetch talks to NCBI). If they are firewalled, the download has to run on a node that is not.
#
# Usage (submit FROM the 00_3_download folder):
#   sbatch submit_download.slurm                                    # all 150 runs, 4 at a time
#   sbatch --array=44,77 submit_download.slurm                      # only those lines of the sample map (re-run of failed/corrupted runs)
# The line of a given SRR, i.e. its array index:
#   grep -n SRR26540980 ../00_2_sample_mapping/sample_map_gex_final.tsv | cut -d: -f1
# Logs (one pair per task) are written next to this script as download_fastq_<jobid>_<taskid>.{out,err}.

set -euo pipefail
# Exit on error, undefined variables, and pipe failures, so a task that fails to download is reported as FAILED by SLURM rather than exiting 0 with nothing on disk.

hostname; date
echo "Job name: $SLURM_JOB_NAME"; echo "Array job: $SLURM_ARRAY_JOB_ID  Task: $SLURM_ARRAY_TASK_ID"
# Log the node, timestamp, and job identifiers into the .out file for traceability.

export DATA_DIR="${DATA_DIR:-/users/genomics/albertoc/Tesi/hopes_and_dreams/datasets}"
DOWNLOAD_ENV="${DOWNLOAD_ENV:-catalano_env}"
# Data root (must hold ~2.4TB): the value used for the thesis run, set here at the entry point and overridable with 'sbatch --export=ALL,DATA_DIR=...'. download_fastq.sh itself keeps requiring DATA_DIR with no default, so the terabytes can never land in a repo clone by accident.
# DOWNLOAD_ENV is the conda environment providing the SRA Toolkit (prefetch/fasterq-dump) and pigz: on the cluster that is catalano_env, built from environments/benchmark-hpc.yml, which declares sra-tools.

module load Miniconda3/202411
eval "$(/soft/system/software/Miniconda3/202411/bin/conda shell.bash hook)"
conda activate "$DOWNLOAD_ENV"
# Load Miniconda and activate the environment holding the SRA Toolkit.

which prefetch fasterq-dump pigz
# Sanity check: confirm the three binaries come from the activated environment before starting a multi-hour transfer.

cd "$SLURM_SUBMIT_DIR"
# Submit with 'sbatch submit_download.slurm' FROM the 00_3_download folder: SLURM copies this file to a spool directory, so "$0" does not point at the repo and cannot be used to locate download_fastq.sh; SLURM_SUBMIT_DIR (the directory sbatch was called from) is the reliable anchor.

SAMPLE_MAP="${SAMPLE_MAP:-$SLURM_SUBMIT_DIR/../00_2_sample_mapping/sample_map_gex_final.tsv}"
SRR=$(awk -v n="$SLURM_ARRAY_TASK_ID" 'NR==n {print $1}' "$SAMPLE_MAP")
[ -n "$SRR" ] || { echo "[error] no SRR on line $SLURM_ARRAY_TASK_ID of $SAMPLE_MAP" >&2; exit 1; }
echo ">>> [$SLURM_ARRAY_TASK_ID] $SRR"
# Map the array index to a run: line N of the sample map (the file has no header, one SRR per line in column 1). An index past the end of the file yields an empty SRR and fails immediately, instead of silently downloading nothing.

SHARD="$(mktemp)"
trap 'rm -f "$SHARD"' EXIT
printf '%s\n' "$SRR" > "$SHARD"
# download_fastq.sh loops over a sample map, so the single run of this task is handed to it as a one-line map on a node-local temp file. This reuses the download logic verbatim (prefetch --max-size 100G, fasterq-dump --include-technical, pigz) rather than duplicating it here, which is what keeps the local and cluster paths identical.

export TMPDIR="$DATA_DIR/fastq_raw"
mkdir -p "$TMPDIR"
# Scratch for prefetch/fasterq-dump: the intermediate .sra and the uncompressed FASTQ of a single run are tens of GB, far more than a node's /tmp usually holds, so point the temp space at the data volume that was sized for this step.

DATA_DIR="$DATA_DIR" SAMPLE_MAP="$SHARD" THREADS="${SLURM_CPUS_PER_TASK:-8}" \
    bash download_fastq.sh
# Run the download for this SRR alone. THREADS follows the allocation (fasterq-dump and pigz both use it), so changing --cpus-per-task above is enough to rescale.

FASTQ_DIR="$DATA_DIR/fastq_raw"
for r in 1 2; do
    f="$FASTQ_DIR/${SRR}_${r}.fastq.gz"
    [ -s "$f" ] || { echo "[error] missing or empty after download: $f" >&2; exit 1; }
done
echo "[ok] $SRR -> $FASTQ_DIR/${SRR}_[12].fastq.gz"
date
# Presence-and-size check so a task that produced nothing (or only R1) is marked FAILED and is easy to spot in 'sacct'. Content checks stay where they belong: verify_integrity.sh (gzip -t, read lengths, TSO) after the whole array has finished.
