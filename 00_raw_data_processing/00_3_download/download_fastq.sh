#!/bin/bash

# download_fastq.sh
# Download the FASTQ files of the human GEX libraries via the SRA Toolkit.
# Total file size is about 2.4 TB! :(
# Why not direct ENA (wget/aria2c): the data were sequenced 151x151bp (non-standard) and released recently; ENA had not processed the runs correctly and EXCLUDED the barcode read (R1) as a "technical read", so bulk ENA downloads were missing the barcode (or corrupted/truncated). The reliable method is prefetch + fasterq-dump with --include-technical.

# INPUT: ../00_2_sample_mapping/sample_map_gex_final.tsv (col 1 = SRR)
# OUTPUT: <SRR>_1.fastq.gz (barcode+UMI, 26bp effective) and <SRR>_2.fastq.gz (cDNA) in $DATA_DIR/fastq_raw
# REQUIREMENTS: sra-tools (prefetch, fasterq-dump), pigz, and ~2.4TB free under $DATA_DIR
# Usage: DATA_DIR=/path/with/space bash download_fastq.sh

set -euo pipefail
# Exit on error, undefined variables, and pipe failures, essential for detecting data corruption or misconfiguration early in the pipeline.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
: "${DATA_DIR:?not set. Export DATA_DIR=<data root with room for ~2.4TB>; the FASTQ are written to \$DATA_DIR/fastq_raw. See README.}"
SAMPLE_MAP="${SAMPLE_MAP:-$SCRIPT_DIR/../00_2_sample_mapping/sample_map_gex_final.tsv}"
FASTQ_DIR="${FASTQ_DIR:-$DATA_DIR/fastq_raw}"
THREADS="${THREADS:-8}"
# Sample map providing the SRR list (a small TSV in the repo, resolved from this script's location so any cwd works), destination for the FASTQ, and number of threads for fasterq-dump/pigz.
# DATA_DIR is required, not defaulted: ~2.4TB must not land inside the repo clone by accident, and the right location is site-specific.

mkdir -p "$FASTQ_DIR"
cd "$FASTQ_DIR"
# Work inside the data folder: prefetch writes its .sra cache into the cwd and fasterq-dump then reads it from there, so the download runs where the space is rather than wherever the script was launched from.

cut -f1 "$SAMPLE_MAP" | while read -r srr; do
    [ -z "$srr" ] && continue
    if [ -f "${srr}_2.fastq.gz" ]; then
        echo "[SKIP] already present: $srr"
        continue
    fi
    echo "=== Downloading $srr ==="
    prefetch --max-size 100G "$srr"
    fasterq-dump --include-technical --split-files -p -e "$THREADS" -O . "$srr"
    pigz -p "$THREADS" "${srr}_1.fastq" "${srr}_2.fastq"
    echo "=== $srr done ==="
done
# For each SRR: skip if already downloaded; prefetch with --max-size 100G (the 20G default skips the 37-40GB files); extract with --include-technical (critical: forces the barcode read R1 to be written, otherwise dropped) and --split-files; then compress with pigz.

echo "Download complete. Run verify_integrity.sh for the checks."
# Remind the user to validate the downloaded files before renaming.
# There should be 2 FASTQ files per SRR, since there are 150 SRR - 1 (empty) = 149 SRR with 2 FASTQ each = 298 FASTQ files in total.
