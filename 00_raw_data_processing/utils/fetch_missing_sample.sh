#!/bin/bash

# fetch_missing_sample.sh
# Re-download a SINGLE corrupted/missing SRR.
# Used for the two samples with a corrupted ENA copy: P07_A_P (SRR26540980) and P10_C_P (SRR26541023). Both have 151bp reads (instead of 101bp) but process correctly with --r1-length=26. Before relaunching Cell Ranger on these, delete the failed run folder: rm -rf cellranger_out/P07_A_P cellranger_out/P10_C_P.

# Usage: DATA_DIR=/path/with/space bash fetch_missing_sample.sh SRR26540980
# On the cluster prefer resubmitting that run's index of the download array, which does the same work inside an allocation and with the same code path as the original download:
#   idx=$(grep -n SRR26540980 ../00_2_sample_mapping/sample_map_gex_final.tsv | cut -d: -f1)
#   cd ../00_3_download && sbatch --array=$idx submit_download.slurm
# (delete the corrupted <SRR>_1/_2.fastq.gz first, otherwise the run is skipped as already present)

set -euo pipefail
# Exit on error, undefined variables, and pipe failures, essential for detecting data corruption or misconfiguration early in the pipeline.

SRR="${1:?no SRR given. Usage: DATA_DIR=<data root> bash fetch_missing_sample.sh SRR26540980}"
: "${DATA_DIR:?not set. Export DATA_DIR=<data root>; the run is re-downloaded into \$DATA_DIR/fastq_raw. See README.}"
FASTQ_DIR="${FASTQ_DIR:-$DATA_DIR/fastq_raw}"
# SRR accession passed as the first argument, and the folder holding the downloaded FASTQ (same convention as 00_3_download/download_fastq.sh, so the re-downloaded run lands next to the others).

mkdir -p "$FASTQ_DIR"
cd "$FASTQ_DIR"
# Work inside the data folder: prefetch caches the .sra in the cwd and fasterq-dump reads it back from there, and the replaced file must overwrite the corrupted one in place.

prefetch --max-size 100G "$SRR"
fasterq-dump --include-technical --split-files -p -e 8 -O . "$SRR"
pigz -p 8 "${SRR}_1.fastq" "${SRR}_2.fastq"
# Download and extract the run: --max-size 100G (the two files are 37-40GB, above prefetch's 20G default); --include-technical keeps the barcode read; then compress with pigz.

gzip -t "${SRR}_1.fastq.gz" && gzip -t "${SRR}_2.fastq.gz" && echo "$SRR OK" || echo "$SRR STILL CORRUPTED"
# Verify gzip integrity of both files.

echo "R1 length check (expected 101 or 151):"
zcat "${SRR}_1.fastq.gz" | head -400 | awk 'NR%4==2{print length($0)}' | sort -u
# Confirm the R1 read length to make sure the barcode read was retrieved.
