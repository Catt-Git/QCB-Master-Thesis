#!/bin/bash

# rename_for_cellranger.sh
# Rename the SRR FASTQ files into the format required by Cell Ranger: <sample>_S1_L001_R1_001.fastq.gz (barcode, was SRR_1) and <sample>_S1_L001_R2_001.fastq.gz (cDNA, was SRR_2). Files are COPIED (rsync), not moved, so the original SRR FASTQ files stay intact.
# Prerequisite: verify_integrity.sh must have confirmed that SRR_1 = barcode (26/101/151bp) and SRR_2 = cDNA. If fasterq-dump produced a different file ordering, fix the _1/_2 mapping below.

# INPUT: ../00_2_sample_mapping/sample_map_gex_final.tsv (col1=SRR, col6=sample_id) and $DATA_DIR/fastq_raw with the SRR*_1/_2.fastq.gz files
# OUTPUT: $DATA_DIR/renamed_fastq/ and $DATA_DIR/renamed_fastq/unique_samples.txt (input to Cell Ranger)
# Usage: DATA_DIR=/path/with/space bash rename_for_cellranger.sh
# On the cluster (how the thesis run was done) this is launched through 'sbatch submit_rename.slurm',
# which sets DATA_DIR and gives the copy its own allocation (~2.4TB of rsync).

set -euo pipefail
# Exit on error, undefined variables, and pipe failures, essential for detecting data corruption or misconfiguration early in the pipeline.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
: "${DATA_DIR:?not set. Export DATA_DIR=<data root>; reads \$DATA_DIR/fastq_raw, writes \$DATA_DIR/renamed_fastq. See README.}"
SAMPLE_MAP="${SAMPLE_MAP:-$SCRIPT_DIR/../00_2_sample_mapping/sample_map_gex_final.tsv}"
SOURCE_DIR="${SOURCE_DIR:-$DATA_DIR/fastq_raw}"
DEST_DIR="${DEST_DIR:-$DATA_DIR/renamed_fastq}"
RENAME_COMMANDS="$DATA_DIR/rename_commands.txt"
# Sample map (a small TSV in the repo, resolved from this script's location so any cwd works), folder holding the original SRR*_1/_2.fastq.gz files, and destination folder for the renamed copies. Since the copies double the footprint, DATA_DIR is required rather than defaulted so they cannot land in the repo clone by accident.

mkdir -p "$DEST_DIR"
# Create the destination folder if it does not exist.

awk -F'\t' -v src="$SOURCE_DIR" -v dest="$DEST_DIR" '{
    srr = $1
    sample = $6
    print src"/"srr"_1.fastq.gz\t"dest"/"sample"_S1_L001_R1_001.fastq.gz"
    print src"/"srr"_2.fastq.gz\t"dest"/"sample"_S1_L001_R2_001.fastq.gz"
}' "$SAMPLE_MAP" > "$RENAME_COMMANDS"
# Build a two-column (source, destination) list mapping each SRR file to its Cell Ranger-compliant name.

# For safety reasons, we do not rename in place; instead, we copy the files to a new folder. This avoids accidental overwriting of the original SRR files and allows for easy reversion if needed.
# Change this command if storage is a concern and you want to move instead of copy.

echo "Commands generated: $(wc -l < "$RENAME_COMMANDS") (expected: 2x number of samples)"
# Report how many copy operations were generated as a sanity check.

RSYNC_OPTS=(-ah)
[ -t 1 ] && RSYNC_OPTS+=(--progress)
# --progress only when stdout is a terminal: it redraws a percentage line with carriage returns, which is useful interactively but turns a SLURM .out file into megabytes of unreadable partial lines.

while IFS=$'\t' read -r original new; do
    echo "$(basename "$original") -> $(basename "$new")"
    rsync "${RSYNC_OPTS[@]}" "$original" "$new"
done < "$RENAME_COMMANDS"
# Copy each file to its new name with rsync (resume-safe), printing one line per file so the job log records what was copied.

cd "$DEST_DIR"
ls *_R1_001.fastq.gz | sed 's/_S1_L001_R1_001.fastq.gz//' | sort -u > unique_samples.txt
echo "Unique samples ready for Cell Ranger: $(wc -l < unique_samples.txt)"
# Derive the unique sample list from the renamed R1 files and report the count.
