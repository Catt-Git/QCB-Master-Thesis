#!/bin/bash

# check_cellranger_completeness.sh
# Verify that EVERY sample produced a Cell Ranger filtered matrix (outs/filtered_feature_bc_matrix.h5).
# Completeness check for the Cell Ranger stage: lists samples whose run did not produce the expected .h5, and cross-references cellranger_out/failed_samples.txt. Run this before build_combined_h5ad.py so the merge does not silently skip samples.

# INPUT: unique_samples.txt (one sample_id per line) and the Cell Ranger output directory
# Usage: bash check_cellranger_completeness.sh [SAMPLES_FILE] [CELLRANGER_DIR]
# Exit code: 0 if all present, 1 if any matrix is missing (so it can be chained with &&).

set -uo pipefail
# Exit on undefined variables and pipe failures; -e omitted so the loop reports ALL missing samples, not just the first.

if [ $# -lt 2 ]; then
    : "${DATA_DIR:?not set. Export DATA_DIR=<data root, holding renamed_fastq/ and cellranger_out/>, or pass both paths as arguments. See README.}"
fi
SAMPLES_FILE="${1:-${DATA_DIR:-}/renamed_fastq/unique_samples.txt}"
CELLRANGER_DIR="${2:-${DATA_DIR:-}/cellranger_out}"
# Sample list and Cell Ranger output directory. The data live outside the repo, so these derive from DATA_DIR using the same convention as 00_5_cellranger/run_cellranger_batch.sh; DATA_DIR is required only when the paths are not passed explicitly as arguments.

expected=0
missing=0
missing_list=()
# Counters and list of samples without a filtered matrix.

while read -r sample; do
    [ -z "$sample" ] && continue
    expected=$((expected + 1))
    h5="$CELLRANGER_DIR/$sample/outs/filtered_feature_bc_matrix.h5"
    if [ ! -f "$h5" ]; then
        missing=$((missing + 1))
        missing_list+=("$sample")
    fi
done < "$SAMPLES_FILE"
# Walk the expected sample list and flag any sample missing its filtered matrix.

present=$((expected - missing))
echo "Expected samples (from $(basename "$SAMPLES_FILE")): $expected"
echo "With filtered_feature_bc_matrix.h5:               $present"
echo "Missing:                                          $missing"
# Report the expected/present/missing counts.

if [ -f "$CELLRANGER_DIR/failed_samples.txt" ]; then
    echo "Note: failed_samples.txt exists, logged failures:"
    sed 's/^/  /' "$CELLRANGER_DIR/failed_samples.txt"
fi
# Surface any failures already logged by the Cell Ranger batch runner.

if [ "$missing" -gt 0 ]; then
    echo "Samples without output (delete the folder and re-run cellranger count):"
    printf '  %s\n' "${missing_list[@]}"
    exit 1
fi
echo "All expected samples produced a filtered matrix."
# List the missing samples and fail (exit 1) if any; otherwise confirm success.
