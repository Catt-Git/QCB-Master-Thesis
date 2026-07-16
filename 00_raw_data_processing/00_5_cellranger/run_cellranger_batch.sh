#!/bin/bash

# run_cellranger_batch.sh
# Run 'cellranger count' on all samples in parallel (5 jobs x 16 cores = 80), inside a single SLURM allocation (see submit_cellranger.slurm).
# Key parameters (verbatim, for Methods): reference refdata-cellranger-GRCh38-3.0.0 (same as the paper); --r1-length=26 (R1 is 101/151bp = 16bp barcode + 10bp UMI + padding; without it Cell Ranger cannot detect the barcode); --chemistry=SC5P-R2 (10x Chromium 5'v1, confirmed by the TSO AAGCAGTGGTATCAAC at the start of R2); --localcores=16 --localmem=64; Cell Ranger 4.0.0 (via 'module load cellranger/4.0.0').

# INPUT: renamed_fastq/ and renamed_fastq/unique_samples.txt
# OUTPUT: cellranger_out/<sample>/outs/filtered_feature_bc_matrix.h5

set -uo pipefail
# Exit on undefined variables and pipe failures; note -e is intentionally omitted so one failing sample does not abort the whole batch (failures are logged instead).

: "${DATA_DIR:?not set. Export DATA_DIR=<data root> and REF=<path to refdata-cellranger-GRCh38-3.0.0>, or edit submit_cellranger.slurm. See README.}"
: "${REF:?not set. Export REF=<path to refdata-cellranger-GRCh38-3.0.0>, or edit submit_cellranger.slurm. See README.}"
FASTQ_DIR="${FASTQ_DIR:-$DATA_DIR/renamed_fastq}"
OUTDIR="${OUTDIR:-$DATA_DIR/cellranger_out}"
PARALLEL_JOBS="${PARALLEL_JOBS:-5}"
# Data root, input FASTQ folder, reference transcriptome, output directory, and number of samples to process concurrently (5 x 16 cores = 80 on a 90-core node).
# DATA_DIR and REF are deliberately REQUIRED rather than defaulted: the ~2.4TB of FASTQ live outside the repo and their location is site-specific, so a hardcoded default would only ever be correct on one machine and would fail late and cryptically everywhere else. Set them in submit_cellranger.slurm (the HPC entry point) or export them before running. The values used for the thesis are recorded in the README.

mkdir -p "$OUTDIR"
cd "$OUTDIR"
# Create and move into the output directory (Cell Ranger writes into the current working directory).

run_one_sample() {
    sample="$1"
    echo "=== Starting $sample ==="
    cellranger count --id="${sample}" \
      --transcriptome="$REF" \
      --fastqs="$FASTQ_DIR" \
      --sample="${sample}" \
      --r1-length=26 \
      --chemistry=SC5P-R2 \
      --localcores=16 --localmem=64

    if [ $? -ne 0 ]; then
        echo "$sample" >> "$OUTDIR/failed_samples.txt"
    fi
}
# Define the per-sample worker: run cellranger count with the fixed parameters and log the sample name to failed_samples.txt on non-zero exit. Note: Cell Ranger refuses to restart with an existing --id, so delete folders of failed runs before relaunching.

export -f run_one_sample
export FASTQ_DIR REF OUTDIR
# Export the function and variables so they are visible to the subshells spawned by xargs.

cat "$FASTQ_DIR/unique_samples.txt" \
  | xargs -P "$PARALLEL_JOBS" -I{} bash -c 'run_one_sample "$@"' _ {}
# Feed the sample list to xargs and run PARALLEL_JOBS workers concurrently, one per sample.

echo "Done. Check failed_samples.txt for any errors."
# Point the user to the failure log after the batch completes.
