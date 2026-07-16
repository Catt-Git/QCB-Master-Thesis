#!/bin/bash
#SBATCH --job-name=h5ad_concat
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=450G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# submit_h5ad_concat.slurm
# SLURM wrapper that activates the scrna environment and runs build_combined_h5ad.py to merge the Cell Ranger .h5 outputs into a single AnnData.

hostname; date
echo "Job name: $SLURM_JOB_NAME"; echo "Job ID: $SLURM_JOB_ID"
# Log the node, timestamp, and job identifiers into the .out file for traceability.

module load Miniconda3/202411
eval "$(/soft/system/software/Miniconda3/202411/bin/conda shell.bash hook)"
conda activate scrna
# Load Miniconda and activate the scrna environment (cellranger module NOT needed here, the merge only reads the .h5 files).

which python3
# Sanity check: confirm python3 comes from the scrna environment.

export DATA_DIR="/users/genomics/albertoc/Tesi/hopes_and_dreams"
# Site-specific data root, set here at the entry point rather than defaulted inside the script. This is the value used for the thesis run; to replicate elsewhere, change this line (or export DATA_DIR before sbatch).

cd "$SLURM_SUBMIT_DIR"
# Submit with 'sbatch submit_h5ad_concat.slurm' FROM the 00_6_build_h5ad folder: SLURM copies this file to a spool directory, so "$0" cannot be used to locate build_combined_h5ad.py; SLURM_SUBMIT_DIR (the directory sbatch was called from) is the reliable anchor.

python3 build_combined_h5ad.py \
    --cellranger-dir "$DATA_DIR/cellranger_out" \
    --samples-file   "$DATA_DIR/renamed_fastq/unique_samples.txt" \
    --metadata-tsv   sample_metadata_final.tsv \
    --out            "$DATA_DIR/all_samples_combined.h5ad"
# Run the merge with python3 (NOT bash: bash would parse the Python syntax as shell commands and raise spurious errors, e.g. ImageMagick 'import').
# The Cell Ranger outputs, the sample list and the resulting .h5ad all live under DATA_DIR (outside the repo); only sample_metadata_final.tsv is a repo file, written next to this script by build_requested_metadata.py.
