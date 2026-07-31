#!/bin/bash
#SBATCH --job-name=scran_norm
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=470G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

hostname; date
echo "Job name: $SLURM_JOB_NAME"; echo "Job ID: $SLURM_JOB_ID"
# Log the node, timestamp, and job identifiers into the .out file for traceability.

module load Miniconda3/202411
eval "$(/soft/system/software/Miniconda3/202411/bin/conda shell.bash hook)"
conda activate catalano_env
# Load Miniconda and activate catalano_env (needs rpy2/anndata2ri + R scran/scib stack, see environments/benchmark-py-r.yml).

export DATA_DIR="${DATA_DIR:-/users/genomics/albertoc/Tesi/hopes_and_dreams/datasets}"
export OMP_NUM_THREADS=90
# Data root: the same single DATA_DIR every phase uses (see datasets/README.md), kept
# overridable with 'sbatch --export=ALL,DATA_DIR=...'. FIG_DIR is deliberately NOT set:
# scran_norm.py then defaults to the repo's own 01_pre_processing/figures/01_3_normalization,
# so a cluster run fills the same versioned folder a local run does instead of dropping
# the two diagnostics somewhere under the data root, outside the clone.

which python3
python3 scran_norm.py