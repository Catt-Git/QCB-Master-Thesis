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

export DATA_DIR=/users/genomics/albertoc/Tesi/hopes_and_dreams/datasets
export FIG_DIR=/users/genomics/albertoc/Tesi/hopes_and_dreams/scib_norm/figures/
export OMP_NUM_THREADS=90

which python3
python3 scran_norm.py