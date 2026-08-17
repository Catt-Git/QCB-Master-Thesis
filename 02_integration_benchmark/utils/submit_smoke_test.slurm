#!/bin/bash
#SBATCH --job-name=smoke_metrics
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

hostname; date
echo "Job name: $SLURM_JOB_NAME"; echo "Job ID: $SLURM_JOB_ID"
# Log the node, timestamp and job identifiers into the .out file for traceability.

module load Miniconda3/202411
eval "$(/soft/system/software/Miniconda3/202411/bin/conda shell.bash hook)"
conda activate catalano_env
# Same environment used in 01_3: it needs the rpy2/anndata2ri stack plus R with
# scran (see environments/benchmark-hpc.yml).

export DATA_DIR=/users/genomics/albertoc/Tesi/hopes_and_dreams/datasets
export OMP_NUM_THREADS=8

# The fixture has to be copied here before launching (it is built locally, where
# shiao.h5ad lives):
#   scp datasets/smoke_fixture.h5ad shiva:$DATA_DIR/
# If it is missing, the script tries to rebuild it from shiao.h5ad, which needs
# far more memory and time than requested above.

which python3
cd "$SLURM_SUBMIT_DIR"
python3 smoke_test_metrics.py
