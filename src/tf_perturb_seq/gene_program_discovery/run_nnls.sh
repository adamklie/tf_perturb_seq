#!/usr/bin/env bash
# submit_nnls.sh — Submit the NNLS projection job to SLURM.
#
# Edit the variables in the CONFIG section, then run:
#   bash submit_nnls.sh

# -- SLURM directives ---------------------------------------------------------
#SBATCH --job-name=nnls_projection
#SBATCH --output=/hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/NNLS/logs/nnls_%j.out       # stdout  (%j = job ID)
#SBATCH --error=/hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/NNLS/logs/nnls_%j.err        # stderr
#SBATCH --time=24:00:00                 # wall-clock limit (HH:MM:SS)
#SBATCH --mem=164G                       # memory per node
#SBATCH --cpus-per-task=8              
#SBATCH --partition=common           

# -- CONFIG — edit these ------------------------------------------------------
INPATH="/hpc/group/gersbachlab/seg95/"
OUTPATH="/hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/NNLS/"
SPECTRA_FILE="tf_perturb_seq/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/PerturbNMF/Result/051026_gersbach_iHep_torch_minibatch/Inference/Inference.gene_spectra_tpm.k_80.dt_2_0.txt"             # relative to INPATH
REF_FILE="helen_data/HumanLiverSeurat_QC_donorBatchC.h5ad"             # relative to INPATH
OUTFILE="usage_matrix_donorBatchC.csv"            # written to OUTPATH

PYTHON_SCRIPT="/hpc/group/gersbachlab/seg95/tf_perturb_seq/src/tf_perturb_seq/gene_program_discovery/run_nnls.py"  # full path to the python script
CONDA_ENV="scanpy_env"            

# Optional NMF settings (comment out to use script defaults)
MAX_ITER=1000
TOL=1e-4

# -- Setup --------------------------------------------------------------------
set -euo pipefail
mkdir -p logs

# Activate environment (adjust for your HPC module system)
# Option A — conda / mamba
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

# Option B — module system (uncomment and adjust if not using conda)
# module load python/3.10 cuda/11.8

# -- Run ----------------------------------------------------------------------
echo "Job started: $(date)"
echo "Node: $(hostname)"

python "${PYTHON_SCRIPT}" \
    --inpath   "${INPATH}"  \
    --outpath  "${OUTPATH}" \
    --spectra  "${SPECTRA_FILE}" \
    --refdata  "${REF_FILE}" \
    --outfile  "${OUTFILE}" \
    --max-iter "${MAX_ITER}" \
    --tol      "${TOL}"

echo "Job finished: $(date)"