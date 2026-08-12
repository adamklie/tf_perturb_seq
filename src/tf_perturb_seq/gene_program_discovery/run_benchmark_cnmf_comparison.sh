#!/bin/bash
#SBATCH --job-name=cnmf_comparison
#SBATCH --output=logs/cnmf_comparison_%j.out
#SBATCH --error=logs/cnmf_comparison_%j.err
#SBATCH --time=4:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4
#SBATCH --partition=common

set -euo pipefail

# - Paths -
SCRIPT_DIR="/hpc/group/gersbachlab/seg95/tf_perturb_seq/src/tf_perturb_seq/gene_program_discovery"
SCRIPT="${SCRIPT_DIR}/benchmark_cnmf_comparison.py"
LOG_DIR="${SCRIPT_DIR}/logs"

mkdir -p "${LOG_DIR}"

# - Environment -
# conda's activation scripts (e.g. activate.d/env_vars.sh) reference variables
# like LD_LIBRARY_PATH without a default, which trips `set -u` if they're not
# already set in the job's environment. Disable nounset just for activation.
set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate torch-cNMF
set -u

echo "============================================"
echo "Job:        ${SLURM_JOB_ID:-local}"
echo "Node:       $(hostname)"
echo "Started:    $(date)"
echo "Python:     $(which python)"
echo "Script:     ${SCRIPT}"
echo "============================================"

# - Run -
python "${SCRIPT}"

echo "--------------------------------------------"
echo "Finished:   $(date)"
echo "============================================"