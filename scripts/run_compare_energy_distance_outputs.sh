#!/bin/bash
#SBATCH --job-name=edist_compare
#SBATCH --output=logs/edist_compare_%j.out
#SBATCH --error=logs/edist_compare_%j.err
#SBATCH --time=1:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1

# ── Configuration ─────────────────────────────────────────────────────────────
# Conda environment name
CONDA_ENV="torch-cNMF"

# Path to the comparison script
SCRIPT="/hpc/group/gersbachlab/seg95/tf_perturb_seq/scripts/compare_energy_distance_outputs.py"

# Root directory containing all dataset folders
DATA_ROOT="/hpc/group/gersbachlab/seg95/tf_perturb_seq/src/tf_perturb_seq/energy_distance"

# Dataset subdirectory names and their short labels (must be same length & order)
DATASETS=(
    data_cardio
    data_de
    data_stem
)
LABELS=(
    cardio
    de
    stem
)

# Output directory
OUTDIR="${DATA_ROOT}/comparison_outputs"

# pval_mean significance cutoff
PVAL_CUTOFF=0.05

# ── Derived args (no edits needed below) ─────────────────────────────────────
DATASET_PATHS=()
for ds in "${DATASETS[@]}"; do
    DATASET_PATHS+=("${DATA_ROOT}/${ds}")
done

# ── Environment ───────────────────────────────────────────────────────────────
mkdir -p logs
mkdir -p "${OUTDIR}"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

echo "=============================================="
echo "  Job:      ${SLURM_JOB_ID}"
echo "  Node:     ${SLURM_NODELIST}"
echo "  Started:  $(date)"
echo "  Datasets: ${DATASETS[*]}"
echo "  Outdir:   ${OUTDIR}"
echo "=============================================="

# ── Run ───────────────────────────────────────────────────────────────────────
python "${SCRIPT}" \
    --datasets "${DATASET_PATHS[@]}" \
    --labels   "${LABELS[@]}" \
    --outdir   "${OUTDIR}" \
    --pval-cutoff "${PVAL_CUTOFF}"

EXIT_CODE=$?

echo "=============================================="
echo "  Finished: $(date)"
echo "  Exit code: ${EXIT_CODE}"
echo "=============================================="

exit ${EXIT_CODE}