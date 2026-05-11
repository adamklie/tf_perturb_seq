#!/bin/bash
# Production cNMF launch for Huangfu HUES8 Embryonic Stem Cell.
# Mirrors Hon's `030726_20iter_5KHVG_torch_halsvar_batch_e7` params exactly
# (halsvar / batch / 20 iter / 5K HVG / density thresholds 0.2 + 2.0 / categorical_key=batch / benchmark k list).
# Same scaffold as the HTv2 testbed (job 10577039) — see
# datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/PerturbNMF/Script/torch-cNMF_batch.sh.
# DO NOT SUBMIT until the HTv2 testbed completes cleanly and the structure has been verified.

#SBATCH --job-name=050926_HuangfuESC_20iter_5KHVG_torch_halsvar_batch
#SBATCH --output=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/PerturbNMF/Result/050926_HuangfuESC_20iter_5KHVG_torch_halsvar_batch/logs/%j.out
#SBATCH --error=/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/PerturbNMF/Result/050926_HuangfuESC_20iter_5KHVG_torch_halsvar_batch/logs/%j.err
#SBATCH --partition=carter-gpu
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --gres=gpu:1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=aklie@ucsd.edu

set -euo pipefail

START_TIME=$(date +%s)
echo "Job: ${SLURM_JOB_ID:-N/A} @ $(hostname) $(date)"
echo "Partition: ${SLURM_JOB_PARTITION:-N/A}"

OUT_DIR="/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/PerturbNMF/Result"
RUN_NAME="050926_HuangfuESC_20iter_5KHVG_torch_halsvar_batch"
LOG_DIR="$OUT_DIR/$RUN_NAME"
DATA_H5AD="/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/PerturbNMF/Data/ESC_sceptre_v1_perturbnmf.h5ad"
PIPELINE_PY="/cellar/users/aklie/projects/tf_perturb_seq/external/PerturbNMF/src/Stage1_Inference/torch-cNMF/Slurm_Version/torch_cnmf_inference_pipeline.py"
VENV_PY=/cellar/users/aklie/projects/tf_perturb_seq/.venv/bin/python

mkdir -p "$LOG_DIR/logs"

export PYTHONPATH="/cellar/users/aklie/projects/tf_perturb_seq/external/PerturbNMF/src:${PYTHONPATH:-}"

echo "Python: $VENV_PY"
$VENV_PY --version
nvidia-smi 2>/dev/null || echo "nvidia-smi not available"

MONITOR_LOG="$LOG_DIR/logs/resource_monitor_${SLURM_JOB_ID:-na}.log"
monitor_resources() {
    while true; do
        {
            date '+%Y-%m-%d %H:%M:%S'
            echo "=== Memory ==="; free -h
            echo "=== GPU ==="
            nvidia-smi --query-gpu=timestamp,name,utilization.gpu,utilization.memory,memory.used,memory.free,temperature.gpu --format=csv 2>/dev/null || echo "n/a"
            echo "---"
        } >> "$MONITOR_LOG"
        sleep 30
    done
}
monitor_resources &
MONITOR_PID=$!

$VENV_PY -u "$PIPELINE_PY" \
    --counts_fn "$DATA_H5AD" \
    --output_directory "$OUT_DIR" \
    --run_name "$RUN_NAME" \
    --species human \
    --algo halsvar \
    --mode batch \
    --init random \
    --tol 1e-7 \
    --use_gpu \
    --batch_max_epoch 1000 \
    --batch_hals_max_iter 1000 \
    --batch_hals_tol 0.005 \
    --numiter 20 \
    --numhvgenes 5000 \
    --sel_thresh 0.2 2.0 \
    --categorical_key batch \
    --gene_names_key symbol \
    --K 5 6 7 8 9 10 11 12 13 14 15 17 19 21 23 25 27 30 35 40 45 50 55 60 70 80 90 100 150 200 \
    --run_factorize --run_refit --run_compile_annotation --run_diagnostic_plots

kill $MONITOR_PID 2>/dev/null || true

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
echo "Total: ${DURATION}s ($((DURATION/3600))h $((DURATION%3600/60))m $((DURATION%60))s)"
echo "Done @ $(date)"
