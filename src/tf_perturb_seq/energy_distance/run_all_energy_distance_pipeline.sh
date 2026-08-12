#!/usr/bin/bash
#SBATCH -J edist_pipeline_step1_2      # Job name
#SBATCH -N 1                          # Total number of nodes requested (16 cores/node)
#SBATCH -t 2:00:00                   # Run time (hh:mm:ss) - 7 day limit
#SBATCH --mem=100G
#SBATCH -p scavenger-gpu
#SBATCH -o run_output_wtc11_hon.out
#SBATCH -e run_output_wtc11_hon.err
#SBATCH --gres=gpu:6000_ada

set -e 

#TARGET_FILE_ID="syn70753570"
#SYNAPSE_TOKEN=${1}

#MUDATA_PATH="/hpc/group/gersbachlab/seg95/CRISPR_Pipeline/helen_output_poolabcd_prod_v10/inference_mudata.h5mu"
#MUDATA_PATH="/hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq_2026_04_09_outs_entertaining_hamster_inference_mudata.h5mu"
#MUDATA_PATH="/hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq_2026_04_13_outs_sceptre_v1_inference_mudata.h5mu"
#MUDATA_PATH="/hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/inference_mudata.h5mu"

# BRIDGE DATASETS
#MUDATA_PATH="/hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/sceptre/inference_mudata.h5mu"
#MUDATA_PATH="/hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/sceptre/inference_mudata.h5mu"
#MUDATA_PATH="/hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Hon_WTC11-benchmark_TF-Perturb-seq/sceptre/inference_mudata.h5mu"
#MUDATA_PATH="/hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Hon_WTC11-benchmark_TF-Perturb-seq/sceptre/inference_mudata.h5mu"
#MUDATA_PATH="/hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Huangfu_WTC11-benchmark_TF-Perturb-seq/sceptre/inference_mudata.h5mu"

CONTAINER_PATH="/hpc/group/gersbachlab/seg95/tf_perturb_seq/src/tf_perturb_seq/energy_distance/edist_pipeline_v2.sif"
CONFIG12_PATH="/hpc/group/gersbachlab/seg95/tf_perturb_seq/src/tf_perturb_seq/energy_distance/data_wtc11_hon/config1_2.json"
CONFIG3_PATH="/hpc/group/gersbachlab/seg95/tf_perturb_seq/src/tf_perturb_seq/energy_distance/data_wtc11_hon/config3.json"
BIN_PATH="/hpc/group/gersbachlab/seg95/tf_perturb_seq/src/tf_perturb_seq/energy_distance/energy_dist_pipeline/bin"
SCRIPT_DIR="/hpc/group/gersbachlab/seg95/tf_perturb_seq/src/tf_perturb_seq/energy_distance"

cd ${SCRIPT_DIR}

export PYTHONPATH="/tmp/muon_deps:/hpc/group/gersbachlab/seg95/tf_perturb_seq/src/tf_perturb_seq/energy_distance/energy_dist_pipeline"

#If data folder doesn't exist, prepare data folder
mkdir -p "./data_wtc11_hon"

# Pull container if not already present
#if [ ! -f "${CONTAINER_PATH}" ]; then
#    echo "Pulling container..."
#    apptainer pull ${CONTAINER_PATH} docker://docker.io/takechikara/energy_distance_env:latest
#fi
#chmod +x edist_pipeline_v2.sif

# Clone pipeline repo if not already present
if [ ! -d "./energy_dist_pipeline" ]; then
    echo "Cloning git repo..."
    git clone https://github.com/Chikara-Takeuchi/energy_dist_pipeline.git
fi

echo "Installing muon..."
apptainer exec --nv ${CONTAINER_PATH} pip install --target=/tmp/muon_deps muon

echo "Preprocess mudata file"
#apptainer exec --nv ${CONTAINER_PATH} python preprocess_mudata.py \
#    --synapse_id ${TARGET_FILE_ID} \
#    --auth_token ${SYNAPSE_TOKEN} \
#    --config1_2_path ${CONFIG12_PATH} \
#    --config3_path ${CONFIG3_PATH}

if [ ! -f "./data_wtc11_hon/pca_dataframe.pickle" ] || [ ! -f "./data_wtc11_hon/gRNA_dict.pickle" ]; then
    echo "Preprocess mudata file"
    apptainer exec --nv --bind $(dirname ${MUDATA_PATH}):$(dirname ${MUDATA_PATH}) ${CONTAINER_PATH} \
        bash -c "PYTHONPATH=/tmp/muon_deps python preprocess_mudata_local.py --mudata_path ${MUDATA_PATH} --config1_2_path ${CONFIG12_PATH} --config3_path ${CONFIG3_PATH}"
else
    echo "Preprocessed files already exist, skipping..."
fi

echo "[Step1] Filtering outlier gRNAs"
apptainer exec --nv --bind /hpc/group:/hpc/group ${CONTAINER_PATH} bash -c "PYTHONPATH=/tmp/muon_deps python ${BIN_PATH}/1_filtering_gRNA.py ${CONFIG12_PATH}"

echo "[Step2] calculate energy distance between targets and non-targeting"
apptainer exec --nv --bind /hpc/group:/hpc/group ${CONTAINER_PATH} bash -c "PYTHONPATH=/tmp/muon_deps python ${BIN_PATH}/2_e_distance_nontargeting.py ${CONFIG12_PATH}"

echo "[Step2_1] visualize results of energy distance analysis"
apptainer exec --nv --bind /hpc/group:/hpc/group ${CONTAINER_PATH} bash -c "PYTHONPATH=/tmp/muon_deps python ${BIN_PATH}/2_1_Plot_figure.py ${CONFIG12_PATH}"

echo "All steps completed successfully."
