#!/usr/bin/bash
#
#SBATCH -J clone_removal      # Job name
#SBATCH -N 1                          # Total number of nodes requested (16 cores/node)
#SBATCH -t 24:00:00                   # Run time (hh:mm:ss) - 20 hrs limit
#SBATCH -p 512GB
#SBATCH -o ./log/clone_removal.out
#SBATCH -e ./log/clone_removal.err

#Note: Use singularity 4.1.0 for this code. for UTSW environment, run "module load singularity/4.1.0" to activate singularity. 
module load singularity/4.1.0

#Define the path to the config file and bin directory
CONTAINER_PATH="/project/GCRB/Hon_lab/s223695/Data_project/Perturb_seq_edist_pipeline/container/edist_pipeline_v01.sif"
CONFIG_PATH="/project/GCRB/Hon_lab/s223695/Data_project/20250606_TF_neuro_reanalysis/sanity_check_pipeline/config.json"
BIN_PATH="/project/GCRB/Hon_lab/s223695/Data_project/20250606_TF_neuro_reanalysis/sanity_check_pipeline/bin"

#Note: When you run this code, please make sure that /pipeline_output/annotation_file_table.csv (or a defined name in config file) exists.

echo "Filtering clonal cells"
singularity exec --nv ${CONTAINER_PATH} python ${BIN_PATH}/clone_removal.py ${CONFIG_PATH}

echo "All steps completed successfully."