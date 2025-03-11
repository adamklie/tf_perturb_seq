#!/bin/bash
#SBATCH --partition=carter-gpu
#SBATCH --account=carter-gpu
#SBATCH --gpus=1
#SBATCH --output=/cellar/users/aklie/data/datasets/tf_perturb_seq/slurm_logs/%x.%A.out
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=14-00:00:00

#####
# USAGE:
# sbatch pipeline.sh
#####

# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configure env
source activate /cellar/users/aklie/opt/miniconda3/envs/rapids_singlecell
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/lib/

# Input file
input_h5ads=(
    /cellar/users/aklie/data/datasets/tf_perturb_seq/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/Preliminary_Analysis/qc.h5ad
    /cellar/users/aklie/data/datasets/tf_perturb_seq/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/Preliminary_Analysis/qc.h5ad
    /cellar/users/aklie/data/datasets/tf_perturb_seq/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/Preliminary_Analysis/qc.h5ad
)
names=(
    Hon_WTC11-cardiomyocyte-differentiation
    Huangfu_HUES8-definitive-endoderm-differentiation
    Huangfu_HUES8-embryonic-stemcell-differentiation
)

# 1. merge.py
script=/cellar/users/aklie/data/datasets/tf_perturb_seq/bin/merge.py
output_h5ad=/cellar/users/aklie/data/datasets/tf_perturb_seq/merge.h5ad
logfile=/cellar/users/aklie/data/datasets/tf_perturb_seq/merge.log
cmd="python $script --inputs ${input_h5ads[@]} --names ${names[@]} --names_key dataset --make_bcs_unique --output $output_h5ad --log_file $logfile"
echo $cmd
#eval $cmd

# 2. preprocess.py
script=/cellar/users/aklie/data/datasets/tf_perturb_seq/bin/preprocess.py
input_h5ad=/cellar/users/aklie/data/datasets/tf_perturb_seq/merge.h5ad
output_h5ad=/cellar/users/aklie/data/datasets/tf_perturb_seq/preprocess.h5ad
plot_dir=/cellar/users/aklie/data/datasets/tf_perturb_seq/plots/preprocess
logfile=/cellar/users/aklie/data/datasets/tf_perturb_seq/preprocess.log
cmd="python $script --input $input_h5ad --output $output_h5ad --target_sum 10000 --n_top_genes 5000 --max_value 10 --plot_dir $plot_dir --log_file $logfile"
echo $cmd
#eval $cmd

# 4. cluster
script=/cellar/users/aklie/data/datasets/tf_perturb_seq/bin/cluster.py
input_h5ad=/cellar/users/aklie/data/datasets/tf_perturb_seq/preprocess.h5ad
output_h5ad=/cellar/users/aklie/data/datasets/tf_perturb_seq/cluster.h5ad
plot_dir=/cellar/users/aklie/data/datasets/tf_perturb_seq/plots/cluster
logfile=/cellar/users/aklie/data/datasets/tf_perturb_seq/cluster.log
cmd="python $script --input $input_h5ad --output $output_h5ad --n_neighbors 15 --n_pcs 50 --integrate dataset --resolution 0.6 --plot_dir $plot_dir --log_file $logfile"
echo $cmd
eval $cmd

# date
date
