#!/bin/bash
#SBATCH --partition=carter-gpu
#SBATCH --account=carter-gpu
#SBATCH --gpus=1
#SBATCH --output=/cellar/users/aklie/data/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/bin/slurm_logs/%x.%A.out
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
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
input_h5=/cellar/users/aklie/data/datasets/tf_perturb_seq/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/Cell_Ranger_Output/filtered_feature_bc_matrix.h5

# 1. make_anndata.py
script=/cellar/users/aklie/data/datasets/tf_perturb_seq/bin/make_anndata.py
output_h5ad=/cellar/users/aklie/data/datasets/tf_perturb_seq/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/Preliminary_Analysis/anndata.h5ad
logfile=/cellar/users/aklie/data/datasets/tf_perturb_seq/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/Preliminary_Analysis/anndata.log
cmd="python $script --input $input_h5 --output $output_h5ad --log_file $logfile"
echo $cmd
eval $cmd

# 2. qc.py
script=/cellar/users/aklie/data/datasets/tf_perturb_seq/bin/qc.py
input_h5ad=/cellar/users/aklie/data/datasets/tf_perturb_seq/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/Preliminary_Analysis/anndata.h5ad
output_h5ad=/cellar/users/aklie/data/datasets/tf_perturb_seq/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/Preliminary_Analysis/qc.h5ad
plot_dir=/cellar/users/aklie/data/datasets/tf_perturb_seq/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/Preliminary_Analysis/plots/qc
logfile=/cellar/users/aklie/data/datasets/tf_perturb_seq/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/Preliminary_Analysis/qc.log
cmd="python $script --input $input_h5ad --output $output_h5ad --max_genes 5000 --max_mt_pct 20 --min_gene_count 3 --plot_dir $plot_dir --log_file $logfile"
echo $cmd
eval $cmd

# 3. preprocess.py
script=/cellar/users/aklie/data/datasets/tf_perturb_seq/bin/preprocess.py
input_h5ad=/cellar/users/aklie/data/datasets/tf_perturb_seq/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/Preliminary_Analysis/qc.h5ad
output_h5ad=/cellar/users/aklie/data/datasets/tf_perturb_seq/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/Preliminary_Analysis/preprocess.h5ad
plot_dir=/cellar/users/aklie/data/datasets/tf_perturb_seq/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/Preliminary_Analysis/plots/preprocess
logfile=/cellar/users/aklie/data/datasets/tf_perturb_seq/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/Preliminary_Analysis/preprocess.log
cmd="python $script --input $input_h5ad --output $output_h5ad --target_sum 10000 --n_top_genes 5000 --max_value 10 --plot_dir $plot_dir --log_file $logfile"
echo $cmd
eval $cmd

# 4. cluster
script=/cellar/users/aklie/data/datasets/tf_perturb_seq/bin/cluster.py
input_h5ad=/cellar/users/aklie/data/datasets/tf_perturb_seq/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/Preliminary_Analysis/preprocess.h5ad
output_h5ad=/cellar/users/aklie/data/datasets/tf_perturb_seq/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/Preliminary_Analysis/cluster.h5ad
plot_dir=/cellar/users/aklie/data/datasets/tf_perturb_seq/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/Preliminary_Analysis/plots/cluster
logfile=/cellar/users/aklie/data/datasets/tf_perturb_seq/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/Preliminary_Analysis/cluster.log
cmd="python $script --input $input_h5ad --output $output_h5ad --n_neighbors 15 --n_pcs 50 --resolution 0.6 --plot_dir $plot_dir --log_file $logfile"
echo $cmd
eval $cmd

# date
date
