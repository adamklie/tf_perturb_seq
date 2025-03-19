#!/bin/bash
#SBATCH --partition=carter-gpu
#SBATCH --account=carter-gpu
#SBATCH --gpus=1
#SBATCH --output=/cellar/users/aklie/data/datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/bin/slurm_logs/%x.%A.out
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
input_h5=/cellar/users/aklie/data/datasets/tf_perturb_seq/datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/Cell_Ranger_Output/filtered_feature_bc_matrix.h5
scripts_dir=/cellar/users/aklie/data/datasets/tf_perturb_seq/scripts
output_dir=/cellar/users/aklie/data/datasets/tf_perturb_seq/datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/Transcriptome_Analysis

# 1. make_anndata.py
script=$scripts_dir/make_anndata.py
output_h5ad=$output_dir/anndata.h5ad
logfile=$output_dir/anndata.log
cmd="python $script --input $input_h5 --output $output_h5ad --log_file $logfile"
echo $cmd
eval $cmd

# 2. qc.py
script=$scripts_dir/qc.py
input_h5ad=$output_dir/anndata.h5ad
output_h5ad=$output_dir/qc.h5ad
plot_dir=$output_dir/plots/qc
logfile=$output_dir/qc.log
cmd="python $script --input $input_h5ad --output $output_h5ad --min_genes 500 --max_mt_pct 20 --min_gene_count 3 --plot_dir $plot_dir --log_file $logfile"
echo $cmd
eval $cmd

# 3. preprocess.py
script=$scripts_dir/preprocess.py
input_h5ad=$output_dir/qc.h5ad
output_h5ad=$output_dir/preprocess.h5ad
plot_dir=$output_dir/plots/preprocess
logfile=$output_dir/preprocess.log
cmd="python $script --input $input_h5ad --output $output_h5ad --target_sum 10000 --n_top_genes 5000 --max_value 10 --plot_dir $plot_dir --log_file $logfile"
echo $cmd
eval $cmd

# 4. cluster
script=$scripts_dir/cluster.py
input_h5ad=$output_dir/preprocess.h5ad
output_h5ad=$output_dir/cluster.h5ad
plot_dir=$output_dir/plots/cluster
logfile=$output_dir/cluster.log
cmd="python $script --input $input_h5ad --output $output_h5ad --n_neighbors 25 --n_pcs 50 --resolution 1 --plot_dir $plot_dir --log_file $logfile"
echo $cmd
eval $cmd

# date
date
