#!/bin/bash
#SBATCH --job-name=1_integrate_gene_technology-benchmark_WTC11_TF-Perturb-seq
#SBATCH --partition=carter-compute
#SBATCH --cpus-per-task=4
#SBATCH --mem=300G
#SBATCH --time=14-00:00:00
#SBATCH --output=/cellar/users/aklie/data/datasets/tf_perturb_seq/datasets/technology-benchmark_WTC11_TF-Perturb-seq/bin/slurm_logs/%x.%A.out

#####
# USAGE:
# sbatch 1_integrate_gene.sh
#####

# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
#source activate /cellar/users/aklie/opt/miniconda3/envs/cellcommander
source activate /cellar/users/aklie/opt/miniconda3/envs/scverse-py39-R431
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/lib/

# Inputs
input_files=(
    /cellar/users/aklie/data/datasets/tf_perturb_seq/datasets/Hon_WTC11-benchmark_TF-Perturb-seq/results/2_qc/2025_12_01/mapping_qc/gene.filtered.h5ad
    /cellar/users/aklie/data/datasets/tf_perturb_seq/datasets/Huangfu_WTC11-benchmark_TF-Perturb-seq/results/2_qc/2025_12_01/mapping_qc/gene.filtered.h5ad
    /cellar/users/aklie/data/datasets/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/results/2_qc/2025_12_01/mapping_qc/gene.filtered.h5ad
    /cellar/users/aklie/data/datasets/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/results/2_qc/2025_12_01/mapping_qc/gene.filtered.h5ad
)
sample_ids=(
    Hon_WTC11-benchmark_TF-Perturb-seq
    Huangfu_WTC11-benchmark_TF-Perturb-seq
    Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3
    Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2
)
outdir_path=/cellar/users/aklie/data/datasets/tf_perturb_seq/datasets/technology-benchmark_WTC11_TF-Perturb-seq/results

# Echo inputs and number of inputs
echo -e "total inputs: ${#input_files[@]}"
echo -e "input_files: ${input_files[@]}"
echo -e "sample_ids: ${sample_ids[@]}"
echo -e "outdir_path: $outdir_path"

# Step 1 -- merge
echo -e "Running step 1 -- Merge data\n"
cmd="cellcommander merge \
--input_paths ${input_files[@]} \
--outdir_path $outdir_path/merge \
--names ${sample_ids[@]} \
--make-bcs-unique"
echo -e "Running:\n $cmd\n"
#eval $cmd
echo -e "Done with step 1\n"

# Step 2 -- normalize
echo -e "Running step 2 -- Normalize data\n"
cmd="cellcommander normalize \
--input_h5ad_path $outdir_path/merge/merge.h5ad \
--outdir_path $outdir_path/normalize \
--methods log1p sctransform \
--vars_to_regress total_gene_umis percent_mito \
--random-state 1234"
echo -e "Running:\n $cmd\n"
#eval $cmd
echo -e "Done with step 2\n"

# Step 3 -- dim reduction
echo -e "Running step 3 -- Dimensionality reduction\n"
cmd="cellcommander reduce-dimensions \
--input_h5ad_path $outdir_path/normalize/normalize.h5ad \
--outdir_path $outdir_path/reduce_dimensions \
--output_prefix seurat_default_pca \
--method seurat_default \
--obsm_key sctransform_scale_data \
--variable-features-key sctransform_genes \
--random-state 1234"
echo -e "Running:\n $cmd\n"
#eval $cmd
echo -e "Done with step 3\n"

# Step 4 -- correct batch effect 
echo -e "Running step 4 -- Integrate data with Harmony\n"
cmd="cellcommander integrate \
--input_h5ad_path $outdir_path/reduce_dimensions/seurat_default_pca.h5ad \
--outdir_path $outdir_path/integrate \
--vars_to_correct sample \
--obsm_key X_seurat_default \
--method harmonyR \
--corrected_obsm_key X_seurat_default_harmony"
echo -e "Running command:\n $cmd\n"
eval $cmd
echo -e "Done with step 4\n"
date
