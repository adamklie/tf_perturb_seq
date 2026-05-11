#!/bin/bash
#SBATCH --qos=1wk
#SBATCH --nodes=1
#SBATCH --mem=300GB
#SBATCH --job-name=run_igvf_crispr_pipeline
#SBATCH --mail-user=seg95@duke.edu
#SBATCH --mail-type=fail,time_limit


# NOTE: prior to running, install Nextflow using conda-forge
# conda create -n nextflow-env -c conda-forge openjdk=17
conda activate nextflow-env

#curl -s https://get.nextflow.io | bash
#mv nextflow /hpc/group/gersbachlab/seg95/.local/bin

cd /hpc/group/gersbachlab/seg95/CRISPR_Pipeline/
chmod +x bin/*

export GOOGLE_APPLICATION_CREDENTIALS="/hpc/group/gersbachlab/seg95/crispr-pipeline-personal/igvf-pertub-seq-pipeline-85ed0d2862cb.json"
export TOWER_ACCESS_TOKEN="eyJ0aWQiOiAxMTY4Nn0uY2UzZDcyNzE2OWNiYzQ4MDE3NDIxOTI1MTQ0YmZiYzA5OTg2YWE3YQ=="
export SINGULARITY_CACHEDIR=/hpc/group/gersbachlab/seg95/CRISPR_Pipeline/singularity-cache
export MPLCONFIGDIR=/hpc/group/gersbachlab/seg95/CRISPR_Pipeline/mpl-config-dir
#export SINGULARITY_TMPDIR=/hpc/group/gersbachlab/seg95/CRISPR_Pipeline/singularity-tmp

export APPTAINER_CACHEDIR=/hpc/group/gersbachlab/seg95/CRISPR_Pipeline/singularity-cache
export APPTAINER_TMPDIR=/hpc/group/gersbachlab/seg95/CRISPR_Pipeline/singularity-tmp
export TMPDIR=/hpc/group/gersbachlab/seg95/CRISPR_Pipeline/tmp

# Ensure config file has the appropriate filenames and paths
# github.com/pinellolab/CRISPR_Pipeline/tree/main/configs

### Test dataset ###
#nextflow run main.nf -profile google --input example-data/helen_samplesheet.tsv --outdir gs://igvf-pertub-seq-pipeline-data/scratch/seg95/helen_run/
#nextflow run main.nf -profile google --input example-data/example_samplesheet.tsv --outdir gs://igvf-pertub-seq-pipeline-data/scratch/seg95/
#nextflow run main.nf -profile slurm --input example-data/example_samplesheet.tsv --outdir example_output

NXF_HOME=/hpc/group/gersbachlab/seg95/CRISPR_Pipeline/.nextflow
nextflow run main.nf -profile slurm \
	--input /hpc/group/gersbachlab/seg95/crispr-pipeline-personal/example-data/helen_samplesheet_seqspec_poolabcdf.tsv \
	--outdir /hpc/group/gersbachlab/seg95/CRISPR_Pipeline/helen_output_poolabcdf \
	-w /hpc/group/gersbachlab/seg95/CRISPR_Pipeline/work_poolabcdf \
	-c /hpc/group/gersbachlab/seg95/CRISPR_Pipeline/nextflow.config 

# Check status: https://cloud.seqera.io/user/saraecamilli/watch
# Link to Google Cloud Bucket: https://console.cloud.google.com/storage/browser/igvf-pertub-seq-pipeline-data;tab=objects?forceOnBucketsSortingFiltering=true&inv=1&invt=AbynVg&project=igvf-pertub-seq-pipeline&prefix=&forceOnObjectsSortingFiltering=false


