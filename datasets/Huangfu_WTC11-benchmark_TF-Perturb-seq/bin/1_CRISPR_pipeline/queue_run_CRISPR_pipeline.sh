#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=20G


nextflow run main.nf -profile google -with-tower --input 2025_11_26_huangfu_cloud_sample_metadata.csv --outdir gs://igvf-pertub-seq-pipeline-data/beer_test/huangfu_WTC11_benchmark_output 

