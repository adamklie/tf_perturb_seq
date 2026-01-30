# Engreitz_WTC11-benchmark_TF-Perturb-seq

## Dataset

WTC11 benchmark TF Perturb-seq dataset from the Engreitz lab.

## Run Pipeline on GCP

```bash
cd $PROJECT_ROOT

# Generate samplesheet
bash datasets/Engreitz_WTC11-benchmark_TF-Perturb-seq/1_generate_per_sample_metadata.sh

# Transfer to GCP -- dry run first
DRY_RUN=true bash datasets/Engreitz_WTC11-benchmark_TF-Perturb-seq/2_upload_to_gcp.sh
bash datasets/Engreitz_WTC11-benchmark_TF-Perturb-seq/2_upload_to_gcp.sh

# Run pipeline on GCP
RUN_IN_BACKGROUND=true bash datasets/Engreitz_WTC11-benchmark_TF-Perturb-seq/3_run_CRISPR_pipeline.sh
```
