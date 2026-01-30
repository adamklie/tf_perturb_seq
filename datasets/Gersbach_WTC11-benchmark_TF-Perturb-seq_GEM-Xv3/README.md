# Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3

## Dataset

WTC11 benchmark TF Perturb-seq dataset from the Gersbach lab (GEM-X v3 chemistry).

## Run Pipeline on GCP

```bash
cd $PROJECT_ROOT

# Generate samplesheet
bash datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/1_generate_per_sample_metadata.sh

# Transfer to GCP -- dry run first
DRY_RUN=true bash datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2_upload_to_gcp.sh
bash datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2_upload_to_gcp.sh

# Run pipeline on GCP
RUN_IN_BACKGROUND=true bash datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/3_run_CRISPR_pipeline.sh
```
