# Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2

## Dataset

WTC11 benchmark TF Perturb-seq dataset from the Gersbach lab (HT v2 chemistry).

## Run Pipeline on GCP

```bash
cd $PROJECT_ROOT

# Generate samplesheet
bash datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/1_generate_per_sample_metadata.sh

# Transfer to GCP -- dry run first
DRY_RUN=true bash datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2_upload_to_gcp.sh
bash datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2_upload_to_gcp.sh

# Run pipeline on GCP
RUN_IN_BACKGROUND=true bash datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/3_run_CRISPR_pipeline.sh

# Run QC pipeline
bash datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/4_run_qc_pipeline.sh
```

Pipeline run reference

red_lettuce -- from Synapse (https://www.synapse.org/Synapse:syn73579418) run by Alex Barrera with Sceptre guide inference
