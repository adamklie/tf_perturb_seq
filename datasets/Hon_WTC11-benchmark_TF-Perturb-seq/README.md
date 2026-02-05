# Hon_WTC11-benchmark_TF-Perturb-seq

## Dataset

WTC11 benchmark TF Perturb-seq dataset from the Hon lab.

## Run Pipeline on GCP

```bash
cd $PROJECT_ROOT

# Generate samplesheet
bash datasets/Hon_WTC11-benchmark_TF-Perturb-seq/1_generate_per_sample_metadata.sh

# Transfer to GCP -- dry run first
DRY_RUN=true bash datasets/Hon_WTC11-benchmark_TF-Perturb-seq/2_upload_to_gcp.sh
bash datasets/Hon_WTC11-benchmark_TF-Perturb-seq/2_upload_to_gcp.sh

# Run pipeline on GCP
RUN_IN_BACKGROUND=true bash datasets/Hon_WTC11-benchmark_TF-Perturb-seq/3_run_CRISPR_pipeline.sh

# 
bash datasets/Hon_WTC11-benchmark_TF-Perturb-seq/4_run_qc_pipeline.sh
```


Pipeline run reference

distracted_swirles -- my own run with an accompanying tower dashboard at friendly_minsky (Sceptre)
elated_almeida -- my own run with an accompanying tower dashboard at elated_almeida (CLEANSER)
honlab_internal -- result from Hon lab internal run