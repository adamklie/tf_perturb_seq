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


## Sample Metadata

| Version | File | Notes |
|---------|------|-------|
| Original | `sample_metadata_gcp_2026_01_18.csv` | Generic seqspecs (rna/guide/hash_seqspec.yml) |
| Patched | `sample_metadata_gcp_2026_01_18_patched.csv` | Decompressed barcode/guide files to `patch/`, cleared lane column |
| New seqspecs | `sample_metadata_gcp_2026_02_15.csv` | Replaced seqspec column with new shared seqspecs from `scratch/bioinfolucas/new_03_seqspec/` |

### 2026-02-15: seqspec update
- Based on `sample_metadata_gcp_2026_01_18_patched.csv`
- Replaced all seqspec paths:
  - scRNA: `gary_IGVFDS4761PYUO_rna.yaml`
  - gRNA: `gary_IGVFDS4761PYUO_guide.yaml`
  - hash: `gary_IGVFDS4761PYUO_hash.yaml`
- All other columns (barcode_onlist, guide_design, etc.) preserved from patched version
- Uploaded to `gs://igvf-pertub-seq-pipeline-data/Hon_WTC11-benchmark_TF-Perturb-seq/`

## Pipeline run reference

distracted_swirles -- my own run with an accompanying tower dashboard at friendly_minsky (Sceptre)
elated_almeida -- my own run with an accompanying tower dashboard at elated_almeida (CLEANSER)
honlab_internal -- result from Hon lab internal run
golden_falcon -- downloaded from Synapse (syn73743227) on 2026-02-17