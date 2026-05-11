# Hon WTC11 Benchmark TF Perturb-seq

WTC11 benchmark TF Perturb-seq dataset from the Hon lab. Data sourced from IGVF portal (analysis set: IGVFDS4761PYUO).

## Files

| File | Description |
|------|-------------|
| `sample_metadata_gcp_2026_02_26.csv` | Latest samplesheet with updated seqspecs. This is the input to the CRISPR pipeline. |
| `Hon_WTC11-benchmark_TF-Perturb-seq_2026_03_11.config` | Nextflow config for the CRISPR pipeline (supports local, SLURM, and GCP profiles). |
| `5_pipeline_comparison.ipynb` | Notebook comparing pipeline runs. |
| `_samplesheets/` | Intermediate samplesheets from earlier pipeline steps. |
| `_configs/` | Old pipeline configs (`beautiful_tornado`, `harmless_teacher`, `introverted_frankenstein`). |
| `_seqspecs/` | Old seqspec yamls (superseded by seqspecs in latest samplesheet). |
| `_cNMF/` | Hon lab internal cNMF results (73 GB). |
| `runs/` | Current pipeline dashboard outputs (`cleanser_unified`, `elated_almeida`, `sceptre_v11`). |

## Scripts (run in order)

| Script | Description |
|--------|-------------|
| `1_generate_per_sample_metadata.sh` | Pulls per-sample metadata from the IGVF portal. |
| `2_upload_to_gcp.sh` | Transfers IGVF files to GCS. Supports `DRY_RUN=true`. |
| `3_run_CRISPR_pipeline.sh` | Runs the CRISPR Nextflow pipeline on GCP. Supports `RUN_IN_BACKGROUND=true`. |
| `4_run_qc_pipeline.sh` | Runs QC on pipeline output (`inference_mudata.h5mu`). Supports `--dry-run`. |
| `5_pipeline_comparison.ipynb` | Compares pipeline run results. |
| `6_run_energy_distance.sh` | Runs energy distance pipeline. Supports `--dry-run`, `--use-harmony`, `--force`. |
| `7_run_cnmf.sh` | Runs cNMF gene program discovery. Supports `--dry-run`, `--slurm`, `--phase`. |

## Running the pipeline

```bash
cd $PROJECT_ROOT

# 1. Generate samplesheet
bash datasets/Hon_WTC11-benchmark_TF-Perturb-seq/1_generate_per_sample_metadata.sh

# 2. Transfer to GCP (dry run first)
DRY_RUN=true bash datasets/Hon_WTC11-benchmark_TF-Perturb-seq/2_upload_to_gcp.sh
bash datasets/Hon_WTC11-benchmark_TF-Perturb-seq/2_upload_to_gcp.sh

# 3. Run CRISPR pipeline
RUN_IN_BACKGROUND=true bash datasets/Hon_WTC11-benchmark_TF-Perturb-seq/3_run_CRISPR_pipeline.sh

# 4. Run QC
bash datasets/Hon_WTC11-benchmark_TF-Perturb-seq/4_run_qc_pipeline.sh

# 6. Run energy distance
bash datasets/Hon_WTC11-benchmark_TF-Perturb-seq/6_run_energy_distance.sh

# 7. Run cNMF
bash datasets/Hon_WTC11-benchmark_TF-Perturb-seq/7_run_cnmf.sh
```

## Samplesheet history

| Version | File | Notes |
|---------|------|-------|
| Original | `sample_metadata_gcp_2026_01_18.csv` | Generic seqspecs (rna/guide/hash_seqspec.yml) |
| Patched | `sample_metadata_gcp_2026_01_18_patched.csv` | Decompressed barcode/guide files to `patch/`, cleared lane column |
| Updated seqspecs | `sample_metadata_gcp_2026_02_15.csv` | Replaced seqspecs with `gary_IGVFDS4761PYUO_{rna,guide,hash}.yaml` |
| Latest | `sample_metadata_gcp_2026_02_26.csv` | Most recent samplesheet |

## Pipeline run reference

- **distracted_swirles** -- own run with tower dashboard at friendly_minsky (Sceptre)
- **elated_almeida** -- own run with tower dashboard at elated_almeida (CLEANSER)
- **honlab_internal** -- result from Hon lab internal run
- **golden_falcon** -- downloaded from Synapse (syn73743227) on 2026-02-17
