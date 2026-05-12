# Engreitz WTC11 Benchmark TF Perturb-seq

WTC11 benchmark TF Perturb-seq from the Engreitz lab. IGVF analysis set `IGVFDS5057HJKP`; Synapse mirror [`syn73615563`](https://www.synapse.org/Synapse:syn73615563) (Sid's GCP run, downloaded 2026-02-17).

## Layout

```
setup/                    # Shared input generation
├── scripts/              # 1_–4_*.sh pipeline drivers
├── configs/              # Base Nextflow .config
└── samplesheets/         # Step 1 → 2 → 3 lineage:
                          #   sample_metadata.csv → ..._gcp_<date>.csv → ..._gcp_<date>_patched.csv
<run>/                    # One per Lucas's parameter sweep
├── crispr_pipeline/
│   ├── pipeline_info/    # params_*.json + versions yml (TRACKED — small, run provenance)
│   └── pipeline_outputs/ / pipeline_dashboard/ / anndata/  (all gitignored)
└── calibration/          # FDR-controlled TSVs (analysis tier; result TSVs gitignored)
```

## Pipeline runs

| Local run | Lucas's GCS source | scrublet | method |
|---|---|---|---|
| `cleanser_500_mito_15pc` | `Benchmark_cleanser_500_mito_15pc` | off | cleanser |
| `cleanser_800_mito_15pc` | `Benchmark_cleanser_800_mito_15pc` | off | cleanser |
| `cleanser_extremes_200_mito_15pc` | `Benchmark_cleanser_extremes_200_sceptre_mito_15pc` | off | sceptre |
| `cleanser_extremes_2000_mito_15pc` | `Benchmark_cleanser_extremes_2000_sceptre_mito_15pc` | off | sceptre |
| `cleanser_knee2_mito_15pc` | `Benchmark_cleanser_knee2_mito_15pc` | off | cleanser |
| `scrublet_off_cleanser_800_mito_15pc` | `Benchmark_cleanser_800_mito_15pc` | off | cleanser |
| `scrublet_on_sceptre_800_mito_15pc` | `REAL_SCRUBLETS_..._scrublet_sceptre_800_mito_15pc` | **on** | sceptre |

Lucas's GCS root: `gs://igvf-pertub-seq-pipeline-data/scratch/bioinfolucas/<lucas_dir>/Engreitz_ccPerturb/`.
Folder names can be misleading — verify params from `<run>/crispr_pipeline/pipeline_info/params_*.json`.

## Reproduce

```bash
DS=datasets/Engreitz_WTC11-benchmark_TF-Perturb-seq
bash $DS/setup/scripts/1_generate_per_sample_metadata.sh  # IGVF portal → setup/samplesheets/sample_metadata.csv
bash $DS/setup/scripts/2_upload_to_gcp.sh                 # → setup/samplesheets/sample_metadata_gcp_<date>.csv
bash $DS/setup/scripts/3_patch_gcp_files.sh               # → setup/samplesheets/sample_metadata_gcp_<date>_patched.csv
bash $DS/setup/scripts/4_run_CRISPR_pipeline.sh           # runs Nextflow on GCP
```
