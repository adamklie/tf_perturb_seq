# Engreitz WTC11 Benchmark TF Perturb-seq

WTC11 benchmark TF Perturb-seq dataset from the Engreitz lab.
Data sourced from IGVF portal (analysis set: IGVFDS5057HJKP) and Synapse (`syn73615563`, Sid's GCP run, downloaded 2026-02-17).

## Layout

```
Engreitz_WTC11-benchmark_TF-Perturb-seq/
├── setup/                              # Shared input generation (Adam's code)
│   ├── scripts/                        # 1_–4_*.sh pipeline drivers
│   ├── configs/                        # Base Nextflow .config
│   └── samplesheets/                   # Canonical samplesheet
│       ├── sample_metadata.csv         # Final (GCS-path patched) — canonical
│       └── meta/                       # Provenance trail: raw IGVF → uploaded → patched
└── <variant>/                          # One per Lucas's Nextflow parameter sweep
    ├── crispr_pipeline/
    │   ├── pipeline_outputs/           # Bulk Nextflow output (gitignored)
    │   ├── pipeline_dashboard/         # Dashboards (gitignored)
    │   └── anndata/                    # h5ad/h5mu (gitignored, where present)
    └── calibration/                    # FDR-controlled cis/trans/pathway/direct_target TSVs
                                        # (analysis tier; result TSVs gitignored)
```

Current variants on disk (Lucas's GCP parameter sweep):
`cleanser_500_mito_15pc`, `cleanser_800_mito_15pc`, `cleanser_extremes_200_mito_15pc`,
`cleanser_extremes_2000_mito_15pc`, `cleanser_knee2_mito_15pc`,
`scrublet_off_cleanser_800_mito_15pc`, `scrublet_on_sceptre_800_mito_15pc`.

## Canonical sources

| Resource | Location |
|---|---|
| Pipeline run (GCS) | `gs://igvf-pertub-seq-pipeline-data/Engreitz_WTC11-benchmark_TF-Perturb-seq/` |
| Pipeline mirror (Synapse) | [`syn73615563`](https://www.synapse.org/Synapse:syn73615563) — Sid's GCP run |
| Base Nextflow config (GCS) | `gs://.../Engreitz_WTC11-benchmark_TF-Perturb-seq.config` |
| Canonical samplesheet (GCS) | `gs://.../sample_metadata_gcp_2026_02_26_patched.csv` |

`setup/configs/Engreitz_WTC11-benchmark_TF-Perturb-seq_2026_03_11.config` is the local copy of the GCS base config.
Lucas's per-variant configs (the parameter sweeps) aren't currently mirrored to GCS — they'd live in his local Nextflow runs.

## Running the pipeline

```bash
cd $PROJECT_ROOT
DS=datasets/Engreitz_WTC11-benchmark_TF-Perturb-seq

# 1. Generate samplesheet
bash $DS/setup/scripts/1_generate_per_sample_metadata.sh

# 2. Transfer fastqs to GCP (dry run first)
DRY_RUN=true bash $DS/setup/scripts/2_upload_to_gcp.sh
bash $DS/setup/scripts/2_upload_to_gcp.sh

# 3. Patch GCP files (decompress seqspecs/onlists/guide-designs)
DRY_RUN=true bash $DS/setup/scripts/3_patch_gcp_files.sh
bash $DS/setup/scripts/3_patch_gcp_files.sh

# 4. Validate GCS paths
uv run python src/tf_perturb_seq/gcp/validate_gcp_paths.py \
    --input $DS/setup/samplesheets/sample_metadata.csv

# 5. Run CRISPR pipeline (varies per variant — Lucas's parameter sweeps)
RUN_IN_BACKGROUND=true bash $DS/setup/scripts/4_run_CRISPR_pipeline.sh
```

## Pipeline run reference

**cobalt_heron** -- downloaded from Synapse (`syn73615563`) on 2026-02-17
