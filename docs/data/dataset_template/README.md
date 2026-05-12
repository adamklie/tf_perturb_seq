# Dataset: {DATASET_NAME}

Replace this header with a one-paragraph description: lab, cell line, condition, IGVF analysis set, Synapse mirror (if any), and links to any active GitHub issues.

## Layout

```
datasets/{DATASET_NAME}/
├── README.md                          # This file
├── setup/                             # Shared input generation
│   ├── scripts/                       # 1_–5_*.sh pipeline drivers
│   ├── configs/                       # Base Nextflow .config(s)
│   └── samplesheets/                  # sample_metadata.csv → _gcp_<date>.csv → _patched.csv
└── <run_name>/                        # One per Nextflow run / parameter set / GCS source
    ├── crispr_pipeline/
    │   ├── pipeline_info/             # params_*.json + versions yml — TRACKED
    │   ├── pipeline_outputs/          # gitignored
    │   ├── pipeline_dashboard/        # gitignored
    │   └── anndata/                   # gitignored
    ├── calibration/                   # FDR-controlled TSVs — gitignored
    ├── qc/                            # mapping_gene / mapping_guide / intended_target — gitignored
    ├── cnmf/<cnmf_run_id>/            # Script/ tracked; Data/, Result/ gitignored
    └── energy_distance/               # configs tracked; image/, logs/, *.csv/.h5mu/.pickle gitignored
```

See [docs/data/DATA.md](../DATA.md) for the conventions (analysis tiers, what's tracked vs gitignored, run-provenance via `pipeline_info/params_*.json`).

## Overview

| Property | Value |
|---|---|
| **Lab** | {LAB_NAME} |
| **Cell Line** | {CELL_LINE} |
| **Differentiation** | {DIFFERENTIATION_STATE} |
| **IGVF Analysis Set** | [{ACCESSION}](https://data.igvf.org/analysis-sets/{ACCESSION}) |
| **Synapse** | {SYNAPSE_ID or N/A} |
| **GCS canonical** | `gs://igvf-pertub-seq-pipeline-data/{DATASET_NAME}/<date>/outs/<run_name>/` |

## Pipeline runs

| Local run | Source | Notes |
|---|---|---|
| `<run_name>` | GCS path / Lucas's sweep label / Synapse ID | scrublet / method / chemistry differences |

Folder names alone are misleading — verify run params from `<run>/crispr_pipeline/pipeline_info/params_*.json`.

## Reproduce

```bash
DS=datasets/{DATASET_NAME}
bash $DS/setup/scripts/1_generate_per_sample_metadata.sh   # IGVF portal → setup/samplesheets/sample_metadata.csv
bash $DS/setup/scripts/2_upload_to_gcp.sh                  # → setup/samplesheets/sample_metadata_gcp_<date>.csv
bash $DS/setup/scripts/3_patch_gcp_files.sh                # → setup/samplesheets/sample_metadata_gcp_<date>_patched.csv
bash $DS/setup/scripts/4_run_CRISPR_pipeline.sh            # runs Nextflow on GCP Batch
bash $DS/setup/scripts/5_run_energy_distance.sh            # production datasets only
```

cNMF runs are kicked off separately under `<run>/cnmf/<cnmf_run_id>/Script/` — see [docs/analysis/cnmf/cNMF.md](../../analysis/cnmf/cNMF.md).

## Notes

{Dataset-specific quirks: seqspec patches, lab-supplied fastq routes, downstream analysis status, open issues.}
