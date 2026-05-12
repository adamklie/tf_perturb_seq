---
name: igvf-portal-staging
description: Stage 1 of the TFP3 pipeline. Query the IGVF portal for an analysis set, transfer fastqs and metadata files from IGVF S3 to GCS, and produce the patched samplesheet that feeds the CRISPR Nextflow pipeline. Triggers on keywords like Stage 1, IGVF portal, analysis set, accession, sample_metadata, sample metadata, samplesheet, GCS upload, GCP upload, S3 to GCS, patch GCS files, decompress tsv.gz, barcode_onlist, guide_design, IGVFDS, onboard dataset, new dataset.
user_invocable: true
---

# Stage 1: IGVF Portal → GCS Staging

You are an interactive assistant for running Stage 1 of the TFP3 pipeline. Stage 1 turns an IGVF analysis set accession into a GCS-resident, patched samplesheet that Stage 2 (CRISPR Nextflow) can consume.

## Pipeline position

```
[1] IGVF Portal + GCS staging  ← you are here
   ↓ produces setup/samplesheets/sample_metadata_gcp_<date>_patched.csv
[2] CRISPR Nextflow on GCP Batch
[3] QC  →  [4] Energy distance  →  [5] cNMF
```

The output of Stage 1 is consumed by `4_run_CRISPR_pipeline.sh` in the same `setup/scripts/` directory.

## Constants

```
REPO_ROOT (local Mac):     /Users/adamklie/Desktop/tfp3/tf_perturb_seq
REPO_ROOT (legacy/HPC):    /Users/adamklie/Desktop/projects/tf_perturb_seq
GCP_PROJECT:               igvf-pertub-seq-pipeline
GCS_BUCKET:                igvf-pertub-seq-pipeline-data
DATASET_DIR:               <REPO_ROOT>/datasets/<DATASET_NAME>
SAMPLESHEETS_DIR:          <DATASET_DIR>/setup/samplesheets/
SCRIPTS_DIR:               <DATASET_DIR>/setup/scripts/
PORTAL_API:                https://api.data.igvf.org
PORTAL_SEARCH:             https://data.igvf.org/search/?type=MeasurementSet&preferred_assay_titles=Perturb-seq&collections=TF+Perturb-seq+Project
```

The existing scripts hard-code `BASE_DIR=/Users/adamklie/Desktop/projects/tf_perturb_seq`. When running on the current machine, the only edit required is `BASE_DIR` at the top of each `setup/scripts/N_*.sh`. Do not change the body of the scripts.

## Samplesheet funnel

Each step produces a new samplesheet in `setup/samplesheets/`; the next step reads it.

| Step | Script | Input | Output |
|---|---|---|---|
| 1 | `1_generate_per_sample_metadata.sh` | IGVF accession (e.g. `IGVFDS4761PYUO`) | `sample_metadata.csv` |
| 2 | `2_upload_to_gcp.sh` | `sample_metadata.csv` | `sample_metadata_gcp_<YYYY_MM_DD>.csv` |
| 3 | `3_patch_gcp_files.sh` | `sample_metadata_gcp_<date>.csv` | `sample_metadata_gcp_<date>_patched.csv` |

For the per-step column changes, read `references/samplesheet-schema.md`.

## Step 0: Identify the substep and target dataset

Ask (or infer from context):

1. **Which substep?** 1 (portal query), 2 (GCS upload), 3 (patch), or all three for a new dataset.
2. **Which dataset?** Path under `datasets/` (e.g., `datasets/Hon_WTC11-benchmark_TF-Perturb-seq/`) or a new dataset name to scaffold.
3. **Is this a brand-new dataset?** If yes, run the scaffolder first (see "Scaffolding a new dataset" below).

Then read the matching reference file before collecting parameters.

| Substep | Reference file |
|---|---|
| 1: portal query | `references/01-portal-query.md` |
| 2: GCS upload | `references/02-gcp-upload.md` |
| 3: patch | `references/03-patch.md` |
| Samplesheet schema | `references/samplesheet-schema.md` |

## Step 1: Portal query

Wraps `src/tf_perturb_seq/portal/generate_per_sample.py`. Pulls the analysis set object from the IGVF portal, walks measurement sets / auxiliary sets / construct library set, and writes one row per (sample, modality).

**Prereqs:** `IGVF_API_KEY` and `IGVF_SECRET_KEY` env vars (or `--keypair` JSON).

**Run:**
```bash
bash <DATASET_DIR>/setup/scripts/1_generate_per_sample_metadata.sh
```

Common edits in the script: `BASE_DIR`, `DATASET_NAME`, `ACCESSION`, and (only if the portal lacks linked seqspecs) `--hash_seqspec / --rna_seqspec / --sgrna_seqspec` flags pointing at local YAMLs under `bin/1_CRISPR_pipeline/seqspec/yaml_files/`.

For the full portal prereq checklist (what the analysis set must contain) and troubleshooting, see `references/01-portal-query.md`.

## Step 2: GCS upload

Wraps `src/tf_perturb_seq/gcp/upload_to_gcp.py`. Transfers IGVF S3 files (by accession) to GCS via gcloud transfer jobs; uploads local files via `gsutil cp`. Rewrites the samplesheet with `gs://...` paths.

**Prereqs:** `gcloud` + `gsutil` installed; `gcloud auth login` done; `AWS_ROLE_ARN` env var set; `requests` Python library.

**Run a dry-run first** to verify which transfers will fire and the GCS destination:
```bash
DRY_RUN=true bash <DATASET_DIR>/setup/scripts/2_upload_to_gcp.sh
```

Then actual run:
```bash
bash <DATASET_DIR>/setup/scripts/2_upload_to_gcp.sh
```

GCS layout: `gs://igvf-pertub-seq-pipeline-data/<DATASET_NAME>/<YYYY_MM_DD>/`.

For auth troubleshooting and what to do if a transfer job hangs, see `references/02-gcp-upload.md`.

## Step 3: Patch (decompress .tsv.gz)

Wraps `src/tf_perturb_seq/gcp/patch_gcp_files.py`. The CRISPR pipeline needs **uncompressed** TSVs for `seqspec`, `barcode_onlist`, `guide_design`, `barcode_hashtag_map`. This step streams each `.tsv.gz` through `gunzip` and writes the uncompressed file under `<gcs_root>/patch/`, then rewrites the samplesheet to point at the patched paths.

**Edit before running:** in `3_patch_gcp_files.sh`, set `INPUT_FILE` and `OUTPUT_FILE` to the dated samplesheet from Step 2.

```bash
DRY_RUN=true bash <DATASET_DIR>/setup/scripts/3_patch_gcp_files.sh   # preview
bash <DATASET_DIR>/setup/scripts/3_patch_gcp_files.sh                # actual
```

If a dataset has no `.tsv.gz` references (already uncompressed on the portal), Step 3 is still safe to run — it will be a no-op aside from rewriting the samplesheet.

For details on which columns get patched and when to override `--patch-dir`, see `references/03-patch.md`.

## Scaffolding a new dataset

For a brand-new dataset, generate the three setup scripts (and the `setup/{configs,samplesheets}` skeleton) before running Step 1:

```bash
python3 .claude/skills/igvf-portal-staging/scripts/scaffold_setup_scripts.py \
  --dataset-name <Lab>_<CellLine>-<cond>_TF-Perturb-seq \
  --accession IGVFDS<XXXXXXXX>
```

The scaffolder is idempotent — it will not overwrite an existing `setup/scripts/N_*.sh` unless `--force` is passed.

After scaffolding:
1. Review the three generated scripts; edit `BASE_DIR` if needed.
2. If the portal accession lacks linked seqspecs, drop YAML fallbacks in `bin/1_CRISPR_pipeline/seqspec/yaml_files/` and uncomment the `--*_seqspec` lines in `1_generate_per_sample_metadata.sh`.
3. Add the dataset's `README.md` from `docs/data/dataset_template/README.md`.

## Handoff to Stage 2

When Step 3 completes, the patched samplesheet path is:

```
<DATASET_DIR>/setup/samplesheets/sample_metadata_gcp_<YYYY_MM_DD>_patched.csv
```

Edit `4_run_CRISPR_pipeline.sh` to point its `--input` (or equivalent) at this file, then hand off to Stage 2.

## Important notes

- `BASE_DIR` is the **only** per-environment edit. Bodies of `N_*.sh` should be identical across datasets.
- All three scripts support `DRY_RUN=true` for steps 2 and 3. Step 1 has no dry-run; it only writes a CSV.
- Run-provenance for downstream stages lives in `<run>/crispr_pipeline/pipeline_info/params_*.json`. Stage 1's samplesheets are tracked under `setup/samplesheets/` so you can replay later.
- Portal prereqs are strict (analysis set must group measurement + auxiliary + construct library sets correctly, every measurement set needs `strand_specificity`, etc.). If Step 1 errors out, fix the portal first — do not work around it.
- For the canonical portal-side reference, see [docs/data/DATA.md](../../../docs/data/DATA.md) §"IGVF portal structure".
