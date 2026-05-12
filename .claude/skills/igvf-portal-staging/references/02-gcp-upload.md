# Step 2: GCS upload

Wraps `src/tf_perturb_seq/gcp/upload_to_gcp.py`. Transfers IGVF S3 files to GCS via `gcloud transfer jobs` (S3→GCS, no local hop), and uploads local files via `gsutil cp`. Rewrites the samplesheet to use `gs://` paths.

## Prereqs

| Tool | Check |
|---|---|
| `gcloud` CLI | `gcloud --version` |
| `gsutil` | `gsutil --version` |
| Python `requests` | `python3 -c "import requests"` |
| GCP auth | `gcloud auth print-access-token` returns a token |
| Active project | `gcloud config set project igvf-pertub-seq-pipeline` (auto-run by the script) |

**Env vars:**
```bash
export AWS_ROLE_ARN="arn:aws:iam::<ACCOUNT_ID>:role/<ROLE_NAME>"
# IGVF_API_KEY / IGVF_SECRET_KEY optional but recommended (used to resolve accessions → S3 URIs)
```

`AWS_ROLE_ARN` is **required** for the S3→GCS transfer job to assume the right IAM role on the IGVF side. If it's not set, accession transfers will fail silently or with a permission error from the transfer job. Source the value from the team password store / Slack pin.

## Run

Always dry-run first:

```bash
DRY_RUN=true bash <DATASET_DIR>/setup/scripts/2_upload_to_gcp.sh
```

Then actual:

```bash
bash <DATASET_DIR>/setup/scripts/2_upload_to_gcp.sh
```

Common per-dataset edits in the script:

| Variable | Edit when |
|---|---|
| `BASE_DIR` | First time on a new machine |
| `DATASET_NAME` | New dataset |
| `COLUMNS` (commented) | Only transfer a subset of file columns (rare) |

The output filename is auto-dated: `sample_metadata_gcp_$(date +%Y_%m_%d).csv`. Running twice in one day overwrites; running on consecutive days produces two files — keep the latest and use it in Step 3.

## upload_to_gcp.py CLI

```
--input <path.csv>           (required) Input samplesheet (output of Step 1)
--output <path.csv>          (required) Output samplesheet with GCS paths
--gcs-bucket <name>          (required) Bucket without gs:// prefix
--gcs-prefix <path/>         (required) Path within bucket, e.g. "<dataset>/<date>/"
--project <gcp-project>      (required) e.g. igvf-pertub-seq-pipeline
--columns R1_path R2_path... (optional) Restrict to specific file columns
--dry-run                    (optional) Print actions without executing
```

Default columns: `R1_path R2_path seqspec barcode_onlist guide_design barcode_hashtag_map`.

## Output

- GCS layout: `gs://igvf-pertub-seq-pipeline-data/<DATASET_NAME>/<YYYY_MM_DD>/...`
- Samplesheet: `setup/samplesheets/sample_metadata_gcp_<YYYY_MM_DD>.csv` — same schema as Step 1, but file references rewritten to `gs://...` paths.

Verify:
```bash
gsutil ls gs://igvf-pertub-seq-pipeline-data/<DATASET_NAME>/<DATE>/
```

## Nuances

- **S3→GCS transfers are async.** `gcloud transfer jobs create` returns immediately; the job runs server-side. The script waits for each job to finish before continuing. For large fastqs this can take 30+ min per file — leave the script running.
- **gsutil cp for local files.** Local-path entries (rare; most files are IGVF accessions) get a synchronous `gsutil cp`. These can hit your network upload limit.
- **Idempotency.** Existing files in the destination GCS prefix are not re-uploaded. To force re-upload, delete the GCS object first (`gsutil rm`).
- **Filename collisions.** Per the user-level guideline, the uploader uses unique filenames (accession or original basename); collisions across columns should be rare but check the dry-run output if anything looks ambiguous.
- **The samplesheet path keeps the original column structure** — only the file-reference cells change. This is what lets Step 3 walk specific columns.

## Common errors

- **`Not authenticated`** — `gcloud auth login` then re-run.
- **Transfer job stays in `IN_PROGRESS` forever** — check `gcloud transfer jobs list` and `gcloud transfer operations list`. Usually `AWS_ROLE_ARN` is wrong or the S3 object doesn't exist.
- **`AccessDenied` on S3** — IGVF role ARN is stale; refresh from team.
- **`Bucket not found`** — wrong project (`gcloud config get-value project`).

## After Step 2

Spot-check a few GCS paths exist:
```bash
gsutil ls gs://igvf-pertub-seq-pipeline-data/<DATASET_NAME>/<DATE>/ | head
```

Then proceed to Step 3 (`03-patch.md`).
