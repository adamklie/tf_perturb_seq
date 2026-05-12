# Step 3: Patch (decompress .tsv.gz on GCS)

Wraps `src/tf_perturb_seq/gcp/patch_gcp_files.py`. The CRISPR Nextflow pipeline requires **uncompressed** TSVs for certain columns. This step streams each `.tsv.gz` through `gunzip` and writes the result back to a `patch/` directory on GCS, then rewrites the samplesheet to point at the uncompressed paths.

## Columns patched (default)

```
seqspec
barcode_onlist
guide_design
barcode_hashtag_map
```

`R1_path` and `R2_path` are intentionally **not** patched — fastqs stay gzipped.

## Prereqs

- Step 2 has completed (`sample_metadata_gcp_<date>.csv` exists with `gs://` paths).
- `gcloud` + `gsutil` authenticated against `igvf-pertub-seq-pipeline`.

## Run

The script's `INPUT_FILE` and `OUTPUT_FILE` are **date-tagged and must be edited each run** (unlike steps 1 and 2 which auto-date). Open `3_patch_gcp_files.sh` and update both to match the Step 2 output filename:

```bash
INPUT_FILE=${DATASET_DIR}/setup/samplesheets/sample_metadata_gcp_<YYYY_MM_DD>.csv
OUTPUT_FILE=${DATASET_DIR}/setup/samplesheets/sample_metadata_gcp_<YYYY_MM_DD>_patched.csv
```

Then dry-run, then actual:

```bash
DRY_RUN=true bash <DATASET_DIR>/setup/scripts/3_patch_gcp_files.sh
bash <DATASET_DIR>/setup/scripts/3_patch_gcp_files.sh
```

## patch_gcp_files.py CLI

```
--input <path.csv>           (required) Step 2 output
--output <path.csv>          (required) Patched output
--columns <col1 col2 ...>    (optional) Override default columns
--patch-dir <gs://...>       (optional) Auto-derived from input paths
--dry-run                    (optional) Print actions, no GCS writes
```

## Patch directory derivation

By default the script derives `<patch-dir>` from the input samplesheet's GCS paths:

```
gs://<bucket>/<dataset>/<date>/...  →  gs://<bucket>/<dataset>/<date>/patch
```

It requires all `.gz` paths to share the same `gs://<bucket>/<dataset>/<date>/` prefix; otherwise it errors out and you must pass `--patch-dir` explicitly. This usually happens if you spliced together samplesheets from two Step 2 runs.

## Mechanics

For each unique `.tsv.gz` path:

```bash
gsutil cat <src.gz> | gunzip | gsutil cp - <patch-dir>/<basename>
```

Streaming — no local disk usage. The samplesheet rewrite replaces every occurrence of `<src.gz>` with `<patch-dir>/<basename>` in the targeted columns.

## Nuances

- **No-op safe.** If a dataset has no `.tsv.gz` references in the patched columns (already uncompressed on the portal), Step 3 still runs and produces `..._patched.csv` — identical to the input but renamed. Stage 2 expects the `_patched` suffix, so always run Step 3.
- **Idempotency.** Re-running checks `gsutil stat` on each destination; existing patched files are not re-decompressed.
- **Streaming cost.** Each decompression is a server-side `gsutil cat | gunzip | gsutil cp`. For typical guide_design / onlist files (few MB each) this is seconds per file.
- **Seqspec patching.** Some datasets need i7/i5 reads stripped from the seqspec YAMLs as part of "patching" (see [docs/data/DATA.md](../../../../docs/data/DATA.md) Step 3 row). That's a manual edit done in the seqspec YAMLs upstream of Step 1; Step 3 itself only does `.gz` decompression.

## Output

`setup/samplesheets/sample_metadata_gcp_<YYYY_MM_DD>_patched.csv` — same schema as Step 2 output but with patched columns rewritten to `gs://.../patch/<file>.tsv`.

This is the file `4_run_CRISPR_pipeline.sh` (Stage 2) consumes.

## Common errors

- **`Cannot auto-derive patch directory`** — input has paths spanning multiple `<dataset>/<date>/` prefixes. Pass `--patch-dir gs://.../patch` explicitly.
- **`gsutil cat ... CommandException: No URLs matched`** — the `.tsv.gz` source path is wrong in the samplesheet (probably edited by hand). Re-run Step 2 or fix the path.
- **Hanging on large file** — the streaming pipe can stall; the 300s timeout in `decompress_gcs_file` will surface an error. For very large files, raise the timeout or use a non-streaming approach.

## After Step 3

```bash
head -1 setup/samplesheets/sample_metadata_gcp_<DATE>_patched.csv  # verify schema unchanged
gsutil ls gs://igvf-pertub-seq-pipeline-data/<DATASET_NAME>/<DATE>/patch/
```

Hand off to Stage 2 by editing `4_run_CRISPR_pipeline.sh` to point at the patched samplesheet.
