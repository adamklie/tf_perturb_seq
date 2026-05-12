# Step 1: Portal query

Wraps `src/tf_perturb_seq/portal/generate_per_sample.py`. Pulls the analysis set object via the IGVF portal API and writes one row per `(sample, modality)` to `setup/samplesheets/sample_metadata.csv`.

## Prereqs (portal side)

The analysis set must already be in a valid state. The Python script asserts most of these — fix on the portal, not in code.

- **Analysis Set**
  - `input_file_sets` lists every measurement set and auxiliary set to process together.
  - `construct_library_sets` has **exactly one** entry (calculated from member sets).
- **Measurement Set(s)** (one per 10x lane/pool)
  - `strand_specificity` populated.
  - `onlist_files` populated.
  - `onlist_method` is `"no combination"`. Other values are not supported by the pipeline.
- **Auxiliary Set(s)** (gRNA, HTO)
  - `measurement_sets` link back to the parent measurement set.
  - If `file_set_type == "cell hashing barcode sequencing"`, `barcode_map` (hashtag→barcode TSV) is required.
- **Construct Library Set**
  - Exactly one `integrated_content_files` with `content_type: "guide RNA sequences"` and `status ∈ {"in progress", "preview", "released"}`. Deprecated guide files must be `revoked` / `archived` / `deleted` / `replaced`.
- **Sequence Files**
  - `content_type: "reads"`, `illumina_read_type ∈ {R1, R2}`, `status ∉ {"deleted", "revoked"}`.
- **Seqspecs**
  - One per modality (`scRNA seq`, `gRNA seq`, optional `cell hashing barcode seq`).
  - Linked to sequence files via `seqspec_of`, `status ∈ {"in progress", "preview", "released"}`, and `upload_status == "validated"`.
  - All sequence files in a modality share the same seqspec read index.
  - If a portal seqspec is missing, pass a fallback YAML via CLI flag (see below).

## Prereqs (local)

```bash
export IGVF_API_KEY=...
export IGVF_SECRET_KEY=...
# or pass --keypair path/to/keypair.json with {"key": "...", "secret": "..."}
```

## Run

```bash
bash <DATASET_DIR>/setup/scripts/1_generate_per_sample_metadata.sh
```

Common per-dataset edits in the script:

| Variable | Edit when |
|---|---|
| `BASE_DIR` | First time on a new machine |
| `DATASET_NAME` | New dataset |
| `ACCESSION` | New dataset (e.g. `IGVFDS4761PYUO`) |
| `--*_seqspec` lines | Portal seqspecs missing for some modalities (point at local YAML under `bin/1_CRISPR_pipeline/seqspec/yaml_files/`) |

## generate_per_sample.py CLI

```
--accession <IGVFDS...>     (required) Analysis Set accession
--output <path.csv>         (required) Output CSV path
--keypair <path.json>       (optional) Keypair JSON {key, secret}
--hash_seqspec <path.yml>   (optional) Fallback for cell hashing barcode modality
--rna_seqspec <path.yml>    (optional) Fallback for scRNA modality
--sgrna_seqspec <path.yml>  (optional) Fallback for gRNA modality
```

## Output

`setup/samplesheets/sample_metadata.csv` — one row per `(sample, modality)`. See `samplesheet-schema.md` for the full column list. File references at this stage are **IGVF accessions** (e.g. `IGVFFI...`) — they get rewritten to `gs://...` in Step 2.

## Common errors

- **`Datasets with multiple guide libraries are not currently supported`** — analysis set has >1 construct library set, or the construct library set has >1 active guide RNA file. Resolve on portal.
- **`Guide libraries are required for running this pipeline`** — no active guide RNA file. Check `status` of the `integrated_content_files`.
- **Missing seqspec for modality X** — pass `--*_seqspec` fallback. If the YAML doesn't exist locally yet, adapt one from a sibling dataset under `datasets/<other>/bin/1_CRISPR_pipeline/seqspec/yaml_files/`.
- **HTTP 401 from portal API** — `IGVF_API_KEY`/`IGVF_SECRET_KEY` not set, or expired. The portal also requires the key to be associated with an account that can view the analysis set.

## After Step 1

Eyeball the output CSV:
- One row per (sample × modality). Most benchmark datasets have 2–4 rows; production datasets have many more.
- `guide_design`, `barcode_onlist`, `barcode_hashtag_map` are populated.
- File reference cells are IGVF accessions, not paths.

Then proceed to Step 2 (`02-gcp-upload.md`).
