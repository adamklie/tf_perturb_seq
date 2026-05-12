# Samplesheet schema

The samplesheet is a CSV with one row per `(sample, modality)`. The schema is identical across all three steps — only the file-reference cells change.

## Schema (columns)

| Column | Type | Set by | Notes |
|---|---|---|---|
| `sample_name` | string | Step 1 | Stable across steps |
| `modality` | enum | Step 1 | `scRNA`, `gRNA`, or `hash` |
| `R1_path` | file ref | Step 1 → Step 2 | IGVFFI accession → `gs://.../R1.fastq.gz` |
| `R2_path` | file ref | Step 1 → Step 2 | IGVFFI accession → `gs://.../R2.fastq.gz` |
| `seqspec` | file ref | Step 1 → Step 2 → Step 3 | IGVF accession → `gs://.../seqspec.yaml.gz` → `gs://.../patch/seqspec.yaml` |
| `barcode_onlist` | file ref | Step 1 → Step 2 → Step 3 | Same funnel as `seqspec` |
| `guide_design` | file ref | Step 1 → Step 2 → Step 3 | Same funnel; the construct library guide RNA file |
| `barcode_hashtag_map` | file ref | Step 1 → Step 2 → Step 3 | Hashtag→barcode TSV; populated for HTO datasets only |
| `strand_specificity` | enum | Step 1 | Pass-through from measurement set |
| `onlist_method` | enum | Step 1 | Always `"no combination"` (other values unsupported) |
| `analysis_set_accession` | string | Step 1 | `IGVFDS...` for traceability |
| `measurement_set_accession` | string | Step 1 | `IGVFDS...` |
| `auxiliary_set_accession` | string | Step 1 | `IGVFDS...`, blank for scRNA rows |
| `construct_library_set_accession` | string | Step 1 | `IGVFDS...` |

Exact column set may drift; check the actual header on a recent dataset:

```bash
head -1 datasets/Hon_WTC11-benchmark_TF-Perturb-seq/setup/samplesheets/sample_metadata_gcp_*_patched.csv
```

## Funnel by file extension

For each file-reference cell:

```
Step 1 output                Step 2 output                            Step 3 output
─────────────                ─────────────                            ─────────────
IGVFFI<ACC>           →     gs://<bucket>/<ds>/<date>/<ACC>.tsv.gz   →  gs://<bucket>/<ds>/<date>/patch/<ACC>.tsv
                            gs://<bucket>/<ds>/<date>/<ACC>.fastq.gz →  (unchanged — fastqs stay gzipped)
```

R1/R2 fastqs stop at Step 2 (gzipped on GCS). Onlist/guide_design/seqspec/hashtag_map continue through Step 3 to their `patch/` siblings.

## Naming convention

```
sample_metadata.csv                                  # Step 1 output
sample_metadata_gcp_<YYYY_MM_DD>.csv                 # Step 2 output (date auto-tagged from upload date)
sample_metadata_gcp_<YYYY_MM_DD>_patched.csv         # Step 3 output (matches Step 2 date)
```

All three live in `<DATASET_DIR>/setup/samplesheets/`. All three are tracked in git (small, plain CSV) and serve as run-provenance for downstream stages.

## Cross-references

- Portal-side requirements that determine which rows/columns exist: `01-portal-query.md`.
- Local file conventions for un-portal'd inputs: `docs/data/DATA.md` §"Per-dataset layout".
