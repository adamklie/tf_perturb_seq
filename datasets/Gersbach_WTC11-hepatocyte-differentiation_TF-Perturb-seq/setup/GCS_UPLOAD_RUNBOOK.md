# Gersbach Hep — GCS upload runbook (Stage 1, Step 2 + 3)

Date started: 2026-06-10. Tracking issue: [#28](https://github.com/adamklie/tf_perturb_seq/issues/28).
Companion: [PORTAL_STATUS_2026_06_10.md](PORTAL_STATUS_2026_06_10.md) (portal state), seqspec validation under `scratch/hep_seqspec_validation/`.

This run was stress-tested before execution (`/grill-me`). The notes below capture the
plan, the failure modes we guarded against, and how to redo a failed/partial run.

## Inputs
- **Samplesheet (Step 2 input):** `setup/samplesheets/sample_metadata_localseqspec_2026_06_10.csv`
  - 1,503 rows = 751 scRNA + 752 gRNA pairs (the one S10 run1/lane3 scRNA R1 is absent on the
    portal and is correctly dropped — see PORTAL_STATUS).
  - `R1_path`/`R2_path`/`barcode_onlist` (`IGVFFI1697XAXI`)/`guide_design` (`IGVFFI8270UPKB`) are IGVF accessions.
  - `seqspec` is the **absolute local path** to the corrected + validated YAMLs
    (`setup/seqspec/{rna,guide}_seqspec.yml`; UMI=12, R1=28, guide construct matches the reads —
    validated with crispr_validator, all regions perfect/close_enough).
- **GCS destination:** `gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/2026_06_10/`
- **Auth:** gcloud active account = `adamklie@igvf-pertub-seq-pipeline.iam.gserviceaccount.com`
  (service account; confirmed it can list/create transfer jobs — prior `igvf-upload-*` jobs succeeded).
  Env: `AWS_ROLE_ARN=arn:aws:iam::…:role/S3toPertubSeqGoogleCloudTransfer`, `IGVF_API_KEY`, `IGVF_SECRET_KEY`.

## Why not just run `2_upload_to_gcp.sh`
The driver hardcodes `BASE_DIR=/carter/...` (NRNB) and `INPUT_FILE=sample_metadata.csv` (the old
mixed-seqspec sheet). We run from the Mac against the validated local-seqspec sheet, so we invoke
`src/tf_perturb_seq/gcp/upload_to_gcp.py` directly with explicit flags instead of editing the driver.

## Known failure modes (guarded)
1. **`--include-prefixes` >1000 per transfer job.** `upload_to_gcp.py` builds ONE
   `gcloud transfer jobs create` per source S3 bucket with one include-prefix per object
   (~3,008 for Hep). Google STS caps includePrefixes at **1,000** → the monolithic job may be
   rejected OR silently truncated. → **Dry-run first to get per-bucket counts; if any bucket >1000,
   chunk into ≤1000-prefix jobs (or canary one sub-pool).** (Counts recorded below.)
2. **Async transfer.** Step 2 fires jobs and returns; transfers run server-side. Do NOT run Step 3
   off the built-in spot-check (it only checks the first file). → Gate on completion + completeness (below).
3. **Step 3 config drift.** `3_patch_gcp_files.sh` defaults to `DATE=2026_06_03` and reads
   `sample_metadata_gcp_${DATE}.csv`. Set `DATE=2026_06_10` and `BASE_DIR` (Mac) before running it.

## Procedure

### Step 2 — S3→GCS transfer (read-only dry-run first)
```bash
DS=datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq
export AWS_ROLE_ARN IGVF_API_KEY IGVF_SECRET_KEY
uv run --no-project --with requests python -u src/tf_perturb_seq/gcp/upload_to_gcp.py \
  --input  $DS/setup/samplesheets/sample_metadata_localseqspec_2026_06_10.csv \
  --output $DS/setup/samplesheets/sample_metadata_gcp_2026_06_10.csv \
  --gcs-bucket igvf-pertub-seq-pipeline-data \
  --gcs-prefix Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/2026_06_10/ \
  --project igvf-pertub-seq-pipeline --no-confirm --dry-run     # drop --dry-run for real run
```
Dry-run = read-only (resolves s3_uris, prints the gcloud command; creates no transfers). Inspect
per-bucket prefix counts; chunk if >1000.

### Verify completion + completeness (the gate to Step 3)
```bash
# 1. all transfer jobs SUCCESS
gcloud transfer operations list --project=igvf-pertub-seq-pipeline
# 2. expected-vs-present over the WHOLE sheet (lists any missing gs:// object)
uv run --no-project python src/tf_perturb_seq/gcp/validate_gcp_paths.py \
  --input $DS/setup/samplesheets/sample_metadata_gcp_2026_06_10.csv
# Expect ~3,010 objects present (3,006 fastqs + onlist + guide + 2 seqspecs). Re-transfer any missing.
```

### Step 3 — patch (decompress tsv.gz/csv.gz, upload local seqspecs, rewrite sheet)
Edit `3_patch_gcp_files.sh`: `BASE_DIR` (Mac), `DATE=2026_06_10`. Then:
```bash
DRY_RUN=true bash $DS/setup/scripts/3_patch_gcp_files.sh   # preview
bash $DS/setup/scripts/3_patch_gcp_files.sh                # writes …/2026_06_10/patch/ + *_patched.csv
```
Step 3 already handles our all-local-seqspec sheet (empty `.yaml.gz` list; rewrites
`rna_/guide_seqspec.yml` → `patch/…`; uploads the 2 local YAMLs) and the `guide_design` `.csv.gz`
(`.csv→.tsv`). Output: `sample_metadata_gcp_2026_06_10_patched.csv` → hand to Stage 2.

## Redo / resume
- Transfers use `--overwrite-when=never` and a dated prefix → re-running is **resumable** (already-present
  objects are skipped). `validate_gcp_paths.py` names exactly which objects are missing; re-transfer just those.
- A botched run can be discarded by deleting the dated prefix:
  `gsutil -m rm -r gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/2026_06_10/`
  (only this dataset's dated folder; verify the path before running).

## Run log
- 2026-06-10: dry-run #1 (`scratch/hep_gcs_upload/dryrun.log`). Resolved 3,008 IGVF files (0 failed) + 2 local seqspecs.
  Per-bucket: **igvf-private 2,997**, igvf-files 9, igvf-public 2. → The igvf-private job would carry
  **2,997 include-prefixes in one command** — over the STS 1,000 limit (and a ~210 KB arg line). **Confirmed the
  chunking bug.**
- 2026-06-10: **Fixed `src/tf_perturb_seq/gcp/upload_to_gcp.py`** — `upload_igvf_files` now chunks include-prefixes
  into ≤1,000-per-job batches (`MAX_INCLUDE_PREFIXES = 1000`), one transfer job per chunk. General fix (benefits any
  >1,000-file dataset). igvf-private → 3 jobs (1000/1000/997); igvf-files, igvf-public → 1 each.
- 2026-06-10: dry-run #2 (`scratch/hep_gcs_upload/dryrun2.log`) — chunking verified: igvf-private → 3 jobs
  (1000/1000/997), igvf-public → 1 (2), igvf-files → 1 (9), + 2 local seqspecs. All jobs ≤1000 prefixes. 0 failed.
- 2026-06-10: **REAL transfer fired** (`scratch/hep_gcs_upload/realrun.log`; 5 jobs, 3,010 refs, 0 create-failures).
  Job outcomes: igvf-public ✅, igvf-private part3 ✅, igvf-private part1/2 IN_PROGRESS, **igvf-files job FAILED
  (PERMISSION_DENIED)**. 2 local seqspec YAMLs uploaded OK.
- **igvf-files permission gap**: 9 fastq reads (status `in progress`, total 21.7 GB) live in `s3://igvf-files/`
  (IGVF upload-staging bucket). The `S3toPertubSeqGoogleCloudTransfer` role can read igvf-private/igvf-public
  but NOT igvf-files → those 9 can't go via STS. They are single mates of 9 lanes (one per sub-pool incl. S10
  gRNA R1). Accessions: IGVFFI1531RPZT, IGVFFI1734TEID, IGVFFI4660TYLY, IGVFFI5688DYDJ, IGVFFI6646SQPQ,
  IGVFFI6683VISA, IGVFFI7202QJLS, IGVFFI7347KQYP, IGVFFI9998TPWM. Expected GCS dest paths recorded in the
  output sheet. WORKAROUND options: (a) portal `@@download` + `gsutil cp` to the expected gs:// paths
  (we have IGVF auth); (b) wait for IGVF to release them to igvf-private, then re-run STS (resumable). _decision TBD_
- 2026-06-10: all 4 STS jobs (igvf-private ×3, igvf-public ×1 = 2,999 fastqs) reached **SUCCESS**; +2 seqspec YAMLs = 3,001 in GCS.
- 2026-06-10: **9 igvf-files reads recovered via NRNB** (`aklie@nrnb-login`, SA `adamklie@…iam`, scratch `/carter/users/aklie/scratch/`).
  Portal `@@download` → verify → `gsutil cp`. 5 matched portal md5 (`recover9.sh`). **4 were `upload_status=invalidated`** with
  stale (byte-rotated) md5 — verified instead by size + gzip + cross-download stability (one checked on Mac: valid 28bp fastq,
  size-exact) and recovered by `recover3.sh` (size+gzip check): IGVFFI1734TEID, IGVFFI6646SQPQ (uploaded from Mac),
  IGVFFI6683VISA, IGVFFI7347KQYP.
- 2026-06-10: **completeness gate PASSED — 3,010/3,010 objects present, 0 missing** (fast check: `gsutil ls -r` vs expected
  gs:// paths from the output sheet; match `\.(gz|yml)$`). Output sheet: `sample_metadata_gcp_2026_06_10.csv`.
  NOTE: `validate_gcp_paths.py` works but is slow (~3,010 sequential `gsutil stat`); the listing-diff is the fast equivalent.

## Caveats carried forward
- **4 recovered reads are `upload_status=invalidated`** on the portal (IGVFFI1734TEID, IGVFFI6646SQPQ, IGVFFI6683VISA,
  IGVFFI7347KQYP). The bytes we uploaded are valid/complete (size-exact, gzip-clean), but IGVF may re-issue corrected versions;
  re-sync those if so. The whole dataset is still portal-`in progress` (Ruhi finalizing) — treat as provisional.
- The igvf-files permission gap (STS role can't read `s3://igvf-files`) will recur for any not-yet-released files on future
  datasets → same NRNB portal-download workaround, or wait for release + re-run STS.

- 2026-06-10: **Step 3 patch DONE** (`3_patch_gcp_files.sh`, BASE_DIR=Mac, DATE=2026_06_10). Decompressed onlist
  `IGVFFI1697XAXI.tsv.gz`→`patch/IGVFFI1697XAXI.tsv` and guide `IGVFFI8270UPKB.csv.gz`→`patch/IGVFFI8270UPKB.tsv`;
  uploaded the 2 local seqspec YAMLs → `patch/`; wrote **`sample_metadata_gcp_2026_06_10_patched.csv`** (1,503 rows;
  seqspec→patch/{rna,guide}_seqspec.yml 751/752, onlist→patch/…tsv, guide→patch/…tsv, fastqs→gs://). All 4 patch/ files
  confirmed present.
- 2026-06-10: validating the patched sheet with `validate_gcp_paths.py` (canonical, per request) → _result TBD_.

- 2026-06-10: Stage 1 validation re-confirmed on the **patched** sheet via fast listing-diff: **3,010/3,010 present, 0 missing**
  (the slow `validate_gcp_paths.py` keeps dying mid-run ~1.5k/3,010 — all-OK up to that point; use the listing-diff instead).
- 2026-06-10: **Stage 2 launched** (`4_run_CRISPR_pipeline.sh`, RUN_LABEL=cleanser_initial, -profile google -with-tower).
  Driver fixed to canonical shape (DATA_DATE=2026_06_10, CONFIG=setup/configs/…, banner) to match Hon CM.
  **First launch failed instantly** at `downloadReference`: `machine type n2-highmem-16 is not compatible with nvidia-tesla-t4`.
  **Config bug** in `…_cleanser_initial.config`: `downloadReference` was listed in BOTH the broad `n2-highmem-16` withName
  block AND its dedicated `n1-highmem-32`+T4 block (T4 requires n1). The working Huangfu/Hon configs only have it in the
  dedicated block. **Fix:** removed `downloadReference` from the broad alternation (line 177). Relaunched same label (-resume).
- 2026-06-10/11: Stage 2 hit two more config ceilings (this `cleanser_initial` config never ran end-to-end before; the
  working Huangfu/Hon configs share these latent issues but never tripped them at ~28 sub-pools — Hep's **47 sub-pools**
  are bigger). Both fixed in `…_cleanser_initial.config`; resource directives don't affect Nextflow's cache, so `-resume`
  keeps all completed mapping work each time:
  1. **Broad-block memory overflow.** `machine type n2-highmem-16 cannot satisfy ... memory 153600 MiB`. The broad withName
     block scales `memory = 50GB×attempt` (→150–200 GB on retries) but hard-pinned `machineType=n2-highmem-16` (128 GB). →
     bumped to **n2-highmem-32** (256 GB).
  2. **Boot-disk exhaustion (the real "another fail").** `anndata_concat` (merges all 47 scRNA sub-pools) writes 47
     `temp_processed/processed_N.h5ad` to the task disk before the final concat; died at file 44/47 with **"No space left
     on device"** (work-dir logs were empty because even `tee` couldn't write). NOT OOM — the memory fix had worked (ran at
     200 GB). → bumped global **`batch.bootDiskSize` 100 → 500 GB**. (Diagnosis came from `.nextflow.log` Command-output; GCS
     `.command.*` were all empty.)

## Status: Stage 1 COMPLETE — Stage-2-ready samplesheet:
`setup/samplesheets/sample_metadata_gcp_2026_06_10_patched.csv`
Hand to Stage 2 (`4_run_CRISPR_pipeline.sh` / crispr-pipeline-runner), pointing `--input` at the patched sheet.
