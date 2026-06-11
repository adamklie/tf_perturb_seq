# Handoff — Gersbach WTC11 hepatocyte TF-Perturb-seq onboarding

Written 2026-06-11. Companion docs (read these first — not duplicated here):
- [GCS_UPLOAD_RUNBOOK.md](GCS_UPLOAD_RUNBOOK.md) — full Stage 1 + Stage 2 run log, every fix with rationale, redo/resume steps.
- [PORTAL_STATUS_2026_06_10.md](PORTAL_STATUS_2026_06_10.md) — portal completeness snapshot + what's ruled out.
- Tracking issue: [#28](https://github.com/adamklie/tf_perturb_seq/issues/28).

## Where things stand

**Stage 1 (portal → GCS): COMPLETE + validated.** 3,010/3,010 objects in
`gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/2026_06_10/`.
Stage-2-ready samplesheet: `setup/samplesheets/sample_metadata_gcp_2026_06_10_patched.csv`
(1,503 rows = 751 scRNA + 752 gRNA; one S10 scRNA R1 was absent on the portal and is correctly dropped).
Seqspecs were corrected + validated against the real reads (UMI=12, R1=28, guide construct lands at R2 63–82) —
`setup/seqspec/{rna,guide}_seqspec.yml`.

**Stage 2 (CRISPR Nextflow on GCP Batch): RUNNING / being -resumed by Adam.**
- Adam launches it himself in a `screen` session (`screen -S hep_crispr`), foreground, via
  `setup/scripts/4_run_CRISPR_pipeline.sh` (RUN_LABEL=`cleanser_initial`, `-profile google -c … -resume -with-tower`).
- Outputs: `gs://…/2026_06_10/outs/cleanser_initial/`. Visible on tower.nf.
- Three config fixes already applied to `setup/configs/…_cleanser_initial.config` (this config had never run end-to-end;
  the 47-sub-pool scale tripped latent issues the ~28-sub-pool Huangfu/Hon configs never hit). See runbook for detail:
  1. `downloadReference` removed from the broad `withName` block (it double-matched → got n2 machineType with a T4 GPU; T4 needs n1).
  2. broad-block `machineType` n2-highmem-16 → **n2-highmem-32** (retry-scaled memory hit 150–200 GB > 128 GB).
  3. global `batch.bootDiskSize` 100 → **500 GB** (`anndata_concat` over 47 sub-pools filled the 100 GB disk — "No space left on device" at temp file 44/47; NOT OOM).
- Last known state: `anndata_concat` was the failing step; the resume should clear it and proceed to
  guide assignment (cleanser) → inference (perturbo + sceptre) → dashboard.

## What the next agent should do

1. **Monitor the resume (GCP-side, non-intrusive — don't touch the screen process or `.nextflow.log` while it runs unless asked):**
   - `gcloud batch jobs list --location=us-central1 --project=igvf-pertub-seq-pipeline` (RUNNING/FAILED states; auth as the SA — see below).
   - `gsutil ls -r gs://…/2026_06_10/outs/cleanser_initial/**` to see published stages.
2. **If another task fails, use the established playbook** (it worked 3×): read the failed Batch job's `.command.*` in its
   `gs://…/work/<hash>/` dir; if those are empty, read `CRISPR_Pipeline/.nextflow.log` "Command output/error" block (the
   authoritative source — GCS work-dir logs were empty on disk-full). Most likely next ceilings are the across-all-sub-pools
   aggregation steps (`mudata_concat`, `inference_mudata`, `mergeMudata`) or perturbo inference — same fix class (size the
   right `withName` block's memory/disk). Resource directives DON'T invalidate Nextflow's cache, so always `-resume`.
3. On completion: mirror `pipeline_info/` locally and hand `pipeline_outputs/inference_mudata.h5mu` to Stage 3 (QC). See
   the crispr-pipeline-runner skill's `references/04-outputs.md`.

## Parked / open items

- **NRNB seqspec-validation harness** — DRAFTED, NOT RUN, under `scripts/seqspec_validation/` (DESIGN.md +
  build/submit/collect scripts). Goal: run the crispr_validator on a representative seqspec per dataset as SLURM array jobs
  on NRNB `carter-compute`. Only 3 of 11 datasets have seqspecs (Gersbach Hep, Hon CM, Hon benchmark). Open Qs in its
  DESIGN.md: IGVF creds on the NRNB login node, whether to `git pull` NRNB first, centralized vs per-dataset wrapper.
  Validator is cloned at `external/crispr_validator/` (gitignored). Known validator quirks already handled in the drafts:
  guide_design `.csv.gz` vs `.tsv.gz` download bug; default `--feature-sample-reads 100000` is too slow (use ~5000).
- **Open decision: sceptre vs cleanser for guide assignment.** Current run uses `cleanser` + `direct-capture` (correct for
  Hep's validated direct-capture chemistry; mirrors Sara's canonical `params_2026-04-01`). Adam questioned whether it should
  be sceptre. Sceptre ↔ CROP-seq (e.g. Huangfu); switching would contradict the validated chemistry, so needs a rationale.
  If a sceptre comparison is wanted, do it as a SEPARATE run (`RUN_LABEL=sceptre_initial` + own config), not by editing
  `cleanser_initial` (cache + label hygiene). Note: inference already runs perturbo+sceptre (`INFERENCE_method='default'`).

## Environment notes (no secrets here)
- GCP: project `igvf-pertub-seq-pipeline`, bucket `igvf-pertub-seq-pipeline-data`. The bucket is writable only by the
  service account `adamklie@igvf-pertub-seq-pipeline.iam.gserviceaccount.com` — `gcloud config set account …` to it
  (the default `adamklie13@gmail.com` gets 403). Same SA is credentialed on NRNB.
- Pipeline clone: `~/Desktop/tfp3/CRISPR_Pipeline` (branch `main`). `NXF_VER=26.03.2-edge` (pinned in the driver).
- `TOWER_ACCESS_TOKEN` is in `~/.zshrc` (a screen shell inherits it). Service-account JSON path is set in the driver.
  These secrets are intentionally NOT reproduced here.
- IGVF portal staging recovery used `aklie@nrnb-login.ucsd.edu` for the 9 GCS-unreachable reads (4 were
  `upload_status=invalidated` with stale md5 — bytes verified valid by size+gzip; IGVF may re-issue them — re-sync if so).

## Suggested skills for the next session
- `crispr-pipeline-runner` — Stage 2 monitoring, resume, failure diagnosis, output mirroring (what the active run is).
- `qc-runner` — Stage 3, once `inference_mudata.h5mu` is produced.
- `igvf-portal-staging` — only if re-syncing the 4 invalidated reads or onboarding another dataset.
- The NRNB harness work is plain scripting (no skill); start from `scripts/seqspec_validation/DESIGN.md`.
