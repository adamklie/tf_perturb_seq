---
name: crispr-pipeline-runner
description: Stage 2 of the TFP3 pipeline. Launch the IGVF CRISPR_Pipeline (Nextflow) on Google Cloud Batch using the patched samplesheet from Stage 1 and a per-run Nextflow .config. Triggers on keywords like Stage 2, CRISPR pipeline, CRISPR_Pipeline, Nextflow, nextflow, nf-core, Google Batch, GCP Batch, -profile google, run pipeline, launch pipeline, sceptre, perturbo, cleanser, guide assignment, inference_mudata, pipeline_dashboard, RUN_LABEL, 4_run_CRISPR_pipeline.
user_invocable: true
---

# Stage 2: CRISPR Pipeline (Nextflow on GCP Batch)

You are an interactive assistant for running Stage 2 of the TFP3 pipeline. Stage 2 takes the patched samplesheet from Stage 1 plus a dataset-specific Nextflow `.config` and runs the IGVF CRISPR_Pipeline on Google Cloud Batch, writing outputs to GCS.

## Pipeline position

```
[1] IGVF portal + GCS staging
   ↓ setup/samplesheets/sample_metadata_gcp_<date>_patched.csv
[2] CRISPR Pipeline (Nextflow on GCP Batch)  ← you are here
   ↓ gs://igvf-pertub-seq-pipeline-data/<DATASET>/<DATA_DATE>/outs/<RUN_LABEL>/...
[3] QC  →  [4] Energy distance  →  [5] cNMF
```

## Constants

```
GCP_PROJECT:               igvf-pertub-seq-pipeline
GCS_BUCKET:                igvf-pertub-seq-pipeline-data
GCS_REGION:                us-central1
PIPELINE_REPO:             https://github.com/IGVF/CRISPR_Pipeline
PIPELINE_LOCAL (typical):  /Users/adamklie/Desktop/tfp3/CRISPR_Pipeline
NEXTFLOW_VERSION (pinned): 26.03.2-edge   (export NXF_VER)
NEXTFLOW_PROFILE:          google
TOWER_URL:                 https://tower.nf
GCS_LAYOUT:                gs://<bucket>/<DATASET>/<DATA_DATE>/outs/<RUN_LABEL>/...
```

Service-account JSON for GCS auth lives outside the repo (typically `~/Desktop/tfp3/igvf-pertub-seq-pipeline-*.json`) — exported via `GOOGLE_APPLICATION_CREDENTIALS`.

## The 4 knobs that differ per run

| Knob | What it controls | Examples |
|---|---|---|
| `DATA_DATE` | Which Stage 1 upload date to consume (`sample_metadata_gcp_<DATA_DATE>_patched.csv`) | `2026_04_15`, `2026_03_11` |
| `RUN_LABEL` | Nextflow run name + GCS output subdir | `seqspec_v3`, `cleanser_800`, `scrublet_on_sceptre_800` |
| `SAMPLE_METADATA` | Path to the patched CSV from Stage 1 | `<BASE_DIR>/setup/samplesheets/sample_metadata_gcp_<DATA_DATE>_patched.csv` |
| `CONFIG` | Per-run Nextflow `.config` | `<BASE_DIR>/setup/configs/<DATASET>_<RUN_LABEL>.config` |

`DATASET_NAME`, `BASE_DIR`, `PIPELINE_PATH`, `OUTDIR` derive from those four.

## Step 0: Identify the substep and target run

Ask (or infer):

1. **What action?**
   - `launch` — start a new run (most common; see Step 1 below)
   - `resume` — re-run with `-resume` after a transient failure (see Step 2 below)
   - `monitor` — tail logs / check Tower / `gsutil ls` outputs (see `references/03-monitoring.md`)
   - `inspect-outputs` — once a run completes, what's in `outs/<RUN_LABEL>/` and what to mirror locally (see `references/04-outputs.md`)
   - `scaffold` — generate a fresh `setup/configs/<DATASET>_<RUN_LABEL>.config` and a `4_run_*.sh` for a new run label (see "Scaffolding a new run" below)
2. **Which dataset?** Path under `datasets/`.
3. **Which `RUN_LABEL` + `DATA_DATE`?** Either an existing one (re-run) or a new one (new parameter set).

Then read the matching reference file before issuing commands.

| Topic | Reference file |
|---|---|
| Nextflow config params | `references/01-config-spec.md` |
| Driver script anatomy | `references/02-driver-script.md` |
| Monitoring + failure modes | `references/03-monitoring.md` |
| Output schema + local mirroring | `references/04-outputs.md` |

## Step 1: Launch a new run

Prereqs (verify before launching — a failed launch on GCP burns time):

- Stage 1 complete: `<BASE_DIR>/setup/samplesheets/sample_metadata_gcp_<DATA_DATE>_patched.csv` exists.
- `setup/configs/<DATASET>_<RUN_LABEL>.config` exists and matches the dataset's chemistry/protocol.
- `CRISPR_Pipeline` is cloned and on the expected branch (see `references/02-driver-script.md`).
- Env: `NXF_VER` pinned, `GOOGLE_APPLICATION_CREDENTIALS` set, `gcloud auth application-default login` done, `TOWER_ACCESS_TOKEN` set (optional but strongly recommended for visibility).

Run:
```bash
RUN_IN_BACKGROUND=true bash <BASE_DIR>/setup/scripts/4_run_CRISPR_pipeline.sh
```

This invokes:
```
nextflow run main.nf \
  -profile google \
  -c <CONFIG> \
  --input <SAMPLE_METADATA> \
  --outdir <OUTDIR> \
  -resume \
  -with-tower
```

For each per-flag detail, see `references/02-driver-script.md`.

## Step 2: Resume a partially-completed run

`-resume` is already in the script — re-running the same `RUN_LABEL` resumes from cache. Two failure modes need different handling:

- **Transient (spot preemption, network, retryable exit codes 137/143/50001/50002/50003/50006):** the Nextflow `errorStrategy` retries automatically (up to `maxRetries=3`). For full-run failures, just re-run the same script — `-resume` picks up where it died.
- **Bad config / bad input:** stop, fix the `.config` or samplesheet, **bump `RUN_LABEL`** to a new value so you don't poison the cache, and launch fresh.

For failure-mode-to-fix mapping, see `references/03-monitoring.md`.

## Scaffolding a new run

Generate a fresh driver script + skeleton `.config` for a new `RUN_LABEL`:

```bash
python3 .claude/skills/crispr-pipeline-runner/scripts/scaffold_run_script.py \
  --dataset-name <DATASET_NAME> \
  --run-label <RUN_LABEL> \
  --data-date <YYYY_MM_DD> \
  --base-config <path/to/sibling-dataset.config>   # optional: copy + edit instead of starting blank
```

The scaffolder writes:
- `datasets/<DATASET>/setup/configs/<DATASET>_<RUN_LABEL>.config` — copy of `--base-config` if given, else a TFP3 default skeleton.
- `datasets/<DATASET>/setup/scripts/4_run_CRISPR_pipeline.sh` — only if missing; otherwise emit a diff-style summary of what to edit (DATA_DATE / RUN_LABEL / SAMPLE_METADATA / CONFIG).

## Handoff to Stage 3 (QC)

When the run completes, the canonical artifacts are:

```
gs://igvf-pertub-seq-pipeline-data/<DATASET>/<DATA_DATE>/outs/<RUN_LABEL>/
├── pipeline_outputs/        # inference_mudata.h5mu — Stage 3 input
├── pipeline_dashboard/      # dashboard.html + additional_qc/ (gene/guide/intended_target metrics)
├── pipeline_info/           # params_*.json + versions yml — mirror locally (small, TRACKED)
└── (many intermediate stages — see references/04-outputs.md)
```

Mirror `pipeline_info/` into the local run dir so the run-provenance commits with the repo:

```bash
mkdir -p <BASE_DIR>/<RUN_LABEL>/crispr_pipeline/pipeline_info
gsutil -m cp -r gs://.../outs/<RUN_LABEL>/pipeline_info/* <BASE_DIR>/<RUN_LABEL>/crispr_pipeline/pipeline_info/
```

Stage 3 (QC) consumes `pipeline_outputs/inference_mudata.h5mu`. See `references/04-outputs.md` for the full mirror recipe.

## Important notes

- **Pin `NXF_VER`.** Nextflow's `-edge` build self-updates aggressively; without pinning, runs are not reproducible across days. The newer drivers export `NXF_VER=26.03.2-edge`.
- **`-c` is non-optional in the current convention.** Some older driver scripts (Engreitz benchmark, Huangfu_WTC11-benchmark, Gersbach benchmarks) **omit** `-c $CONFIG` — those runs pick up pipeline defaults, which is rarely what we want. When updating an older dataset, add `-c $CONFIG`.
- **Cache poisoning.** Re-launching with the same `RUN_LABEL` but a changed `.config` can produce stale cached steps. If you change `.config`, bump `RUN_LABEL` and start fresh (Nextflow caches per work dir; the script writes work to `gs://<bucket>/work`).
- **Spot preemption is normal.** `batch.spot = true` is on by default; tasks get preempted and retry. Don't panic at warnings in the log.
- **One run = one `RUN_LABEL` = one `outs/<label>/` subdir.** Tech-benchmark datasets do 7 parameter-sweep runs; cardiomyocyte does fewer larger ones. Outputs do **not** share directories.
- **Configs live in `setup/configs/`.** Never under `<run>/crispr_pipeline/configs/`. The Hon cardio script has a known path bug pointing under `<run>/` — copy the canonical `setup/configs/` version when adapting.
- Run-provenance is `pipeline_info/params_*.json`. Folder names are misleading — always read the params JSON to know what differs between two runs.
- For the full CRISPR pipeline docs (not project-specific): [docs/analysis/crispr_pipeline/CRISPR_PIPELINE.md](../../../docs/analysis/crispr_pipeline/CRISPR_PIPELINE.md) and [CRISPR_PIPELINE_OUTPUTS.md](../../../docs/analysis/crispr_pipeline/CRISPR_PIPELINE_OUTPUTS.md).
