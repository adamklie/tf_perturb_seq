# `4_run_CRISPR_pipeline.sh` anatomy

Driver script at `datasets/<DATASET>/setup/scripts/4_run_CRISPR_pipeline.sh`. The canonical (newer) shape — used by Hon cardio, Huangfu HUES8 DE/ESC, Gersbach hepatocyte.

## Required env (top of script)

```bash
export NXF_VER=26.03.2-edge
export GOOGLE_APPLICATION_CREDENTIALS=/path/to/igvf-pertub-seq-pipeline-<keyid>.json
# Optional but recommended:
export TOWER_ACCESS_TOKEN=...   # from tower.nf
```

- **`NXF_VER`** — pin the Nextflow version. Without this, `-edge` self-updates and runs lose reproducibility across days.
- **`GOOGLE_APPLICATION_CREDENTIALS`** — service-account JSON for `igvf-pertub-seq-pipeline`. Lives outside the repo. Without it, the GCS work dir is unwritable and the run fails at step 1.
- **`TOWER_ACCESS_TOKEN`** — if set, `tower {}` in the config enables Tower; `-with-tower` then streams run state to tower.nf. Without it the flag is silently ignored.

## Script variables

```bash
DATASET_NAME=...                 # e.g. Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq
BASE_DIR=<repo>/datasets/${DATASET_NAME}
DATA_DATE=YYYY_MM_DD             # matches Stage 1 upload date
RUN_LABEL=...                    # per-run label, e.g. seqspec_v3, cleanser_800
SAMPLE_METADATA=$BASE_DIR/setup/samplesheets/sample_metadata_gcp_${DATA_DATE}_patched.csv
PIPELINE_PATH=<path to IGVF/CRISPR_Pipeline clone>
CONFIG=$BASE_DIR/setup/configs/${DATASET_NAME}_${RUN_LABEL}.config
OUTDIR=gs://igvf-pertub-seq-pipeline-data/${DATASET_NAME}/${DATA_DATE}/outs/${RUN_LABEL}
LOG_FILE=$BASE_DIR/logs/${DATASET_NAME}_${RUN_LABEL}_$(date +%Y%m%d_%H%M%S).log
RUN_IN_BACKGROUND=${RUN_IN_BACKGROUND:-false}
```

**Per-run edits:** `DATA_DATE`, `RUN_LABEL`. Sometimes `SAMPLE_METADATA` if you're using a non-standard patched CSV variant (e.g. Hon cardio's `_patched_v3.csv`).

## The Nextflow command

```bash
cd $PIPELINE_PATH
nextflow run main.nf \
    -profile google \
    -c $CONFIG \
    --input $SAMPLE_METADATA \
    --outdir $OUTDIR \
    -resume \
    -with-tower
```

Flag-by-flag:

| Flag | Required? | Notes |
|---|---|---|
| `-profile google` | yes | Selects the Google Batch executor + retry strategy + spot config from `profiles { google { ... } }` |
| `-c $CONFIG` | yes (canonical) | Loads dataset-specific params. Older scripts omit this — they pick up pipeline defaults, rarely correct |
| `--input $SAMPLE_METADATA` | yes | Stage 1's patched CSV |
| `--outdir $OUTDIR` | yes | GCS path; **include `RUN_LABEL`** so runs don't collide |
| `-resume` | almost always | Re-uses cached work for unchanged tasks. Safe to leave on |
| `-with-tower` | optional | Streams to tower.nf if `TOWER_ACCESS_TOKEN` is set; no-op otherwise |
| `-w gs://.../work` | implicit | Comes from `workDir` in config |
| `-r <branch>` | sometimes | Pin pipeline branch/tag; relevant if `PIPELINE_PATH` is a checkout that you don't manage |
| `-params-file <yaml>` | not used | TFP3 keeps everything in `.config`. The "loads dataset-specific parameters from YAML" comment in old scripts is misleading |

## Foreground vs background

```bash
bash <DATASET_DIR>/setup/scripts/4_run_CRISPR_pipeline.sh                         # foreground
RUN_IN_BACKGROUND=true bash <DATASET_DIR>/setup/scripts/4_run_CRISPR_pipeline.sh  # nohup + PID printed
```

Background is the right default — runs are typically 4–8 hours wall, sometimes longer for full TF library production datasets.

## Logs

Local: `$BASE_DIR/logs/${DATASET_NAME}_${RUN_LABEL}_<TS>.log`. Created by `mkdir -p $BASE_DIR/logs` in the script. These are big (10s–100s of MB after a full run) and gitignored.

GCS work dir: `gs://<bucket>/work/<hash>/` — each task has its own subdirectory with `.command.sh`, `.command.out`, `.command.err`, `.exitcode`. Useful when a specific task fails (see `references/03-monitoring.md`).

## The `CRISPR_Pipeline` clone

`PIPELINE_PATH` points at a local checkout of https://github.com/IGVF/CRISPR_Pipeline. Nextflow runs `main.nf` from there; the work dir is on GCS so the local checkout just provides `main.nf` and the modules.

- Keep this checkout outside the repo (e.g. `~/Desktop/tfp3/CRISPR_Pipeline`). Don't commit it as a submodule — versions are pinned via `-r` or just by current HEAD.
- Pull updates before launching: `cd $PIPELINE_PATH && git pull`. The IGVF repo evolves; running today's container set against last-month's `main.nf` can surface mismatches.
- Branch matters. Check `git status` / `git log -1` before launching. The standard branch is `main` unless the team has explicitly pinned to a feature branch.

## Older script shape (Engreitz, Huangfu_WTC11-benchmark, Gersbach benchmarks)

Those scripts omit:
- `export NXF_VER` (no version pin)
- `export GOOGLE_APPLICATION_CREDENTIALS` (assumes ambient auth)
- `-c $CONFIG` (uses pipeline defaults)
- `-resume`

When updating these to launch a new run, port them to the canonical shape above. Don't run them as-is.

## Common edits when scaffolding a new run

1. Bump `RUN_LABEL` to a fresh value.
2. Update `DATA_DATE` if consuming a re-uploaded samplesheet.
3. Ensure `$BASE_DIR/setup/configs/${DATASET_NAME}_${RUN_LABEL}.config` exists.
4. Spot-check `OUTDIR` includes the new `RUN_LABEL` (the script auto-derives this; just verify).
5. If the chemistry-specific flags in the config need adjusting (sweep), do it in the `.config`, not in the driver.
