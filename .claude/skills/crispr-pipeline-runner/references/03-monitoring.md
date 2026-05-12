# Monitoring + failure modes

A typical TFP3 CRISPR run is 4–8 hours wall on GCP Batch. Watch the run in **three** places.

## 1. Local log file

```bash
tail -f <BASE_DIR>/logs/<DATASET>_<RUN_LABEL>_*.log
```

What to look for:

- `executor >  google-batch (N)` — N tasks submitted to Batch. Should climb over the first 5–10 minutes.
- `[xx/abc123] process > NAME (M)` — progress per process. Anything in `[100%]` is done.
- `WARN: Spot VM preempted task ...` — normal; Nextflow retries.
- `ERROR ~ ...` followed by `WORK_DIR: gs://.../work/<hash>/...` — terminal task failure. Inspect the work dir for `.command.err`.

## 2. Nextflow Tower (recommended)

If `TOWER_ACCESS_TOKEN` was exported and `-with-tower` is in the command:

- Visit https://tower.nf → your workspace → recent runs.
- One row per `RUN_LABEL`. Click in to see per-task state, retries, resource use.
- Use the **Resources** tab to spot OOM (memory usage at 100% of allocation = bump the `process { withName: X { memory = ... } }` block in the next run).

If you forgot to set the token, Tower is silently absent — there's no error.

## 3. GCS output dir

```bash
gsutil ls gs://igvf-pertub-seq-pipeline-data/<DATASET>/<DATA_DATE>/outs/<RUN_LABEL>/
```

Each completed stage publishes a directory there (see `references/04-outputs.md` for the full layout). When `pipeline_outputs/inference_mudata.h5mu` and `pipeline_dashboard/dashboard.html` exist, the run is done.

```bash
# Quick "is it done" check
gsutil ls gs://.../outs/<RUN_LABEL>/pipeline_outputs/inference_mudata.h5mu && echo DONE
```

## Failure modes

### Transient (auto-retry, no action)

| Symptom | Cause | What to do |
|---|---|---|
| `task.exitStatus = 137` | OOM-kill, container exceeded memory | Auto-retried; Nextflow scales memory by `attempt` |
| `task.exitStatus = 143` | SIGTERM (preemption) | Auto-retried |
| `50001..50006` | Google Batch transient errors | Auto-retried up to `maxRetries=3` |
| Tower shows status `RUNNING` but log is quiet | Spot preemption + retry; backoff | Wait |

### Run-killers (require fixing then re-launch)

| Symptom | Cause | Fix |
|---|---|---|
| `Cannot invoke method toLowerCase() on null object` | Samplesheet malformed (TSV instead of CSV, blank `file_modality`, wrong header) | Re-check Stage 1 output |
| `No such variable: <param>` | Typo or missing param in `.config` | Diff against the IGVF pipeline's `nextflow.config` |
| `Process ... terminated with an error exit status (1)` after `maxRetries` | Real bug in the task; check `WORK_DIR/.command.err` | Read err, fix config or input |
| `403 Forbidden` on a `gs://` path | Wrong service-account key in `GOOGLE_APPLICATION_CREDENTIALS` | Re-issue / re-point |
| `Pipeline execution stopped` immediately at launch | Java version (need 17+) or `NXF_VER` mismatch with `main.nf` syntax | `java -version`; check Nextflow release notes |
| `Cannot find ... reference ...` in `downloadReference` | GENCODE / IGVF reference URL stale | Update `REFERENCE_gtf_download_path` in `.config` |
| Run hangs in `downloadReference` for hours | `aria2c` container can't reach EBI mirror | Switch the URL to a different mirror or use `REFERENCE_gtf_local_path` |
| Run completes but `pipeline_dashboard/` is empty | `createDashboard` failed silently | Check Tower / log for `createDashboard` warnings |

### When to bump `RUN_LABEL` vs `-resume` in place

- **Same RUN_LABEL + `-resume`**: change is to driver env (`NXF_VER`, auth), pipeline checkout, or pure infrastructure (a process's resource block). Cached results stay valid.
- **Bump RUN_LABEL**: change is to any `params {}` in the `.config`, or to the input samplesheet. Old cache becomes invalid; using `-resume` would silently keep stale results.

When in doubt, bump.

## Inspecting a failed task's work dir

When the log says:

```
ERROR ~ Error executing process > 'guide_assignment_sceptre (xx)'

  Caused by:
    Process `guide_assignment_sceptre (xx)` terminated with an error exit status (1)

  Work dir:
    gs://igvf-pertub-seq-pipeline-data/work/aa/bbbbbbbb...
```

Inspect:

```bash
WORK=gs://igvf-pertub-seq-pipeline-data/work/aa/bbbbbbbb
gsutil ls $WORK
gsutil cat $WORK/.command.sh         # the actual command nextflow ran
gsutil cat $WORK/.command.err        # stderr
gsutil cat $WORK/.command.log        # nextflow wrapper log
gsutil cat $WORK/.exitcode           # numeric exit code
```

`.command.sh` is reproducible — you can copy it to a local machine and re-run for debugging if the inputs are small.

## Killing a run

If the script is running in the foreground: Ctrl-C. Nextflow handles SIGINT by canceling Batch tasks.

If `RUN_IN_BACKGROUND=true`:

```bash
ps aux | grep nextflow                       # find PID
kill -INT <PID>                              # graceful (cancels Batch jobs)
gcloud batch jobs list --location=us-central1 | grep <pattern>   # verify cleanup
```

For a truly stuck run, cancel jobs by job name:

```bash
gcloud batch jobs delete <job-name> --location=us-central1
```

## Cost monitoring

`batch.spot = true` keeps cost down but doesn't make runs free. Order-of-magnitude:

- Benchmark dataset (50-gene library, ~5 lanes): tens of $.
- Full TF library production (~2000 targets, ~30 lanes): low hundreds of $.

Watch `gcloud billing` if a run loops on retries (often a config bug causing the same task to fail-retry indefinitely). Spot retries are capped at `batch.maxSpotAttempts = 5`; non-spot retries at `maxRetries = 3`. A pathological case shouldn't exceed those.
