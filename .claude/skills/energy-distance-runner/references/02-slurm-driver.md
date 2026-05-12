# SLURM driver (`5_run_energy_distance.sh`)

The per-dataset SLURM submitter wraps `scripts/run_energy_distance_pipeline.sh` with the right resources, env, and source flag.

## Canonical shape (Hon CM, 2026_04_19_no_spacer)

```bash
#!/usr/bin/env bash
#SBATCH --job-name=edist_hon_cm
#SBATCH --partition=carter-gpu
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --gres=gpu:a30:1
#SBATCH --output=/cellar/.../results/energy_distance/<RUN>/logs/%j.out
#SBATCH --error=/cellar/.../results/energy_distance/<RUN>/logs/%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=aklie@ucsd.edu

set -euo pipefail

module load apptainer            # not on default PATH on nrnb compute nodes

SYNAPSE_ID="synXXXXXXXX"          # or use --gcs-mudata-path / --mudata-path
OUTPUT_DIR="/cellar/.../results/energy_distance/<RUN>"
mkdir -p "${OUTPUT_DIR}/logs"

[[ -z "${SYNAPSE_AUTH_TOKEN:-}" ]] && { echo "SYNAPSE_AUTH_TOKEN not set" >&2; exit 2; }

bash /cellar/users/aklie/projects/tf_perturb_seq/scripts/run_energy_distance_pipeline.sh \
  --synapse-id "$SYNAPSE_ID" \
  --output-dir "$OUTPUT_DIR"
```

## Per-dataset edits

| Field | What it controls |
|---|---|
| `--job-name` | Short tag for `squeue` (`edist_<dataset_short>`) |
| `--output` / `--error` | Per-job log paths; reuse `${OUTPUT_DIR}/logs/%j.{out,err}` |
| Source line | One of `SYNAPSE_ID=...` / `GCS_MUDATA=gs://...` / `MUDATA_PATH=/cellar/...` |
| `OUTPUT_DIR` | `datasets/<DS>/<RUN>/energy_distance/` (canonical) or `datasets/<DS>/results/energy_distance/<RUN>/` (legacy; in use today) |
| Runner flag | Pass `--synapse-id "$SYNAPSE_ID"` / `--gcs-mudata-path "$GCS_MUDATA"` / `--mudata-path "$MUDATA_PATH"` |

## SLURM resources (don't usually tune)

- `--partition=carter-gpu` — the partition with A30s on UCSD nrnb. Fallback: any GPU partition with an A30/L4/L40s.
- `--gres=gpu:a30:1` — one A30. Step 2's PyTorch perm test fits comfortably.
- `--cpus-per-task=8` — used by some parallel ops in step 2.
- `--mem=200G` — production datasets with full TF library peak around 150 GB; benchmark datasets fit in 64 GB.
- `--time=2-00:00:00` — 48 h cap. Most runs finish in 6–18 h.

If the wall time looks like it'll bust:
- Confirm GPU is actually attached: `nvidia-smi` at top of the job (the canonical wrapper does this).
- Halve `permute_per_bg` in `config1_2.json` from 1000 → 500. Doubles speed at the cost of p-value resolution.

## Source flags — when to use which

| Source | Flag | Typical use |
|---|---|---|
| Already on HPC | `--mudata-path /cellar/...` | After Stage 2 mirror or a prior run already downloaded it |
| GCS | `--gcs-mudata-path gs://igvf-pertub-seq-pipeline-data/...` | Most datasets (canonical IGVF CRISPR pipeline outputs land here) |
| Synapse | `--synapse-id synXXXXX` | When the dataset hasn't been fully mirrored to GCS (currently: Hon CM) |

Idempotent: the runner checks `<OUTPUT_DIR>/inference_mudata.h5mu` first and skips download if present. Wall time savings on a re-submit are substantial.

## Submitting

```bash
cd /cellar/users/aklie/projects/tf_perturb_seq
sbatch datasets/<DS>/<RUN>/energy_distance/scripts/5_run_energy_distance.sh
squeue -u $USER
```

## Monitoring

```bash
JID=<jobid>
sacct -j $JID --format=JobID,JobName,State,Elapsed,MaxRSS,MaxVMSize,ExitCode
tail -f /cellar/.../<RUN>/logs/${JID}.out

# Periodic GPU check (need to ssh to the compute node)
ssh <node> "nvidia-smi --query-gpu=memory.used,memory.total,utilization.gpu --format=csv"
```

## Resubmitting selectively

The runner's idempotency makes resubmits cheap:

- If only step 3 needs to run, the runner skips steps 0/1/2 (outputs exist) and runs step 3 alone — fast (minutes for benchmark; an hour or two for production).
- If step 2 needs to re-run with different perm budget: delete `pval_edist_full.csv` and re-submit; runner detects it's missing and re-runs from step 2.
- If everything needs to re-run: delete `<OUTPUT_DIR>` (or move it aside) and submit fresh.

## Failure modes

| Symptom | Cause | Fix |
|---|---|---|
| `module: command not found` | `module load` not available | Check that the SLURM job is actually on nrnb |
| `apptainer: command not found` | Module isn't loaded | Add `module load apptainer` (canonical wrapper has it) |
| `nvidia-smi: command not found` or `No devices found` | GPU not actually attached | Check `--gres=gpu:a30:1` and squeue for "(Resources)" reason |
| Job dies at step 0 with KeyError | MuData missing required columns | See `01-inputs-config.md` §"Required Stage 2 columns" |
| Pickle deserialization error in step 2 | numpy version mismatch (see `04-container-quirks.md`) | Don't pip-install muon into /tmp/muon_deps |
| Step 1 script not found | Submodule pin or rename mismatch | Runner falls back from `1_filtering_gRNA.py` to `1_filtereing_gRNA.py`; if both missing, bump submodule |
| OOM on GPU in step 2 | Target_cell_num_max too high for available VRAM | Halve `target_cell_num_max` in `config1_2.json` |
