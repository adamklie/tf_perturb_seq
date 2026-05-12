# SLURM array invocation

The QC driver is a SLURM array script. One array task = one row of `samples.tsv`.

## Submitting

```bash
N=$(($(wc -l < samples.tsv) - 1))   # exclude header
sbatch --array=1-${N} \
  /cellar/users/aklie/projects/tf_perturb_seq/scripts/qc_array.sh \
  /path/to/samples.tsv \
  /cellar/users/aklie/projects/tf_perturb_seq \
  [--dry-run]
```

Positional args to `qc_array.sh`:

1. **Samples TSV path** — absolute path; must be readable from compute nodes.
2. **Project root** — repo root containing `src/tf_perturb_seq/qc/`, `.venv/`, and `scripts/run_qc_pipeline.sh`.
3. **`--dry-run`** (optional) — passed through to `run_qc_pipeline.sh`. The three `python ...` invocations are printed but not executed.

## Built-in resource request

The `#SBATCH` header in `qc_array.sh` is:

```
#SBATCH -J qc_array
#SBATCH -c 4
#SBATCH --mem=200G
#SBATCH --partition=carter-compute
#SBATCH -t 06:00:00
#SBATCH -o .../scratch/qc_array_logs/qc_array.%A_%a.out
#SBATCH -e .../scratch/qc_array_logs/qc_array.%A_%a.err
```

These are sized for production datasets (full TF library × many cells). Override per-submission if needed:

```bash
# Smaller resources for benchmark dataset
sbatch --array=1-${N} --mem=64G --cpus-per-task=2 -t 02:00:00 \
  qc_array.sh samples.tsv $REPO_ROOT

# Different partition (if carter-compute is queued long)
sbatch --array=1-${N} --partition=general-compute \
  qc_array.sh samples.tsv $REPO_ROOT
```

## Log paths

```
SLURM_LOG_DIR=/cellar/users/aklie/projects/tf_perturb_seq/scratch/qc_array_logs/
  qc_array.<JOBID>_<TASKID>.out    # stdout
  qc_array.<JOBID>_<TASKID>.err    # stderr
```

The script creates `SLURM_LOG_DIR` if missing. **`%A_%a`** in the SBATCH `-o`/`-e` expands to `<JOBID>_<TASKID>`, one log pair per task.

## Monitoring

```bash
# Live queue
squeue -u $USER -t RUNNING,PENDING --array

# Per-task status table
sacct -j <JOBID> --format=JobID,JobName,State,Elapsed,MaxRSS,MaxVMSize,ExitCode

# Tail one task's stdout
tail -f /cellar/users/aklie/projects/tf_perturb_seq/scratch/qc_array_logs/qc_array.<JOBID>_<TASKID>.out

# Which tasks failed?
sacct -j <JOBID> --format=JobID,State | grep -v COMPLETED
```

## Re-running a single failed task

If task `<TASKID>` failed (OOM, bad input, etc.) but others completed, re-run just that index:

```bash
sbatch --array=<TASKID> \
  qc_array.sh samples.tsv $REPO_ROOT
```

You can also pass a comma-separated list or range: `--array=3,7,12` or `--array=10-20`.

## Common failure modes

| Symptom | Cause | Fix |
|---|---|---|
| Task exits immediately with `ERROR: No TSV line N` | Array index doesn't match a TSV row (too high) | Recompute N as lines - 1 (header) |
| `ERROR: Input does not exist:` | `inference_mudata.h5mu` missing or unreadable | Mirror from GCS; verify path is on a shared FS, not local to a node |
| Task gets `MaxRSS` close to `--mem` and exits 137 (OOM) | Real OOM | Bump `--mem` to 300G or 400G; rare unless dataset is unusually large |
| `ImportError` from one of the Python modules | `.venv` wasn't built or is stale | On HPC: `cd $REPO_ROOT && uv sync` |
| `KeyError: 'num_expressed_genes'` in `mapping_gene` | Older pipeline output without that column | The module has a fallback to `log1p_n_genes_by_counts` — if that's also missing, regenerate `inference_mudata.h5mu` from Stage 2 |
| Task hangs at startup | NFS lag finding the venv | Rerun; if persistent, mount-check |
| All tasks queue forever | `carter-compute` partition full | Override `--partition=general-compute` |

## After a successful run

Verify each output dir got populated:

```bash
awk -F'\t' 'NR>1 {print $2}' samples.tsv | while read OUTDIR; do
  echo "$OUTDIR:"
  ls -1 "$OUTDIR/mapping_gene/" "$OUTDIR/mapping_guide/" "$OUTDIR/intended_target/" 2>/dev/null | wc -l
done
```

Expect roughly 12+ files per OUTDIR (3 subdirs × 3–5 outputs each). For per-module output detail, see `03-outputs.md`.

## Quirks

- The `--dry-run` flag, when passed to `qc_array.sh`, is forwarded to `run_qc_pipeline.sh`, which `printf`'s the would-be commands instead of running them. Use it on a fresh `samples.tsv` to verify pathing before burning compute.
- `qc_array.sh` exits 1 if any of: TSV missing, `run_qc_pipeline.sh` not executable, line N empty after the array index. Failures surface in the `.err` file with explicit `ERROR:` prefixes.
- The `SLURM_LOG_DIR` is hardcoded in `qc_array.sh` to `/cellar/users/aklie/...`. Different HPC user → patch the script (copy to `datasets/<ds>/bin/` and edit there per the frozen-scripts rule).
