---
name: qc-runner
description: Stage 3 of the TFP3 pipeline. Run local QC on inference_mudata.h5mu outputs from Stage 2 — gene mapping QC, guide mapping QC, intended-target knockdown QC. Submits one SLURM array task per (run × dataset) on the UCSD carter-compute partition. Triggers on keywords like Stage 3, QC, qc, qc_array, mapping_gene, mapping_guide, intended_target, knockdown QC, AUROC, AUPRC, knee plot, guide capture, run_qc_pipeline, samples.tsv, SLURM array, carter-compute.
user_invocable: true
---

# Stage 3: Local QC (SLURM array)

You are an interactive assistant for running Stage 3 of the TFP3 pipeline. Stage 3 reads each run's `inference_mudata.h5mu` (mirrored from Stage 2) and produces per-batch and overall QC tables + plots for gene expression, guide capture, and intended-target knockdown. This runs on UCSD HPC under SLURM.

## Pipeline position

```
[1] Portal → [2] CRISPR Nextflow → [3] QC ← you are here → [4] Energy distance → [5] cNMF
                ↓ pipeline_outputs/inference_mudata.h5mu
           [QC]
                ↓ <run>/qc/{mapping_gene,mapping_guide,intended_target}/<prefix>_*.{tsv,png}
```

## Constants

```
REPO_ROOT (HPC):           /cellar/users/aklie/projects/tf_perturb_seq
REPO_ROOT (local):         /Users/adamklie/Desktop/tfp3/tf_perturb_seq
ARRAY_DRIVER:              <REPO_ROOT>/scripts/qc_array.sh         (frozen)
PER_SAMPLE_DRIVER:         <REPO_ROOT>/scripts/run_qc_pipeline.sh  (frozen)
PY_MODULES:                <REPO_ROOT>/src/tf_perturb_seq/qc/{mapping_gene,mapping_guide,intended_target}.py
VENV:                      <REPO_ROOT>/.venv     (uv sync; activated by run_qc_pipeline.sh)
SLURM_PARTITION:           carter-compute
SLURM_RESOURCES:           4 CPUs, 200G mem, 6h per task
SLURM_LOG_DIR:             /cellar/users/aklie/projects/tf_perturb_seq/scratch/qc_array_logs/
OUTPUT TIER (per run):     <BASE_DIR>/<RUN_LABEL>/qc/
```

QC runs on **HPC** (not local Mac). The driver scripts under `scripts/` are frozen — per memory `[[feedback_dataset_local_scripts]]`, copy variants into `datasets/<ds>/bin/` if you need a per-dataset tweak.

## Two-tier driver

```
qc_array.sh  (SLURM array wrapper)
   reads:   samples.tsv  (3 cols: INPUT, OUTDIR, RUN_NAME, with header)
   per row: bash run_qc_pipeline.sh --project-root ... --input ... --outdir ... --run-name ...
            └─ source .venv/bin/activate
            └─ python mapping_gene.py     → mapping_gene/
            └─ python mapping_guide.py    → mapping_guide/
            └─ python intended_target.py  → intended_target/
```

## The `samples.tsv`

3 columns, tab-separated, header row mandatory:

```tsv
input	outdir	run_name
/cellar/.../datasets/Hon_WTC11-benchmark_TF-Perturb-seq/cleanser_800/crispr_pipeline/pipeline_outputs/inference_mudata.h5mu	/cellar/.../datasets/Hon_WTC11-benchmark_TF-Perturb-seq/cleanser_800/qc	Hon_WTC11-benchmark_cleanser_800
/cellar/.../datasets/Huangfu_WTC11-benchmark_TF-Perturb-seq/cleanser_800/crispr_pipeline/pipeline_outputs/inference_mudata.h5mu	/cellar/.../datasets/Huangfu_WTC11-benchmark_TF-Perturb-seq/cleanser_800/qc	Huangfu_WTC11-benchmark_cleanser_800
...
```

For schema details and how to build this, see `references/01-samplesheet.md`. A scaffolder is provided:

```bash
python3 .claude/skills/qc-runner/scripts/build_samples_tsv.py \
  --repo-root <REPO_ROOT> \
  --output samples.tsv \
  [--datasets Hon_WTC11-benchmark_TF-Perturb-seq Huangfu_WTC11-benchmark_TF-Perturb-seq ...] \
  [--runs cleanser_800 seqspec_v3 ...]
```

## Step 0: Identify the run set

Ask (or infer):

1. **Which datasets + runs?** A single (dataset, run) pair, or a sweep across multiple. The skill scales the same way — N rows in `samples.tsv`, `sbatch --array=1-N`.
2. **Are Stage 2 outputs on HPC?** The `INPUT` column must be a path readable from carter-compute. If `inference_mudata.h5mu` is only on GCS, mirror to HPC first (see `references/01-samplesheet.md` §"Mirroring from GCS to HPC").
3. **Dry-run first?** `qc_array.sh ... --dry-run` runs the array but each task prints its `python ...` invocation instead of executing.

Then read the matching reference file.

| Topic | Reference file |
|---|---|
| samples.tsv schema + builder | `references/01-samplesheet.md` |
| SLURM array invocation | `references/02-slurm-array.md` |
| Output schema + interpretation | `references/03-outputs.md` |
| What the three modules actually compute | `references/04-module-internals.md` |

## Step 1: Build `samples.tsv`

Either by hand (small N) or with the builder:

```bash
python3 .claude/skills/qc-runner/scripts/build_samples_tsv.py \
  --repo-root /cellar/users/aklie/projects/tf_perturb_seq \
  --output /cellar/users/aklie/projects/tf_perturb_seq/scratch/qc_array_logs/samples_<date>.tsv
```

By default it scans every `datasets/*/<run>/crispr_pipeline/pipeline_outputs/inference_mudata.h5mu` and emits one row per find. Filter with `--datasets` and/or `--runs` to scope.

Verify the file:

```bash
head -3 samples_<date>.tsv
wc -l samples_<date>.tsv   # N = lines minus 1 (header)
```

## Step 2: Dry-run the array

```bash
N=$(($(wc -l < samples_<date>.tsv) - 1))   # exclude header
sbatch --array=1-${N} \
  <REPO_ROOT>/scripts/qc_array.sh \
  samples_<date>.tsv \
  <REPO_ROOT> \
  --dry-run
```

Check the resulting `scratch/qc_array_logs/qc_array.<JOBID>_*.out` files — each should show the three `python ...` invocations that would have run.

## Step 3: Real submission

```bash
sbatch --array=1-${N} \
  <REPO_ROOT>/scripts/qc_array.sh \
  samples_<date>.tsv \
  <REPO_ROOT>
```

Watch progress:

```bash
squeue -u $USER -t RUNNING,PENDING --states=all
sacct -j <JOBID> --format=JobID,JobName,State,Elapsed,MaxRSS
tail -f /cellar/users/aklie/projects/tf_perturb_seq/scratch/qc_array_logs/qc_array.<JOBID>_<TASKID>.out
```

For the full SLURM cookbook (re-running a single task, partition fallbacks, memory bumps), see `references/02-slurm-array.md`.

## Step 4: Sanity-check outputs

For each `(dataset, run)` row:

```bash
RUN_QC=<BASE_DIR>/<RUN_LABEL>/qc
ls ${RUN_QC}/mapping_gene/${RUN_NAME}_gene_*
ls ${RUN_QC}/mapping_guide/${RUN_NAME}_guide_*
ls ${RUN_QC}/intended_target/${RUN_NAME}_intended_target_*
```

Each module produces a `_metrics.tsv` (the headline numbers) and at least one plot. The most informative cross-sample comparison is **intended_target's AUROC** — it answers "do the guides knock down their target genes?"

```bash
# Quick AUROC table across a sweep:
for f in datasets/*/*/qc/intended_target/*_intended_target_metrics.tsv; do
  echo "$f"
  awk -F'\t' 'NR==2{print "  AUROC=" $X "  AUPRC=" $Y}' "$f"  # adjust column indices to schema
done
```

For full output interpretation, see `references/03-outputs.md`.

## Important notes

- **HPC-only.** The driver scripts hard-code carter-compute paths (`/cellar/users/aklie/...`) and assume a SLURM environment. Don't try to run them on a Mac.
- **QC outputs are gitignored.** `**/qc/**/*.{tsv,png,pdf,h5ad}` per `.gitignore`. They're regenerable from `inference_mudata.h5mu`.
- **Per-task resources (4 CPUs / 200 GB / 6h)** are sized for production datasets with the full TF library. Benchmark datasets fit in <30 min but pay no penalty for the over-allocation.
- **Engreitz bandaid.** `mapping_guide.py` copies `gene` batch labels onto `guide` obs when there's zero overlap (RNA/guide libraries in different IGVF accessions). If you see `n_batches=1` in a guide metrics table where the gene side has multiple batches, the bandaid fired — see `references/04-module-internals.md`.
- **`intended_target.py` uses `trans_per_guide_results`, not cis.** This is deliberate — including NTC guides as negative controls is needed for AUROC. The intended-target table excludes NTC rows but the metrics file reports overall + per-batch knockdown counts using the full distribution.
- The samples.tsv is **not** tracked. Build it on demand under `scratch/qc_array_logs/`.
- For the underlying design rationale, see `src/tf_perturb_seq/qc/plan.md`.
