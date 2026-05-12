---
name: deg-calibration
description: Calibrate PerTurbo per-element differential-expression (DEG) results from the IGVF CRISPR pipeline using an empirical null distribution built from non-targeting controls (NTCs), then apply BH FDR correction. Produces calibrated trans, cis, direct-target, and trans-only result tables. Triggers on keywords like calibration, calibrate, DEG calibration, differential expression calibration, empirical null, NTC, non-targeting controls, FDR, BH correction, PerTurbo, perturbo, perturbo_trans, perturbo_cis, calibrated_trans_results, calibrated_direct_target, p-value calibration, hit calling, t-fit, ecdf, run_calibration.sh.
user_invocable: true
---

# DEG Calibration

You are an interactive assistant for calibrating PerTurbo per-element DEG results from Stage 2 (the CRISPR Nextflow pipeline). The IGVF pipeline emits posterior p-values from a negative binomial model; these are **not** well-calibrated for frequentist hit calling. Calibration derives empirical p-values by comparing each targeting element's test statistic against the NTC distribution, then BH-corrects.

## Pipeline position

```
[1] Portal → [2] CRISPR Nextflow → [3] QC → [4] Energy distance → [5] cNMF
                ↓ pipeline_outputs/perturbo_trans_per_element_output.tsv.gz
                ↓ pipeline_outputs/perturbo_cis_per_element_output.tsv.gz
                ↓ pipeline_outputs/inference_mudata.h5mu
           [Calibration]   ← you are here (parallel tier to QC / Energy distance / cNMF)
                ↓ <run>/calibration/<prefix>_calibrated_*.tsv
```

This skill is **DEG calibration** (PerTurbo per-element). Do **not** confuse with **U-test perturbation calibration** under `<run>/cnmf/Script/U-test_perturbation_calibration.sh` — that's a Stage 5 (cNMF) substep that calibrates per-perturbation NMF program scores, not per-element DEGs.

## Constants

```
REPO_ROOT (local):     /Users/adamklie/Desktop/tfp3/tf_perturb_seq
DRIVER:                <REPO_ROOT>/scripts/run_calibration.sh    (frozen — do not edit per-dataset)
PY_TOOL:               <REPO_ROOT>/src/tf_perturb_seq/inference/calibrate.py
VENV:                  <REPO_ROOT>/.venv     (uv sync provisions this; the driver activates it)
OUTPUT TIER (per run): <BASE_DIR>/<RUN_LABEL>/calibration/
```

`scripts/run_calibration.sh` at the repo root is the canonical entry point. Per memory `[[feedback_dataset_local_scripts]]`, the scripts under `scripts/` are **frozen** — if a dataset needs a variant, copy to `datasets/<ds>/bin/` and edit there, don't edit in place.

## Inputs

| Input | From | Path |
|---|---|---|
| `perturbo_trans_per_element_output.tsv.gz` | Stage 2 | `<run>/crispr_pipeline/pipeline_outputs/` (or `gs://.../outs/<run>/pipeline_outputs/`) |
| `perturbo_cis_per_element_output.tsv.gz` | Stage 2 | same dir (optional; some calls use it implicitly) |
| `inference_mudata.h5mu` | Stage 2 | same dir |

You must mirror Stage 2's `pipeline_outputs/` locally before calibrating — see `crispr-pipeline-runner` skill's `references/04-outputs.md` for the recipe.

## Outputs

Four TSVs written to `<outdir>/<prefix>_*.tsv`:

| Output | Contents |
|---|---|
| `<prefix>_calibrated_trans_results.tsv` | Master table: every (targeting element × gene) test, calibrated p, BH-adjusted q |
| `<prefix>_calibrated_direct_target_results.tsv` | Tests of each guide's direct target gene (TF on its own gene) |
| `<prefix>_calibrated_cis_results.tsv` | Cis-window tests (annotated using `--cis-window`, default 100 kb) |
| `<prefix>_calibrated_trans_only_results.tsv` | Trans tests excluding the cis window |

For the schema of each, see `references/01-inputs-outputs.md`.

## Step 0: Identify the target run

Ask (or infer):

1. **Which dataset + run?** `datasets/<DATASET>/<RUN_LABEL>/`
2. **Where are Stage 2 outputs?** Local mirror at `<BASE_DIR>/<RUN_LABEL>/crispr_pipeline/pipeline_outputs/`, or still on GCS (need to mirror first).
3. **Which null method?** `t-fit` (default; better in low-NTC regime) or `ecdf` (rank-based, distribution-free). See `references/02-null-methods.md`.
4. **Anything non-default?** Custom cis window, one-sided test, non-standard NTC label.

Then read the matching reference file.

| Topic | Reference file |
|---|---|
| Input/output schema | `references/01-inputs-outputs.md` |
| Choosing `--null-method` | `references/02-null-methods.md` |
| Full CLI option reference | `references/03-cli-options.md` |

## Step 1: Mirror Stage 2 outputs locally (if not already)

```bash
DATASET=<...>
DATA_DATE=<YYYY_MM_DD>
RUN_LABEL=<...>
BASE_DIR=<REPO_ROOT>/datasets/${DATASET}
LOCAL=${BASE_DIR}/${RUN_LABEL}/crispr_pipeline/pipeline_outputs
GCS=gs://igvf-pertub-seq-pipeline-data/${DATASET}/${DATA_DATE}/outs/${RUN_LABEL}/pipeline_outputs

mkdir -p ${LOCAL}
gsutil -m cp ${GCS}/inference_mudata.h5mu ${LOCAL}/
gsutil -m cp ${GCS}/perturbo_trans_per_element_output.tsv.gz ${LOCAL}/
gsutil -m cp ${GCS}/perturbo_cis_per_element_output.tsv.gz ${LOCAL}/
```

## Step 2: Run calibration

### Default path — submit as a SLURM job (recommended for production datasets)

Calibration loads the full MuData (16–30 GB) + the 20M-row trans TSV into memory and runs single-process for 30 min–2 h. **Do not run inline in the foreground or as a background shell task** — the harness will reap it after ~15 s. Use `sbatch`:

```bash
DS=<DATASET>
RUN=<RUN_LABEL>
REPO=<REPO_ROOT>
mkdir -p ${REPO}/scratch/calibration_logs
mkdir -p ${REPO}/datasets/${DS}/${RUN}/calibration

sbatch \
    --job-name=calib_${RUN} \
    --partition=carter-compute \
    --cpus-per-task=8 \
    --mem=256G \
    --time=06:00:00 \
    --output=${REPO}/scratch/calibration_logs/calib_${RUN}.%j.out \
    --error=${REPO}/scratch/calibration_logs/calib_${RUN}.%j.err \
    --wrap="bash ${REPO}/scripts/run_calibration.sh \
      --project-root ${REPO} \
      --trans-results ${REPO}/datasets/${DS}/${RUN}/crispr_pipeline/pipeline_outputs/perturbo_trans_per_element_output.tsv.gz \
      --mudata        ${REPO}/datasets/${DS}/${RUN}/crispr_pipeline/pipeline_outputs/inference_mudata.h5mu \
      --outdir        ${REPO}/datasets/${DS}/${RUN}/calibration \
      --prefix        ${DS}_${RUN}"
```

Track progress: `squeue -u $USER` and `tail -f ${REPO}/scratch/calibration_logs/calib_${RUN}.<jobid>.out`.

**Resource notes**:
- `--mem=256G` — MuData read in full + pandas overhead on 20M+ row TSV. Lower mem can OOM mid-load.
- `--cpus-per-task=8` — `calibrate.py` is mostly single-threaded but pandas/numpy spawn helper threads; over-allocating costs nothing.
- `--time=06:00:00` — most runs land in 30 min–2 h; 6 h is a safety buffer.
- `--partition=carter-compute` — CPU-only, no GPU needed.

### Submit all production datasets at once

When mirroring calibration across all 4 production datasets, loop over `(DATASET, RUN_LABEL)` pairs and `sbatch` each — independent jobs so they run in parallel:

```bash
declare -a DATASETS=(
  "Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq:muddy_penguin"
  "Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq:sceptre_v1"
  "Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq:2026_04_19_no_spacer"
  "Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq:sara_synapse_syn74842722"
)

for entry in "${DATASETS[@]}"; do
  DS="${entry%:*}"; RUN="${entry##*:}"
  # sbatch ... --wrap="bash run_calibration.sh ..."   (as above)
done
```

### Fallback — direct shell (dev / dry-run only)

```bash
bash ${REPO}/scripts/run_calibration.sh \
  --project-root ${REPO} \
  --trans-results ... --mudata ... --outdir ... --prefix ... \
  --dry-run    # always start with --dry-run; once happy, drop the flag
```

Inline foreground runs are fine for `--dry-run` validation but **not** for the real call — even a small dataset takes longer than the harness's bg-task reaper allows.

## Step 3: Run calibration with non-default options

The wrapper `scripts/run_calibration.sh` only exposes `--trans-results`, `--mudata`, `--outdir`, `--prefix`, `--null-method`, `--dry-run`. It does **not** expose `--cis-window`, `--non-targeting-label`, `--guide-mod-key`, or `--one-sided`. To use any of those, call `calibrate.py` directly:

```bash
source <REPO_ROOT>/.venv/bin/activate
python <REPO_ROOT>/src/tf_perturb_seq/inference/calibrate.py \
  --trans-results <...> \
  --mudata <...> \
  --outdir <...> \
  --prefix <...> \
  --null-method t-fit \
  --cis-window 250000 \
  --non-targeting-label non_targeting \
  --guide-mod-key guide \
  [--one-sided]
```

The wrapper also parses `--cis-results` but silently drops it; that flag is **not** an argument to `calibrate.py` either (cis results are inferred by the Python code from related paths/mudata). If a future PR exposes `--cis-results`, update this skill.

## Step 4: Sanity-check outputs

```bash
LS=<BASE_DIR>/<RUN_LABEL>/calibration
ls $LS/<prefix>_calibrated_*.tsv

# Quick histogram of calibrated p-values (uniform under null = good calibration)
head -1 $LS/<prefix>_calibrated_trans_results.tsv
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="empirical_pvalue")c=i; next} c{print $c}' \
    $LS/<prefix>_calibrated_trans_results.tsv | \
  python3 -c "import sys,numpy as np; p=np.array([float(x) for x in sys.stdin if x.strip()]); print(f'n={len(p)} mean={p.mean():.3f} median={np.median(p):.3f}')"

# Hit counts at common FDR cutoffs
for q in 0.01 0.05 0.1; do
  awk -F'\t' -v q=$q 'NR>1 && $NF<q' \
      $LS/<prefix>_calibrated_trans_results.tsv | wc -l \
      | xargs -I{} echo "  q < $q: {} tests"
done
```

If the p-value histogram is **strongly non-uniform under the null** (heavy enrichment near 0 or 1), revisit `--null-method` and check the NTC count — see `references/02-null-methods.md`.

## Important notes

- Calibration **must** run on per-element TSVs from `pipeline_outputs/`, not on the per-guide TSVs or on MuData `.uns` tables (which are per-guide). Using the wrong inputs silently produces nonsense.
- `<run>/calibration/` is **gitignored** per `docs/data/DATA.md`. Calibrated TSVs are regenerable from the pipeline outputs; don't commit them.
- Default `--non-targeting-label` is `"non_targeting"` and `--guide-mod-key` is `"guide"`. These match the IGVF pipeline output schema; only override if a dataset has a non-standard label in `guide.var["type"]`.
- The `.venv` the driver activates is provisioned by `uv sync` at the repo root. If activation fails, run `uv sync` in `<REPO_ROOT>` first.
- For Stage 5 (cNMF) U-test perturbation calibration, use the `cnmf-runner` skill (planned), **not** this one.
- For the design rationale of empirical-null calibration, see `src/tf_perturb_seq/inference/plan.md`.
