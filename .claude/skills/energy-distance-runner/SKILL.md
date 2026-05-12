---
name: energy-distance-runner
description: Stage 4 of the TFP3 pipeline. Run the containerized energy-distance pipeline (Chikara Takeuchi's) on a dataset's inference_mudata.h5mu — gRNA filtering, permutation test for per-target energy distance vs non-targeting, optional pairwise target matrix + 2D embedding. SLURM GPU on UCSD carter-gpu using apptainer. Triggers on keywords like Stage 4, energy distance, edistance, e-distance, energy_dist_pipeline, energy_dist_TFperturb, Chikara, Takeuchi, gRNA filtering, DISCO test, non-targeting permutation, target-by-target matrix, edist, edistance clustering, 5_run_energy_distance, carter-gpu, apptainer, edist_pipeline.sif.
user_invocable: true
---

# Stage 4: Energy Distance (containerized, SLURM GPU)

You are an interactive assistant for running Stage 4 of the TFP3 pipeline. Stage 4 takes a single dataset's `inference_mudata.h5mu` (Stage 2 output) and computes per-perturbation energy distance vs the non-targeting baseline, with permutation-based significance. Optionally builds a pairwise target × target distance matrix and a 2D embedding for downstream clustering.

The pipeline is Chikara Takeuchi's; we wrap it in `scripts/run_energy_distance_pipeline.sh` and submit one SLURM job per dataset.

## Pipeline position

```
[1] Portal → [2] CRISPR Nextflow → [3] QC → [4] Energy distance ← you are here → [5] cNMF
                                              ↓ pval_edist_full.csv, target_by_target_matrix.csv,
                                                edist_embedding_info.csv
```

## Constants

```
REPO_ROOT (HPC):           /cellar/users/aklie/projects/tf_perturb_seq
RUNNER:                    <REPO_ROOT>/scripts/run_energy_distance_pipeline.sh   (frozen)
PER_DATASET_DRIVER:        datasets/<DS>/<RUN>/energy_distance/scripts/5_run_energy_distance.sh   (SLURM wrapper)
PIPELINE_BIN (submodule):  external/energy_dist_pipeline/bin/
PREPROCESS:                src/tf_perturb_seq/crispr_pipeline/preprocess_mudata_local.py
CONTAINER (apptainer):     /cellar/users/aklie/opt/containers/edist_pipeline.sif
CONTAINER_IMAGE:           docker://docker.io/takechikara/energy_distance_env:latest
SLURM_PARTITION:           carter-gpu
SLURM_RESOURCES:           1 GPU (a30), 8 CPU, 200G, 2 days
```

The energy-distance pipeline runs on **HPC only** — needs apptainer + GPU. Don't try locally.

## Pipeline stages (inside the container)

| Step | Script | Inputs | Outputs |
|---|---|---|---|
| 0 | `preprocess_mudata_local.py` | `inference_mudata.h5mu` | `preprocessed.h5ad`, `gRNA_dict.pickle`, `pca_dataframe.pickle`, `annotation_table.csv` |
| 1 | `1_filtering_gRNA.py` | step 0 outputs + `config1_2.json` | `targeting_outlier_table.csv`, `non_targeting_outlier_table.csv`, `discordance_gRNA.csv` |
| 2 | `2_e_distance_nontargeting.py` | step 1 outputs | `pval_edist_full.csv` — **the headline output** |
| 2.1 | `2_1_Plot_figure.py` | step 2 outputs | diagnostic PNGs |
| 3 (skipped by default) | `3_e_distance_among_regions.py` | + `config3.json` cutoffs | `target_by_target_matrix.csv`, `edist_embedding_info.csv` |

Step 3 is intentionally skipped because it depends on choosing a significance + edist cutoff. Run it after eyeballing `pval_edist_full.csv`.

## Step 0: Identify dataset + run + MuData source

Ask (or infer):

1. **Which dataset + run?** `datasets/<DS>/<RUN>/`.
2. **Where is `inference_mudata.h5mu`?** One of three sources:
   - `--mudata-path` (already on HPC) — fastest
   - `--gcs-mudata-path gs://...` — runner downloads to `OUTPUT_DIR/`
   - `--synapse-id synXXXXX` — runner downloads via synapseclient (needs `SYNAPSE_AUTH_TOKEN`)
3. **Run step 3 too?** Default is skip; pick cutoffs first.

Then read the matching reference file.

| Topic | Reference file |
|---|---|
| Inputs + config1_2.json knobs | `references/01-inputs-config.md` |
| SLURM driver + GPU resources | `references/02-slurm-driver.md` |
| Outputs + step 3 (cutoffs) | `references/03-outputs-step3.md` |
| Container + submodule quirks | `references/04-container-quirks.md` |

## Step 1: Submit a run

Most datasets already have a `5_run_energy_distance.sh` SLURM wrapper. Submit it:

```bash
cd /cellar/users/aklie/projects/tf_perturb_seq
sbatch datasets/<DS>/<RUN>/energy_distance/scripts/5_run_energy_distance.sh
```

The wrapper:
- Requests `--partition=carter-gpu --gres=gpu:a30:1 --cpus-per-task=8 --mem=200G --time=2-00:00:00`.
- `module load apptainer`.
- Calls `scripts/run_energy_distance_pipeline.sh` with the right source flag and OUTPUT_DIR.

For a new dataset/run with no wrapper yet, copy from a sibling:

```bash
cp datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/2026_04_19_no_spacer/energy_distance/scripts/5_run_energy_distance.sh \
   datasets/<NEW_DS>/<NEW_RUN>/energy_distance/scripts/5_run_energy_distance.sh
# Edit: --job-name, log paths, OUTPUT_DIR, source flag (--synapse-id / --gcs-mudata-path / --mudata-path)
```

## Step 2: Monitor

```bash
squeue -u $USER
sacct -j <JOBID> --format=JobID,State,Elapsed,MaxRSS,NodeList
tail -f <OUTPUT_DIR>/logs/<JOBID>.out
nvidia-smi -L   # confirm GPU on the assigned node
```

Typical wall time: 4–12 hours for benchmark datasets, 12–36 hours for production datasets with the full TF library. Wall is dominated by step 2's permutation budget (default 20 backgrounds × 1000 permutations each).

## Step 3: Inspect step 2 outputs

```bash
cd <OUTPUT_DIR>
head -3 pval_edist_full.csv
ls -1 *.png    # step 2.1 diagnostic plots
```

`pval_edist_full.csv` columns (subset): `intended_target_name`, `edist_mean`, `edist_std`, `pvalue`, plus per-permutation breakdowns. See `references/03-outputs-step3.md` for the full schema.

## Step 4: Run step 3 separately (optional)

After picking cutoffs (typically `pvalue<0.05`, `edist>0.5`), edit `<OUTPUT_DIR>/config3.json` then submit step 3:

```bash
# Re-run the runner with --run-step3 (skips steps 0/1/2 if outputs exist; idempotent)
bash /cellar/users/aklie/projects/tf_perturb_seq/scripts/run_energy_distance_pipeline.sh \
  --mudata-path <OUTPUT_DIR>/inference_mudata.h5mu \
  --output-dir  <OUTPUT_DIR> \
  --run-step3
```

This is idempotent — preprocess + steps 1 / 2 / 2.1 outputs already exist and get skipped. Only step 3 runs fresh.

Step 3 outputs:
- `target_by_target_matrix.csv` — pairwise e-distance between significant targets
- `edist_embedding_info.csv` — 2D t-SNE coordinates

## Important notes

- **Submodule must be initialized.** `external/energy_dist_pipeline` is pinned. After `git submodule update --init external/energy_dist_pipeline`, verify `external/energy_dist_pipeline/bin/` exists with the step scripts.
- **Step 1 filename typo.** At the currently pinned submodule commit, the script is `1_filtereing_gRNA.py` (typo). Upstream main has `1_filtering_gRNA.py`. The runner checks both — if you bump the submodule past the rename commit, no action needed.
- **DON'T pip-install muon into `/tmp/muon_deps`.** Inside the container, use the pre-installed `muon` from `/app/.venv`. Installing a newer muon pulls numpy ≥2.0 whose pickled artifacts can't be deserialized by the container's numpy 1.26.4 in steps 1/2/2.1 — silent corruption.
- **Output dir convention is in flux.** Old layout: `datasets/<DS>/results/energy_distance/<RUN>/`. New (canonical, per `docs/data/DATA.md`): `datasets/<DS>/<RUN>/energy_distance/`. Existing SLURM wrappers still use the old `results/...` path; new ones should use the canonical layout. The runner doesn't care — it writes wherever `--output-dir` points.
- **`5_run_energy_distance.sh` lives under `<RUN>/energy_distance/scripts/`**, not `<DS>/setup/scripts/`. This is per-run because the OUTPUT_DIR + source flags differ per run.
- **`config1_2.json` and `config3.json` are auto-generated** by the runner in `<OUTPUT_DIR>` each invocation. You can edit them between step 2 and step 3 (cutoffs); the runner won't overwrite if you re-invoke with `--run-step3` because it only regenerates when not present — actually it does regenerate; see `references/04-container-quirks.md` for the workaround.
- **`SYNAPSE_AUTH_TOKEN`** must already be in the environment if using `--synapse-id`. Per memory `[[synapse_project]]`, it's in `~/.zshrc`.
- **`<run>/energy_distance/configs/` is tracked**; `image/`, `logs/`, `*.csv`, `*.tsv`, `*.h5mu`, `*.pickle` are gitignored.
- Detailed pipeline docs: [docs/analysis/energy_dist/ENERGY_DISTANCE.md](../../../docs/analysis/energy_dist/ENERGY_DISTANCE.md), output schema: [ENERGY_DISTANCE_OUTPUTS.md](../../../docs/analysis/energy_dist/ENERGY_DISTANCE_OUTPUTS.md).
