---
name: cnmf-runner
description: Stage 5 of the TFP3 pipeline. Run cNMF gene-program discovery via PerturbNMF (torch-cNMF inference + perturbation evaluation + U-test calibration + K-selection + per-TF/program PDFs + Excel summary). This is a TFP3-specific wrapper over the external EngreitzLab PerturbNMF skill — it handles h5mu→h5ad conversion, the three h5mu prep steps, the required local fix branch, per-dataset run directory conventions, and the compute budget. Triggers on keywords like Stage 5, cNMF, PerturbNMF, torch-cNMF, sk-cNMF, gene programs, K-selection, U-test calibration, perturbation association, h5mu prep, prepare_h5mu_for_eval, inject_umap, guide_annotation, K=200, compile excel summary, program analysis, perturbed gene analysis.
user_invocable: true
---

# Stage 5: cNMF Gene Programs (PerturbNMF wrapper)

You are an interactive assistant for running Stage 5 of the TFP3 pipeline. Stage 5 takes `inference_mudata.h5mu` from Stage 2 and produces a multi-stage cNMF analysis: K-sweep inference → perturbation/GO/geneset/trait evaluation → U-test FDR calibration → K-selection panel (group decision) → per-target PDFs → Excel summary.

**This skill delegates the stage mechanics to the external PerturbNMF skill.** It exists to encode the TFP3-specific bits that don't live upstream: input conversion, h5mu prep convention, required local fix branch, dataset/run dir layout, and the empirical compute budget.

## Pipeline position

```
[1] Portal → [2] CRISPR Nextflow → [3] QC → [4] Energy distance → [5] cNMF ← you are here
                ↓ inference_mudata.h5mu
           [TFP3 input conversion]  ← h5mu_to_perturbnmf_h5ad.py
                ↓ <run>/cnmf/<run_name>/Data/*.h5ad
           [External PerturbNMF runner takes over here]
                ↓ Stages 1 (inference) → 2a (eval) → 2b (U-test) → 3a (K-sel) → 3c (perturbed gene) → 3e (Excel)
                ↓ Stage 3b (program analysis) deferred (upstream OOM issue #7)
```

## Constants

```
EXTERNAL_SKILL (upstream): https://github.com/EngreitzLab/PerturbNMF/tree/main/.claude/skills/perturbNMF-runner
EXTERNAL_SKILL (local):   external/PerturbNMF/.claude/skills/perturbNMF-runner/SKILL.md
                          (NOTE: external/PerturbNMF/ is gitignored — only present on HPC
                           or after manual clone; see references/04-external-handoff.md)
LOCAL_CLONE_PATH:         external/PerturbNMF
REQUIRED_BRANCH:          fix/utest-oom-leak  (2 commits ahead of origin/main; see references/02-required-fixes.md)
CONVERTER:                src/tf_perturb_seq/cnmf/h5mu_to_perturbnmf_h5ad.py
PREP_SCRIPTS (legacy):    src/tf_perturb_seq/cnmf/{prepare_h5mu_for_eval,inject_umap_into_h5mu}.py
                          (use only for h5mu's produced before --compute_umap convention)
NEW CONVENTION:           pass --compute_umap to Convert_file_adata.py so UMAP propagates from Stage 1 input
DATASET_RUN_DIR:          datasets/<DS>/<RUN>/cnmf/<cnmf_run_name>/
                          ├── Data/                # h5ad input + guide_annotation.tsv
                          ├── Script/              # per-stage SLURM submit scripts
                          ├── README.md            # selected K + stage status
                          └── Result/<run_name>/   # all outputs
EXTERNAL_RESOURCES:       external/PerturbNMF/src/Stage2_Evaluation/Resources/{OpenTargets_L2G_Filtered.csv.gz, hocomoco_meme.meme}
SYNAPSE_TARGET:           syn64423137/2026_UTSW/datasets/<DS>/cnmf/   (per-run upload via Script/upload_to_synapse.py)
```

## Step 0: Identify substep + target run

Ask (or infer):

1. **What action?**
   - `prepare` — convert `inference_mudata.h5mu` → PerturbNMF AnnData input (see Step 1 below).
   - `inference` — Stage 1 (GPU torch-cNMF K-sweep). **External skill territory** — invoke that.
   - `prep-for-eval` — three TFP3-specific prep steps after Stage 1 (skip if you ran Stage 1 with `--compute_umap`).
   - `evaluate / calibrate / k-select / plot / summarize` — Stages 2/3. **External skill territory.**
   - `upload-synapse` — push the run dir to Synapse `syn64423137/2026_UTSW/`.
2. **Which dataset + run + cNMF run name?** `<DS>/<RUN>/cnmf/<cnmf_run_name>/`.
3. **Selected K?** Required for Stages 3a/3c/3e (and the post-K-sel UMAP injection if needed).

Then read the matching reference file.

| Topic | Reference file |
|---|---|
| h5mu→h5ad conversion + 3 prep steps + going-forward convention | `references/01-h5mu-prep.md` |
| Required PerturbNMF fix branch + PRs #4 #5 | `references/02-required-fixes.md` |
| Per-stage compute budget (RAM/CPU/GPU/time) | `references/03-compute-budget.md` |
| How to invoke the external perturbNMF-runner skill | `references/04-external-handoff.md` |

## Step 1: Convert MuData → PerturbNMF AnnData

Stage 2 produces `inference_mudata.h5mu` (separate `gene` + `guide` modalities). PerturbNMF expects a single AnnData with the gene matrix as `X` and guide info packed into `obsm`/`uns`.

```bash
source <REPO_ROOT>/.venv/bin/activate
python <REPO_ROOT>/src/tf_perturb_seq/cnmf/h5mu_to_perturbnmf_h5ad.py \
  --in_h5mu  <REPO_ROOT>/datasets/<DS>/<RUN>/crispr_pipeline/pipeline_outputs/inference_mudata.h5mu \
  --out_h5ad <REPO_ROOT>/datasets/<DS>/<RUN>/cnmf/<cnmf_run_name>/Data/<DS>_<RUN>_perturbnmf.h5ad
```

For the column-mapping spec (gene var must have `symbol`, guide var needs `guide_id` / target columns), see `references/01-h5mu-prep.md`.

**Going forward:** use `Convert_file_adata.py --compute_umap` (the per-dataset variant under `<run>/cnmf/<cnmf_run_name>/Script/`). This bakes a UMAP into the input AnnData so it propagates through Stage 1 inference; obviates the post-hoc `inject_umap_into_h5mu` prep step. Reference implementation: `datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/.../Script/Convert_file_adata.py`.

## Step 2: Stage 1 inference

**Always start from a known-working SLURM script and adapt — don't compose one from the pipeline's `--help`.** The pipeline has footguns (CPU fallthrough on a wrong kwarg signature, prepare-only-no-factorize mode) that the canonical scripts already side-step.

### Canonical reference scripts (use as templates)

| Scale | Reference script | Notes |
|---|---|---|
| Production (~270 k cells, e.g. DE / ESC) | [`datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/muddy_penguin/cnmf/Script/050926_HuangfuDE_20iter_5KHVG_torch_halsvar_batch.sh`](../../../datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/muddy_penguin/cnmf/Script/050926_HuangfuDE_20iter_5KHVG_torch_halsvar_batch.sh) | 128 GB / 72 h / A30 |
| Benchmark (~50–100 k cells) | [`datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/cleanser_800_mito_15pc/cnmf/Script/torch-cNMF_batch.sh`](../../../datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/cleanser_800_mito_15pc/cnmf/Script/torch-cNMF_batch.sh) | 128 GB / 48 h / GPU; produced consensus.txt May 10 2026 |
| Large (1 M+ cells, e.g. Hon CM / Hep) | Adapt the production template; bump `--mem` to 256–384 GB and `--time` to 24–36 h. |

### Required flag set (copy from the templates)

The pipeline's CPU path uses a re-dispatched `run_nmf(...)` whose kwargs **do not** match what `prepare()` builds. Always pass `--use_gpu` plus the exact kwarg set the templates carry — they are not optional:

```
--algo halsvar --mode batch
--init random --tol 1e-7 --use_gpu
--batch_max_epoch 1000 --batch_hals_max_iter 1000 --batch_hals_tol 0.005
--numiter 20 --numhvgenes 5000
--categorical_key batch --gene_names_key symbol
--sel_thresh 2.0                      # production; HTv2 used `0.2 2.0`
--K 30 50 60 80 100 200 250 300       # TFP3 production K sweep
--run_factorize --run_refit --run_compile_annotation
```

SLURM-side, mirror the templates: `--partition=carter-gpu --gres=gpu:a30:1 --cpus-per-task=4 --mem={128..384}G --time={12..72}h`, write logs to `<run>/cnmf/Result/<run_name>/Inference/logs/%j.{out,err}`.

### Known footguns (each cost ~10–20 min of compute when triggered)

1. **Missing `--use_gpu`** → `TypeError: run_nmf() got an unexpected keyword argument 'batch_max_epoch'` at the `factorize` step. The CPU path re-dispatches to the installed `nmf` package (PyPI nmf-torch 0.1.1, shimmed under `nmf_torch/`) whose `run_nmf` doesn't accept the kwargs `prepare()` saved. Do NOT try to fix this by reinstalling nmf-torch or editing the shim — add `--use_gpu`. See [[nmf-torch-env-regression]] memory for the full mis-diagnosis history.
2. **Missing `--run_factorize --run_refit --run_compile_annotation`** → pipeline does prepare only, writes `Inference.norm_counts.h5ad` + `Inference.tpm.h5ad` + `cnmf_tmp/`, prints `Pipeline finished.` in ~8 min, no spectra. Add the three flags; prepare reruns idempotently.
3. **Old pipeline path (`src/Inference/torch-cNMF/...`)** → `python: can't open file ...`. Upstream restructured to `src/Stage1_Inference/torch-cNMF/...` (May 2026). Patch: `sed -i 's|src/Inference/torch-cNMF|src/Stage1_Inference/torch-cNMF|g' <script>.sh`. Same restructure applies to Stage 2 (`src/Stage2_Evaluation/...`) and Stage 3 (`src/Stage3_Interpretation/...`).

For the per-stage compute budget see `references/03-compute-budget.md`. After a successful submit, check `squeue` after ~30 min — if `Elapsed` < 15 min in `sacct -j <jobid>`, suspect footgun #1 or #2 above and re-check the script before assuming the run is fast.

## Step 3: h5mu prep (if you didn't run --compute_umap upstream)

Three small jobs that make Stage 2/3 work on h5mu's produced before the new convention:

```bash
sbatch <run>/cnmf/<cnmf_run_name>/Script/prepare_h5mu_for_eval.sh    # obs['sample']='all' + remap NT guide_targets
sbatch <run>/cnmf/<cnmf_run_name>/Script/inject_umap_into_h5mu.sh    # scanpy UMAP from rna → obsm in both rna and cNMF modalities (selected K only)
# Generate Data/guide_annotation.tsv from ref/guide_libraries/harmonized/harmonized_guide_file_poolabcd.tsv
# (rename guide_id → guide_names)
```

All three are required for Stage 2b U-test to work. See `references/01-h5mu-prep.md`.

## Step 4: Stage 2/3 (external skill)

Each Stage is a separate external-skill invocation:

- Stage 2a (evaluation): `--stage evaluation`
- Stage 2b (U-test calibration): `--stage u-test-calibration` — **requires fix branch** (see `references/02-required-fixes.md`)
- Stage 3a (K-selection): `--stage k-selection`
- Stage 3c (perturbed gene PDFs): `--stage perturbed-gene`
- Stage 3e (Excel summary): `--stage excel-summary`
- Stage 3b (program analysis): deferred per upstream issue #7

For invocation examples see `references/04-external-handoff.md`.

## Step 5: Synapse upload

After all stages complete:

```bash
[[ -z "${SYNAPSE_AUTH_TOKEN:-}" ]] && { echo "set SYNAPSE_AUTH_TOKEN"; exit 2; }
python <run>/cnmf/<cnmf_run_name>/Script/upload_to_synapse.py \
  --dataset <DE|ESC|CM|...> \
  --upload-only-k200-h5mu       # skip non-selected K h5mu's to stay under ~8 GB/dataset
# --dry-run first to preview
```

Destination on Synapse: `syn64423137/2026_UTSW/datasets/<DS>/cnmf/` (flat — no `<cnmf_run_name>/` nesting in the uploaded layout). Schema: [docs/jamborees/2026_UTSW/schemas/cnmf.json](../../../docs/jamborees/2026_UTSW/schemas/cnmf.json).

## Important notes

- **External submodule must be on `fix/utest-oom-leak`.** `cd external/PerturbNMF && git log --oneline -3` should show commits `9665129` and `520cb15` on top of upstream `origin/main`. Both PRs (#4 and #5) are pending review; until merged, we run on our branch. See `references/02-required-fixes.md`.
- **Stage 3b is deferred.** Upstream issue [#7](https://github.com/EngreitzLab/PerturbNMF/issues/7) — precompute OOMs at production scale, per-program plotting glacially slow. Coverage gap is filled by Stage 3e Excel + Stage 3c per-TF PDFs.
- **Per-dataset symlink workaround.** Upstream issue [#6](https://github.com/EngreitzLab/PerturbNMF/issues/6) — Stage 3c expects `<run>/adata` but Stage 1 writes `<run>/Inference/adata`. Maintain `<run>/adata → Inference/adata` symlink per dataset.
- **`merge_pdfs_in_folder` hangs** on thousands of PDFs (upstream issue #8). Use `pdfunite` instead for the per-TF merged combined PDF.
- **`<run>/cnmf/<cnmf_run_name>/{Script,README.md}` is tracked**; `Data/`, `Result/` are gitignored (large). Per `[[feedback_track_only_scripts_and_readmes_in_cnmf]]`.
- **Per-K h5mu's are large** (~5–6 GB at K=200; ~30+ GB across the full K-sweep). Default Synapse upload skips all except selected-K to stay under ~8 GB/dataset.
- **Trait enrichment requires** `external/PerturbNMF/src/Stage2_Evaluation/Resources/OpenTargets_L2G_Filtered.csv.gz`. Download from [EngreitzLab/gene_network_evaluation](https://github.com/EngreitzLab/gene_network_evaluation/tree/main/smk/resources).
- **Motif enrichment is currently skipped** — needs `hg38.fa` + a cell-type-specific enhancer-gene linking file (scE2G), unavailable for HUES8.
- For the full project doc: [docs/analysis/cnmf/PerturbNMF.md](../../../docs/analysis/cnmf/PerturbNMF.md). For output schema: [cNMF_OUTPUTS.md](../../../docs/analysis/cnmf/cNMF_OUTPUTS.md).
