# Huangfu HUES8 Definitive Endoderm — cNMF

**Status**: 🔄 Mirroring the existing `042926_huangfu_de_torchcnmf_KskillA` run to Synapse (in flight 2026-05-10). Selected k = 200, density threshold = 2.0.

## What's being mirrored

| | |
|---|---|
| Run name | `042926_huangfu_de_torchcnmf_KskillA` |
| Source | HPC: `/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/PerturbNMF/Result/042926_huangfu_de_torchcnmf_KskillA/` (37 GB) |
| Selected k | **200** |
| Density threshold | **2.0** (single — this run did not sweep dt=0.2 + 2.0) |
| Mirror target | `2026_UTSW/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/cnmf/042926_huangfu_de_torchcnmf_KskillA/` |
| Bundle size | ~9 GB after curation |
| Mirror command | `nohup .venv/bin/python docs/jamborees/2026_UTSW/scripts/mirror_cnmf_outputs.py --dataset Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq --source-dir /cellar/.../042926_huangfu_de_torchcnmf_KskillA --selected-k 200 &` |
| Log on HPC | `docs/jamborees/2026_UTSW/mirror_huangfu_de_cnmf.log` |

## ⚠ Layout deviations from the schema

This run was generated with a slightly older PerturbNMF tooling that produces a layout different from the Hon benchmark / canonical schema. The mirror script auto-detects and uploads as-is:

- File prefix is `Inference.` (instead of `<run_name>.`) — e.g. `Inference.gene_spectra_score.k_200.dt_2_0.txt`
- Flat files live inside `Inference/` subdir (instead of run-dir top level)
- Folder is `Evaluation/` (instead of `Eval/`)
- `Annotation/` lives at `Inference/Annotation/`
- Single density threshold `dt_2_0` only (no `dt_0_2` since `--sel_thresh 2.0` was the run config)

Per-(k,dt) Eval/ folder + per-K loadings are still present (8 K values: 30, 50, 60, 80, 100, 200, 250, 300 — a curated subset rather than the full benchmark sweep).

## Run parameters (differ from Hon's "guiding light")

(from `042926_huangfu_de_torchcnmf_KskillA_inference.sh`)

| Param | Huangfu DE 042926 | Hon benchmark (reference) |
|---|---|---|
| K sweep | 30, 50, 60, 80, 100, 200, 250, 300 (8 values) | 5–200 (30 values) |
| numiter | 10 | 20 |
| numhvgenes | 2000 | 5000 |
| sel_thresh | 2.0 (single) | 0.2 + 2.0 |
| seed | 14 | (default 123) |
| tol | 1e-4 | 1e-7 |
| algo / mode / categorical_key | halsvar / batch / batch | same |

The differences are **acceptable for the jamboree** — Adam confirmed these runs are "pretty much done minus a few issues." For cross-dataset comparisons, anchor on shared k values (50, 60, 80, 100, 200) and on `gene_spectra_score` z-scored loadings (which are normalization-invariant).

## Schema + walkthrough

- Machine-readable schema: [`schemas/cnmf.json`](../../../schemas/cnmf.json)
- Analysis-level walkthrough: [`docs/analysis/cNMF_OUTPUTS.md`](../../../../../analysis/cNMF_OUTPUTS.md), [`docs/analysis/cNMF.md`](../../../../../analysis/cNMF.md)

## Pre-staged Hon-mirroring run (held)

A new Hon-mirroring SLURM script (`050926_HuangfuDE_20iter_5KHVG_torch_halsvar_batch.sh`) is also pre-staged on HPC if we ever want to re-run with Hon's exact params. **Holding submission** since the 042926 run is sufficient. See [Issue 5](../../../issues/htv2-cnmf-testbed.md).
