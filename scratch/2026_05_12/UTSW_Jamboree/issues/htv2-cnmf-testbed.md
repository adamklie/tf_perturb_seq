# Issue 5 — HTv2 cNMF testbed verification → production launch

**Status**: 🔄 SLURM job 10577039 running on `carter-gpu-02` (NVIDIA A30); persistent monitor `bk2ncl52e` armed. Will fire when it leaves RUNNING.

**Owner of fix**: us. Once the testbed completes cleanly, launch production cNMF on the datasets that have full CRISPR bundles.

## TL;DR

HTv2 (Gersbach WTC11 benchmark, ~37k cells × 13k genes) is being used as a small testbed to verify the cNMF pipeline structure end-to-end on UCSD nrnb before we burn GPU-hours on production datasets (Huangfu DE/ESC are ~7× larger). Hon's `030726_20iter_5KHVG_torch_halsvar_batch_e7` run is the parameter "guiding light" — we mirrored those params exactly.

## What we set up to get the testbed running

This took 4 sbatch attempts to debug. Captured here so we don't repeat for production runs.

| Try | SLURM ID | Outcome | Root cause | Fix |
|---|---|---|---|---|
| 1 | 10576011 | Failed in 8s | Apptainer container's `cnmf 1.7.0` declares `nmf-torch` as a dep but it's missing inside the image | Drop the container; use the project venv at `/cellar/users/aklie/projects/tf_perturb_seq/.venv` (has `torch_cnmf` + `nmf_torch` already pip-installed) |
| 2 | 10576878 | Failed in 9s | `source .venv/bin/activate` from SLURM resolved Python to base conda's `/cellar/users/aklie/opt/miniconda3/bin/python` — venv didn't actually activate | Use absolute path `/cellar/users/aklie/projects/tf_perturb_seq/.venv/bin/python` |
| 3 | 10576918 | Failed in 14s | Argparse error: `--batch_max_iter` was renamed to `--batch_max_epoch` in the new pipeline; `--use_gpu` was missing from our script | Rename flag + add `--use_gpu` |
| 4 | 10577039 | RUNNING (current) | — | — |

## Current run details

| | |
|---|---|
| SLURM ID | 10577039 |
| Node | carter-gpu-02 (NVIDIA A30, 24 GB) |
| Run name | `050926_HTv2_20iter_5KHVG_torch_halsvar_batch` |
| Output dir | `/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/PerturbNMF/Result/050926_HTv2_20iter_5KHVG_torch_halsvar_batch/` |
| Time limit | 48h (carter-gpu has no time cap) |
| Memory | 128G (using ~6 GB so far — very comfortable) |
| Persistent monitor | `bk2ncl52e` (will ping when state ≠ RUNNING) |

### Hon-matched params (the "guiding light" set)

```
--algo halsvar --mode batch --init random --tol 1e-7 --use_gpu
--batch_max_epoch 1000 --batch_hals_max_iter 1000 --batch_hals_tol 0.005
--numiter 20 --numhvgenes 5000
--sel_thresh 0.2 2.0 --categorical_key batch --gene_names_key symbol
--K 5 6 7 8 9 10 11 12 13 14 15 17 19 21 23 25 27 30 35 40 45 50 55 60 70 80 90 100 150 200
--run_factorize --run_refit --run_compile_annotation --run_diagnostic_plots
```

### Latest progress snapshot (from earlier check-in)

- 6h 35m elapsed; 10 of 30 K values in progress (k=5–15, 17, 19, 21, 23 with `.npz` snapshots in `cnmf_tmp/`)
- 484 spectra `.npz` snapshots saved
- GPU: 28–40% utilization, 1.3 GB / 24 GB memory, 42°C
- Memory: 6.3 GB / 503 GB used
- Output dir: 1.0 GB

Realistic ETA: ~12–24 h more for factorization (higher K values are slower); then refit + compile + diagnostic plots add a couple more hours. Well within the 48h time limit.

## Acceptance criteria for the testbed

- [ ] Job exits with `State=COMPLETED` (not FAILED, OOM, or TIMEOUT).
- [ ] Output directory layout matches the structure documented in [`docs/analysis/cNMF_OUTPUTS.md`](../../../analysis/cNMF_OUTPUTS.md) — i.e., flat per-(k,dt) outputs (`gene_spectra_score`, `gene_spectra_tpm`, `usages.consensus`, `clustering.k_<X>.dt_<Y>.png`) + `Eval/<k>_<dt>/` subdirs + `adata/cNMF_<k>_<dt>.h5mu` per-(k,dt) MuDatas.
- [ ] At least one of each per-k file class is present and non-trivially sized (e.g. `gene_spectra_score.k_50.dt_2_0.txt` is a parseable TSV with k=50 program rows).
- [ ] Diagnostic plots written to `Plot/k_selection_*/`.

## After the testbed succeeds — production launch

Datasets ready to launch immediately (full CRISPR bundles already on Synapse):
- [ ] **Huangfu HUES8 Definitive Endoderm** — input MuData on Synapse [`syn74834952`](https://www.synapse.org/Synapse:syn74834952) (or use HPC local at `runs/muddy_penguin/pipeline_dashboard/inference_mudata.h5mu`).
- [ ] **Huangfu HUES8 Embryonic Stem Cell** — input MuData on Synapse [`syn74835010`](https://www.synapse.org/Synapse:syn74835010).

Datasets gated on other issues:
- Hon WTC11 Cardiomyocyte — [Issue #2](hon-cm-crispr-bundle.md) (Weizhou's CRISPR bundle).
- Gersbach WTC11 Hepatocyte — [Issue #3](gersbach-hep-deliverables.md) (Sara's deliverables).
- Engreitz WTC11 Endothelial — [Issue #4](engreitz-no-data.md) (no portal data).

Per-dataset launch will follow the HTv2 pattern: copy `Convert_file_adata.py` + `torch-cNMF_batch.sh` from HTv2's `PerturbNMF/Script/` to each dataset's `PerturbNMF/Script/`, change paths, sbatch.

## Pointers

| Object | Path |
|---|---|
| Run scripts | `/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/PerturbNMF/Script/{Convert_file_adata.py,torch-cNMF_batch.sh}` (also synced to local repo `datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/PerturbNMF/Script/`) |
| Pipeline source | `external/PerturbNMF/src/Stage1_Inference/torch-cNMF/Slurm_Version/torch_cnmf_inference_pipeline.py` |
| cNMF schema | [`schemas/cnmf.json`](../schemas/cnmf.json) |
| cNMF analysis walkthrough | [`docs/analysis/cNMF_OUTPUTS.md`](../../../analysis/cNMF_OUTPUTS.md), [`docs/analysis/cNMF.md`](../../../analysis/cNMF.md) |
| HPC log dir | `/cellar/.../HTv2/PerturbNMF/Result/050926_HTv2_20iter_5KHVG_torch_halsvar_batch/logs/` |
| Hon "guiding light" reference run | HPC `/cellar/users/aklie/projects/tf_perturb_seq/datasets/Hon_WTC11-benchmark_TF-Perturb-seq/PerturbNMF/Result/030726_20iter_5KHVG_torch_halsvar_batch_e7/` — used to verify the schema in [`schemas/cnmf.json`](../schemas/cnmf.json) |
