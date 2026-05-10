# Gersbach WTC11 benchmark HTv2 — cNMF (testbed)

**Status**: 🔄 SLURM job 10577039 running on `carter-gpu-02`. Persistent monitor `bk2ncl52e` armed.

## Purpose

Verify the cNMF pipeline structure end-to-end (script ports, env activation, container choice, output layout, mirror script wiring) before launching production cNMF on Huangfu DE + ESC. HTv2 was chosen because it's small enough to complete in hours, not days.

## Run

| | |
|---|---|
| Run name | `050926_HTv2_20iter_5KHVG_torch_halsvar_batch` |
| Output dir (HPC) | `/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/PerturbNMF/Result/050926_HTv2_20iter_5KHVG_torch_halsvar_batch/` |
| SLURM script | [`PerturbNMF/Script/torch-cNMF_batch.sh`](../../../../../datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/PerturbNMF/Script/torch-cNMF_batch.sh) |
| Convert script | [`PerturbNMF/Script/Convert_file_adata.py`](../../../../../datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/PerturbNMF/Script/Convert_file_adata.py) |
| Pipeline source | `/cellar/users/aklie/projects/tf_perturb_seq/external/PerturbNMF/src/Stage1_Inference/torch-cNMF/Slurm_Version/torch_cnmf_inference_pipeline.py` |
| Params | mirrors Hon's `030726_20iter_5KHVG_torch_halsvar_batch_e7` exactly (halsvar / batch / 20 iter / 5K HVG / density thresholds 0.2 + 2.0 / `categorical_key=batch` / benchmark k list 5–200) |

## After the testbed completes

Acceptance criteria + the production-launch plan are in [Issue 5: HTv2 cNMF testbed verification → production launch](../../../issues/htv2-cnmf-testbed.md).

## Schema + walkthrough

- Machine-readable schema: [`schemas/cnmf.json`](../../../schemas/cnmf.json)
- Analysis-level walkthrough: [`docs/analysis/cNMF_OUTPUTS.md`](../../../../../analysis/cNMF_OUTPUTS.md), [`docs/analysis/cNMF.md`](../../../../../analysis/cNMF.md)
