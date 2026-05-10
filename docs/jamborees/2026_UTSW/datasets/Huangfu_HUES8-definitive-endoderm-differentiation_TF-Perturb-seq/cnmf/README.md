# Huangfu HUES8 Definitive Endoderm — cNMF

**Status**: ⏳ **Pre-staged, ready to launch.** SLURM script written + scp'd to HPC. Holding submission until the HTv2 testbed (job 10577039) verifies the pipeline structure cleanly. See [Issue 5](../../../issues/htv2-cnmf-testbed.md).

## Pre-staged

| | |
|---|---|
| Run name | `050926_HuangfuDE_20iter_5KHVG_torch_halsvar_batch` |
| Input `.h5ad` | `/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/PerturbNMF/Data/DE_muddy_penguin_perturbnmf.h5ad` (1.3 GB; pre-converted from MuData) |
| SLURM script | [`PerturbNMF/Script/050926_HuangfuDE_20iter_5KHVG_torch_halsvar_batch.sh`](../../../../../datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/PerturbNMF/Script/050926_HuangfuDE_20iter_5KHVG_torch_halsvar_batch.sh) |
| Params | mirrors Hon's `030726_20iter_5KHVG_torch_halsvar_batch_e7` (halsvar / batch / 20 iter / 5K HVG / density thresholds 0.2 + 2.0 / `categorical_key=batch` / benchmark k list 5–200) |

## To launch (when HTv2 testbed verifies)

```bash
ssh aklie@nrnb-login.ucsd.edu
sbatch /cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/PerturbNMF/Script/050926_HuangfuDE_20iter_5KHVG_torch_halsvar_batch.sh
```

Time limit set to 72h (Huangfu DE is ~7× larger than HTv2 testbed; expect ~24-48h actual wall time).

## Synapse target (when run completes)

`2026_UTSW/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/cnmf/<run_name>/` — populated by [`scripts/mirror_cnmf_outputs.py`](../../../scripts/mirror_cnmf_outputs.py) (curation rule: `schemas/cnmf.json`).

## Schema + walkthrough

- Machine-readable schema: [`schemas/cnmf.json`](../../../schemas/cnmf.json)
- Analysis-level walkthrough: [`docs/analysis/cNMF_OUTPUTS.md`](../../../../../analysis/cNMF_OUTPUTS.md), [`docs/analysis/cNMF.md`](../../../../../analysis/cNMF.md)
