# Huangfu HUES8 Embryonic Stem Cell — cNMF

**Status**: ⏳ Not run yet. **Unblocked** — full canonical CRISPR bundle is on Synapse ([`syn74835010`](https://www.synapse.org/Synapse:syn74835010)) and the inference MuData is ready. Will launch as soon as the HTv2 testbed (job 10577039) verifies the pipeline structure cleanly. See [Issue 5](../../../issues/htv2-cnmf-testbed.md).

## Plan

Follow the HTv2 testbed pattern:

1. Copy [`Convert_file_adata.py`](../../Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/cnmf/README.md) and `torch-cNMF_batch.sh` from the HTv2 testbed `PerturbNMF/Script/` to this dataset's `PerturbNMF/Script/` (HPC).
2. Update paths + `RUN_NAME` (suggested: `<DATE>_HuangfuESC_20iter_5KHVG_torch_halsvar_batch`).
3. Set `--counts_fn` to the inference MuData. Source options:
   - GCS: `gs://igvf-pertub-seq-pipeline-data/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/2026_04_13/outs/sceptre_v1/pipeline_dashboard/inference_mudata.h5mu`
   - Synapse: under [`syn74835010`](https://www.synapse.org/Synapse:syn74835010)
4. `sbatch`.

Use the same Hon-mirroring params (`halsvar` / `batch` / 20 iter / 5K HVG / density thresholds 0.2 + 2.0 / `categorical_key=batch` / benchmark k list).

## Synapse target (when run completes)

`2026_UTSW/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/cnmf/<run_name>/` — to be populated by `scripts/mirror_cnmf_outputs.py` (TBD).

## Schema + walkthrough

- Machine-readable schema: [`schemas/cnmf.json`](../../../schemas/cnmf.json)
- Analysis-level walkthrough: [`docs/analysis/cNMF_OUTPUTS.md`](../../../../../analysis/cNMF_OUTPUTS.md), [`docs/analysis/cNMF.md`](../../../../../analysis/cNMF.md)
