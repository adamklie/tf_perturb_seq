# Hon WTC11 Cardiomyocyte — cNMF

**Status**: ⏳ Not run yet. Gated on the full CRISPR pipeline bundle from **Weizhou** (see [Issue 2](../../../issues/hon-cm-crispr-bundle.md)) — once `pipeline_info/` is added to the canonical bundle and the source MuData is settled, this can launch.

## Plan when unblocked

Once the CRISPR bundle is canonical, follow the HTv2 testbed pattern (see [Issue 5](../../../issues/htv2-cnmf-testbed.md) for the full recipe):

1. Copy [`Convert_file_adata.py`](../../Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/cnmf/README.md) and `torch-cNMF_batch.sh` from the HTv2 testbed under `datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/PerturbNMF/Script/` to this dataset's `PerturbNMF/Script/` (HPC).
2. Update paths + `RUN_NAME` (suggested: `<DATE>_HonCM_20iter_5KHVG_torch_halsvar_batch`).
3. Set `--counts_fn` to the inference MuData (currently on Synapse [`syn74522725`](https://www.synapse.org/Synapse:syn74522725)).
4. `sbatch`.

Use the same Hon-mirroring params (`halsvar` / `batch` / 20 iter / 5K HVG / density thresholds 0.2 + 2.0 / `categorical_key=batch` / benchmark k list).

## Synapse target (when run completes)

`2026_UTSW/datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/cnmf/<run_name>/` — to be populated by `scripts/mirror_cnmf_outputs.py` (TBD; not yet built — see [`schemas/cnmf.json`](../../../schemas/cnmf.json) `bundle_inclusion_rule` for the curation spec).

## Schema + walkthrough

- Machine-readable schema: [`schemas/cnmf.json`](../../../schemas/cnmf.json)
- Analysis-level walkthrough: [`docs/analysis/cNMF_OUTPUTS.md`](../../../../../analysis/cNMF_OUTPUTS.md), [`docs/analysis/cNMF.md`](../../../../../analysis/cNMF.md)
