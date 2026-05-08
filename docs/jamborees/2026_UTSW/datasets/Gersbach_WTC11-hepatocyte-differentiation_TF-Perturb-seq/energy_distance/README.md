# Energy distance — Gersbach WTC11 Hepatocyte

| | |
|---|---|
| Dataset ID | `Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq` |
| Schema | [`../../../schemas/energy_distance.json`](../../../schemas/energy_distance.json) |
| Run label | TBD (canonical CRISPR pipeline run not yet picked) |
| Source MuData | Synapse [`syn74728027`](https://www.synapse.org/Synapse:syn74728027) — `gersbach_iPSC_iHep_04.06.2026_cleanser_inference_mudata.h5mu` (30.31 GB, latest cleanser run) |
| Status | Not run yet — awaiting Gersbach team's canonical run + 3-folder bundle layout |
| HPC output dir | `/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/results/energy_distance/<run>/` (TBD) |
| Synapse | _not yet created_ |

## How to run (when ready)

Dataset doesn't yet have a `5_run_energy_distance.sh`. Once the canonical inference MuData is settled:

```bash
# pull from Synapse on the HPC
sbatch /cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/5_run_energy_distance.sh
```

Runner supports `--synapse-id syn74728027`.

## Notes

- Gersbach hepatocyte uses the cleanser / direct-capture pipeline path (vs. sceptre / CROP-seq for Huangfu and Hon CM); guide assignment differences may affect e-distance comparability across datasets.
- The Synapse parent folder `syn70518849` has multiple MuData versions (v1, v2, latest cleanser) — pick a canonical one before running.
- 47 measurement sets, all `in progress` on the IGVF portal as of 2026-05-07.
