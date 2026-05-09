# Energy distance — Huangfu HUES8 Embryonic Stem Cell

| | |
|---|---|
| Dataset ID | `Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq` |
| Schema | [`../../../schemas/energy_distance.json`](../../../schemas/energy_distance.json) |
| Run label | `sceptre_v1` |
| Source MuData | `gs://igvf-pertub-seq-pipeline-data/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/2026_04_13/outs/sceptre_v1/inference_mudata.h5mu` (17.52 GB) |
| Status | ✅ Complete (2026-05-09; SLURM job 10547705, 10h 0m on carter-gpu-02) |
| HPC output dir | `/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/results/energy_distance/sceptre_v1/` |
| Synapse | [`syn74883475`](https://www.synapse.org/Synapse:syn74883475) |
| Validation | All 4 layers PASS (file presence, schema, value ranges, schema-identity vs HTv2 reference). 2267 targets × 46 cols. ⚠ Calibration concern, see below. |

## How to run / re-run

From the HPC:

```bash
ssh aklie@nrnb-login.ucsd.edu
sbatch /cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/5_run_energy_distance.sh
```

Wraps the shared runner `/cellar/users/aklie/projects/tf_perturb_seq/scripts/run_energy_distance_pipeline.sh`.

## Mirror outputs to Synapse (after run completes)

```bash
# from the HPC, in the project root
.venv/bin/python docs/jamborees/2026_UTSW/scripts/mirror_edistance_outputs.py \
  --dataset Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq \
  --source-dir /cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/results/energy_distance/sceptre_v1
```

## Notes

- Companion run to the Huangfu DE e-distance; the two share construct library (`IGVFDS3299AXST`) and guide file (`IGVFFI8270UPKB`), differing only in the differentiation state.
- ESC distances are ~5.7× larger than DE (median 527 vs 93) — undifferentiated stem cells appear to have stronger TF-perturbation effects (or higher baseline cell-state variance — see calibration concern).

## ⚠ Calibration concern (2026-05-09)

**Same calibration issue as the DE sibling — see [DE README](../../Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/energy_distance/README.md#calibration-concern-2026-05-09) for full diagnosis.**

ESC-specific numbers:
- All 100 negative-control targets have `pval_mean = 0`
- NC distance_mean median = 528.3 vs targeting distance_mean median = 527.3 (ratio 1.00× — identical)
- HTv2 reference distance_mean median: 1.78
- ESC's distances are even larger than DE's (median 527 vs 93) — same root causes apply, amplified by stem cells' more variable transcriptional state.

P-values not reliable; use raw distance as effect-size proxy instead (e.g., distance > NC max ≈ 602 for ESC).
