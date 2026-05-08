# Energy distance — Huangfu HUES8 Embryonic Stem Cell

| | |
|---|---|
| Dataset ID | `Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq` |
| Schema | [`../../../schemas/energy_distance.json`](../../../schemas/energy_distance.json) |
| Run label | `sceptre_v1` |
| Source MuData | `gs://igvf-pertub-seq-pipeline-data/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/2026_04_13/outs/sceptre_v1/inference_mudata.h5mu` (17.52 GB) |
| Status | ⏳ Queued on UCSD nrnb HPC 2026-05-08 (steps 1, 2, 2.1; step 3 deferred) |
| HPC output dir | `/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/results/energy_distance/sceptre_v1/` |
| Synapse | _pending — will populate `synapse_paths.tsv` after mirror_ |

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

- Companion run to the Huangfu DE e-distance; the two share construct library (`IGVFDS3299AXST`) and guide file (`IGVFFI8270UPKB`), differing only in the differentiation state. Worth a paired comparison once both runs complete.
