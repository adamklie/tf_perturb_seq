# Energy distance — Huangfu HUES8 Definitive Endoderm

| | |
|---|---|
| Dataset ID | `Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq` |
| Schema | [`../../../schemas/energy_distance.json`](../../../schemas/energy_distance.json) |
| Run label | `muddy_penguin` |
| Source MuData | `gs://igvf-pertub-seq-pipeline-data/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/2026_04_09/outs/muddy_penguin/inference_mudata.h5mu` (18.02 GB) |
| Status | ⏳ Queued on UCSD nrnb HPC 2026-05-08 (steps 1, 2, 2.1; step 3 deferred) |
| HPC output dir | `/cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/results/energy_distance/muddy_penguin/` |
| Synapse | _pending — will populate `synapse_paths.tsv` after mirror_ |

## How to run / re-run

From the HPC:

```bash
ssh aklie@nrnb-login.ucsd.edu
sbatch /cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/5_run_energy_distance.sh
```

Wraps the shared runner `/cellar/users/aklie/projects/tf_perturb_seq/scripts/run_energy_distance_pipeline.sh`. Container `edist_pipeline.sif` is shared across datasets; pulled on first run.

## Mirror outputs to Synapse (after run completes)

```bash
# from the HPC, in the project root
.venv/bin/python docs/jamborees/2026_UTSW/scripts/mirror_edistance_outputs.py \
  --dataset Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq \
  --source-dir /cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/results/energy_distance/muddy_penguin
```

## Notes

- Same source MuData as the Huangfu ESC e-distance run uses for ESC; both share the construct library (`IGVFDS3299AXST`) and guide file (`IGVFFI8270UPKB`).
- `muddy_penguin` is the canonical run because of its 13 bp `spacer_tag` (`GAGTACATGGGGG`); the earlier `entertaining_hamster` (12 bp) had ~6× lower guide UMI capture and is not used.
