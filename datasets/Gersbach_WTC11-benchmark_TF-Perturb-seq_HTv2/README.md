# Gersbach WTC11 Benchmark TF Perturb-seq (HTv2)

WTC11 benchmark TF Perturb-seq from the Gersbach lab (HTv2 chemistry). IGVF analysis set `IGVFDS6237URFJ`. Synapse: [`syn73712262`](https://www.synapse.org/Synapse:syn73712262) (Alex's local run, 2026-02-17).

## Layout

```
setup/                    # Shared input generation
├── scripts/              # 1_–4_*.sh pipeline drivers + 5_run_energy_distance.sh
├── configs/              # Base Nextflow .config
└── samplesheets/         # Canonical sample_metadata.csv
<run>/                    # Per Lucas's parameter sweep (CRISPR) or analysis run
├── crispr_pipeline/
│   ├── pipeline_info/    # params_*.json + versions yml (TRACKED — small, run provenance)
│   └── pipeline_outputs/ / pipeline_dashboard/ / anndata/  (all gitignored)
├── calibration/          # FDR-controlled TSVs (analysis tier; result TSVs gitignored)
├── energy_distance/      # ED config*.json TRACKED; result CSVs/PDFs/logs gitignored
└── cnmf/<cnmf_run>/      # cNMF run (only on cleanser_800_mito_15pc so far)
    ├── Script/           # TRACKED — Stage 1/2/3 wrappers
    └── Data/ / Result/   # gitignored (bulk)
```

## Pipeline runs

| Local run | Lucas's GCS source | scrublet | method |
|---|---|---|---|
| `cleanser_500_mito_15pc` | `Benchmark_cleanser_500_mito_15pc` | off | cleanser |
| `cleanser_800_mito_15pc` | `Benchmark_cleanser_800_mito_15pc` | off | cleanser |
| `cleanser_extremes_200_mito_15pc` | `Benchmark_cleanser_extremes_200_sceptre_mito_15pc` | off | sceptre |
| `cleanser_extremes_2000_mito_15pc` | `Benchmark_cleanser_extremes_2000_sceptre_mito_15pc` | off | sceptre |
| `scrublet_off_cleanser_800_mito_15pc` | `Benchmark_cleanser_800_mito_15pc` | off | cleanser |
| `scrublet_on_sceptre_800_mito_15pc` | `REAL_SCRUBLETS_..._scrublet_sceptre_800_mito_15pc` | **on** | sceptre |

Lucas's GCS root: `gs://igvf-pertub-seq-pipeline-data/scratch/bioinfolucas/<lucas_dir>/Gersbach_HTV2/`.
Verify params from `<run>/crispr_pipeline/pipeline_info/params_*.json`.

Other runs (non-CRISPR-Nextflow):
- `2026_04_08_prelim/energy_distance/` — preliminary energy-distance run
- `cleanser_800_mito_15pc/cnmf/050926_HTv2_20iter_5KHVG_torch_halsvar_batch/` — cNMF (torch-halsvar batch, 20 iter, 5K HVG, k-skill A)

## Reproduce

```bash
DS=datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2
bash $DS/setup/scripts/1_generate_per_sample_metadata.sh
bash $DS/setup/scripts/2_upload_to_gcp.sh
bash $DS/setup/scripts/3_patch_gcp_files.sh
bash $DS/setup/scripts/4_run_CRISPR_pipeline.sh
bash $DS/setup/scripts/5_run_energy_distance.sh
```
