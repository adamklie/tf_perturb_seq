# Huangfu WTC11 Benchmark TF Perturb-seq

WTC11 benchmark TF Perturb-seq from the Huangfu lab (Gary, IGVF). IGVF analysis set `IGVFDS5057HJKP`. Synapse: [`syn72386406`](https://www.synapse.org/Synapse:syn72386406) (Stanford-side run, uploaded by Gary).

## Layout

```
setup/                    # Shared input generation
├── scripts/              # 1_–4_*.sh pipeline drivers
├── configs/              # Base Nextflow .config
└── samplesheets/         # Step 1 → 2 → 3 lineage:
                          #   sample_metadata.csv → ..._gcp_2026_01_30.csv → ..._gcp_2026_01_30_patched.csv (+ _v2)
<run>/                    # Per Lucas's parameter sweep
├── crispr_pipeline/
│   ├── pipeline_info/    # params_*.json + versions (TRACKED — small)
│   └── pipeline_outputs/ / pipeline_dashboard/ / anndata/  (all gitignored)
└── calibration/          # FDR-controlled TSVs (analysis tier; result TSVs gitignored)
gary_syn72386406/         # Synapse provenance for the Stanford-side Yan Mo h5mu
└── cnmf/                 # cNMF run on /oak/stanford/.../IGVF_Huangfu_WTC11/Data/inference_mudata.h5ad
    ├── Script/           # TRACKED
    └── RUN_NAME.txt      # 030726_20iter_5KHVG_torch_halsvar_batch_e7
    └── Data/ / Result/   # gitignored (bulk)
```

## Canonical runs (CRISPR paper)

As of 2026-06-04, the canonical runs for Figure 1 are **one `cleanser` + one `sceptre` run**, differing *only* in `GUIDE_ASSIGNMENT_method` (all other params identical). Synced locally to `basic_threshold_{cleanser,sceptre}/crispr_pipeline/` — the three terminal folders `pipeline_dashboard/`, `pipeline_info/`, `pipeline_outputs/`. Full registry (params, sizes, lab/contact): [pipeline_runs.tsv](../../docs/manuscripts/CRISPRi_tech_benchmark/docs/pipeline_runs.tsv).

| Local run | GUIDE_ASSIGNMENT | GCS source (under `scratch/bioinfolucas/`) | Local status |
|---|---|---|---|
| `basic_threshold_sceptre` | sceptre | `Benchmark_basic_run_threshold/HuangFuDataset/` | ✓ synced |
| `basic_threshold_cleanser` | cleanser | `Benchmark_basic_run_threshold_cleanser/HuangFuDataset/` | ✓ synced |

These supersede the 7-run parameter sweep below (retained but no longer canonical).

## Pipeline runs

| Local run | Lucas's GCS source | scrublet | method |
|---|---|---|---|
| `cleanser_500_mito_15pc` | `Benchmark_cleanser_500_mito_15pc` | off | cleanser |
| `cleanser_800_mito_15pc` | `Benchmark_cleanser_800_mito_15pc` | off | cleanser |
| `cleanser_extremes_200_mito_15pc` | `Benchmark_cleanser_extremes_200_sceptre_mito_15pc` | off | sceptre |
| `cleanser_extremes_2000_mito_15pc` | `Benchmark_cleanser_extremes_2000_sceptre_mito_15pc` | off | sceptre |
| `cleanser_knee2_mito_15pc` | `Benchmark_cleanser_knee2_mito_15pc` | off | cleanser |
| `scrublet_off_cleanser_800_mito_15pc` | `Benchmark_cleanser_800_mito_15pc` | off | cleanser |
| `scrublet_on_sceptre_800_mito_15pc` | `REAL_SCRUBLETS_..._scrublet_sceptre_800_mito_15pc` | **on** | sceptre |

Lucas's GCS root: `gs://igvf-pertub-seq-pipeline-data/scratch/bioinfolucas/<lucas_dir>/HuangFuDataset/`.

## Reproduce

```bash
DS=datasets/Huangfu_WTC11-benchmark_TF-Perturb-seq
bash $DS/setup/scripts/1_generate_per_sample_metadata.sh
bash $DS/setup/scripts/2_upload_to_gcp.sh
bash $DS/setup/scripts/3_patch_gcp_files.sh
bash $DS/setup/scripts/4_run_CRISPR_pipeline.sh
```
