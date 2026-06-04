# Hon WTC11 Benchmark TF Perturb-seq

WTC11 benchmark TF Perturb-seq from the Hon lab (Gary). IGVF analysis set `IGVFDS4761PYUO`. Synapse: [`syn73743227`](https://www.synapse.org/Synapse:syn73743227) (Lucas's GCP run).

## Layout

```
setup/                    # Shared input generation
├── scripts/              # 1_–6_*.sh pipeline drivers + 5_pipeline_comparison.ipynb
├── configs/              # Base Nextflow .config
└── samplesheets/         # Canonical sample_metadata.csv
seqspec/                  # Hon lab seqspec yamls (rna/guide/hash)
<run>/                    # Per Lucas's parameter sweep
├── crispr_pipeline/
│   ├── pipeline_info/    # params_*.json + versions yml (TRACKED — small)
│   └── pipeline_outputs/ / pipeline_dashboard/ / pipeline_qc/ / pipeline_comparisons/ / anndata/  (all gitignored)
└── calibration/          # FDR-controlled TSVs (analysis tier; result TSVs gitignored)
synapse_syn73743227/      # Synapse mirror provenance
└── cnmf/030726_20iter_5KHVG_torch_halsvar_batch_e7/  # cNMF run from Synapse-uploaded h5mu
elated_almeida/           # Older Adam-local Nextflow run (not in Lucas's parameter sweep)
```

## Canonical runs (CRISPR paper)

As of 2026-06-04, the canonical runs for Figure 1 are **one `cleanser` + one `sceptre` run**, differing *only* in `GUIDE_ASSIGNMENT_method` (all other params identical). Synced locally to `basic_threshold_{cleanser,sceptre}/crispr_pipeline/` — the three terminal folders `pipeline_dashboard/`, `pipeline_info/`, `pipeline_outputs/`. Full registry (params, sizes, lab/contact): [pipeline_runs.tsv](../../docs/manuscripts/CRISPRi_tech_benchmark/docs/pipeline_runs.tsv).

| Local run | GUIDE_ASSIGNMENT | GCS source (under `scratch/bioinfolucas/`) | Local status |
|---|---|---|---|
| `basic_threshold_sceptre` | sceptre | `Benchmark_basic_run_threshold/GaryHonDataset/` | ✓ synced |
| `basic_threshold_cleanser` | cleanser | `Benchmark_basic_run_threshold_cleanser/GaryHonDataset/` | ✓ synced |

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
| `elated_almeida` | older Adam-local Nextflow run (no Lucas mapping) | n/a | n/a |

Lucas's GCS root: `gs://igvf-pertub-seq-pipeline-data/scratch/bioinfolucas/<lucas_dir>/GaryHonDataset/`.
Verify params from `<run>/crispr_pipeline/pipeline_info/params_*.json`.

## Reproduce

```bash
DS=datasets/Hon_WTC11-benchmark_TF-Perturb-seq
bash $DS/setup/scripts/1_generate_per_sample_metadata.sh
bash $DS/setup/scripts/2_upload_to_gcp.sh
bash $DS/setup/scripts/3_run_CRISPR_pipeline.sh   # Hon's Nextflow runner (no patch step)
bash $DS/setup/scripts/6_run_energy_distance.sh
```
