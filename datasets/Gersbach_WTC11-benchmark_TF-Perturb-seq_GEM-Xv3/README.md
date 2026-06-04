# Gersbach WTC11 Benchmark TF Perturb-seq (GEM-X v3)

WTC11 benchmark TF Perturb-seq from the Gersbach lab (GEM-X v3 chemistry). IGVF analysis set `IGVFDS6673ZFFG`. Synapse: [`syn73712373`](https://www.synapse.org/Synapse:syn73712373) (Alex's local run `jade_dolphin`, downloaded 2026-02-17); also `frisky_rabbit` ([`syn73579845`](https://www.synapse.org/Synapse:syn73579845), sceptre).

## Layout

```
setup/                    # Shared input generation
├── scripts/              # 1_–4_*.sh pipeline drivers
│   └── cnmf/             # 6_run_cnmf.sh (gene-program discovery template)
├── configs/              # Base Nextflow .config
└── samplesheets/         # Step 1 → 2 → 3 lineage + IGVFDS6673ZFFG.updated.v4.csv (Alex's curated)
<run>/                    # One per Lucas's parameter sweep
├── crispr_pipeline/
│   ├── pipeline_info/    # params_*.json + versions yml (TRACKED — small, run provenance)
│   └── pipeline_outputs/ / pipeline_dashboard/ / anndata/  (all gitignored)
└── calibration/          # FDR-controlled TSVs (analysis tier; result TSVs gitignored)
```

## Canonical runs (CRISPR paper)

As of 2026-06-04, the canonical runs for Figure 1 are **one `cleanser` + one `sceptre` run**, differing *only* in `GUIDE_ASSIGNMENT_method` (all other params identical). Synced locally to `basic_threshold_{cleanser,sceptre}/crispr_pipeline/` — the three terminal folders `pipeline_dashboard/`, `pipeline_info/`, `pipeline_outputs/`. Full registry (params, sizes, lab/contact): [pipeline_runs.tsv](../../docs/manuscripts/CRISPRi_tech_benchmark/docs/pipeline_runs.tsv).

| Local run | GUIDE_ASSIGNMENT | GCS source (under `scratch/bioinfolucas/`) | Local status |
|---|---|---|---|
| `basic_threshold_sceptre` | sceptre | `Benchmark_basic_run_threshold/Gersbach_gemX_v3/` | ✓ synced |
| `basic_threshold_cleanser` | cleanser | `Benchmark_basic_run_threshold_cleanser/Gersbach_gemX_v3/` | ✓ synced |

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

Lucas's GCS root: `gs://igvf-pertub-seq-pipeline-data/scratch/bioinfolucas/<lucas_dir>/Gersbach_gemX_v3/`.
Verify params from `<run>/crispr_pipeline/pipeline_info/params_*.json`.

## Reproduce

```bash
DS=datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3
bash $DS/setup/scripts/1_generate_per_sample_metadata.sh  # → setup/samplesheets/sample_metadata.csv
bash $DS/setup/scripts/2_upload_to_gcp.sh                 # → setup/samplesheets/sample_metadata_gcp_<date>.csv
bash $DS/setup/scripts/3_patch_gcp_files.sh               # → setup/samplesheets/sample_metadata_gcp_<date>_patched.csv
bash $DS/setup/scripts/4_run_CRISPR_pipeline.sh           # runs Nextflow on GCP
bash $DS/setup/scripts/cnmf/6_run_cnmf.sh                 # cNMF gene-program discovery (no run yet)
```
