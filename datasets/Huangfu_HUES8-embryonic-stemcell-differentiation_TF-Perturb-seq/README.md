# Huangfu HUES8 Embryonic Stem Cell TF Perturb-seq (ESC)

Production TF Perturb-seq from the Huangfu lab — HUES8 embryonic stem cells. Synapse mirror: [`syn74835010`](https://www.synapse.org/Synapse:syn74835010).

## Layout

```
setup/                          # Shared input generation
├── scripts/                    # 1_–5_*.sh pipeline drivers
├── configs/                    # Base Nextflow .config (sceptre_v1)
└── samplesheets/               # Canonical sample_metadata.csv
sceptre_v1/                     # Production CRISPR pipeline run
├── crispr_pipeline/
│   ├── pipeline_info/          # params_*.json + versions (TRACKED — small)
│   └── pipeline_outputs/ / pipeline_dashboard/ / anndata/  (all gitignored)
├── cnmf/042926_huangfu_esc_torchcnmf_KskillA/   # cNMF Stage 1+ (Apr 2026, KskillA pattern)
│   ├── Script/                 # TRACKED — Stage 1/2/3 wrappers
│   └── Data/ / Result/         # gitignored (bulk)
├── qc/                         # QC outputs (gitignored)
└── energy_distance/            # ED outputs (configs TRACKED; logs/image gitignored)
```

## Pipeline runs

| Local run | Source | Notes |
|---|---|---|
| `sceptre_v1` | GCS `2026_04_13/outs/sceptre_v1/` | Production CRISPR + cNMF + ED + QC |

GCS canonical: `gs://igvf-pertub-seq-pipeline-data/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/2026_04_13/outs/sceptre_v1/`.

cNMF run `042926_huangfu_esc_torchcnmf_KskillA` is mirrored on Synapse at [`syn74893846`](https://www.synapse.org/Synapse:syn74893846).

## Reproduce

```bash
DS=datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq
bash $DS/setup/scripts/1_generate_per_sample_metadata.sh
bash $DS/setup/scripts/2_upload_to_gcp.sh
bash $DS/setup/scripts/3_patch_gcp_files.sh
bash $DS/setup/scripts/4_run_CRISPR_pipeline.sh
bash $DS/setup/scripts/5_run_energy_distance.sh
```
