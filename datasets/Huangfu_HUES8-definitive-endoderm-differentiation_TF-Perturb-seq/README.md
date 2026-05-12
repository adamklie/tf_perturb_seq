# Huangfu HUES8 Definitive Endoderm TF Perturb-seq (DE)

Production TF Perturb-seq from the Huangfu lab — HUES8 definitive endoderm differentiation. IGVF analysis set `IGVFDS5057HJKP`. Synapse mirror: [`syn74834952`](https://www.synapse.org/Synapse:syn74834952).

## Layout

```
setup/                          # Shared input generation
├── scripts/                    # 1_–5_*.sh pipeline drivers
├── configs/                    # Base Nextflow .config (muddy_penguin)
└── samplesheets/               # Canonical sample_metadata.csv
muddy_penguin/                  # Production CRISPR pipeline run
├── crispr_pipeline/
│   ├── pipeline_info/          # params_*.json + versions (TRACKED — small)
│   └── pipeline_outputs/ / pipeline_dashboard/ / anndata/  (all gitignored)
├── cnmf/042926_huangfu_de_torchcnmf_KskillA/   # cNMF Stage 1+ (Apr 2026, KskillA pattern)
│   ├── Script/                 # TRACKED — Stage 1/2/3 wrappers
│   └── Data/ / Result/         # gitignored (bulk)
├── qc/                         # QC outputs (gitignored)
└── energy_distance/            # ED outputs (configs TRACKED; results/logs/image gitignored)
```

## Pipeline runs

| Local run | Source | Notes |
|---|---|---|
| `muddy_penguin` | GCS `2026_04_09/outs/muddy_penguin/` | Production CRISPR + cNMF + ED + QC |

GCS canonical: `gs://igvf-pertub-seq-pipeline-data/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/2026_04_09/outs/muddy_penguin/`. An earlier `entertaining_hamster` run also exists on GCS but isn't mirrored locally.

cNMF run `042926_huangfu_de_torchcnmf_KskillA` is mirrored on Synapse at [`syn74893844`](https://www.synapse.org/Synapse:syn74893844).

## Reproduce

```bash
DS=datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq
bash $DS/setup/scripts/1_generate_per_sample_metadata.sh
bash $DS/setup/scripts/2_upload_to_gcp.sh
bash $DS/setup/scripts/3_patch_gcp_files.sh
bash $DS/setup/scripts/4_run_CRISPR_pipeline.sh
bash $DS/setup/scripts/5_run_energy_distance.sh
```
