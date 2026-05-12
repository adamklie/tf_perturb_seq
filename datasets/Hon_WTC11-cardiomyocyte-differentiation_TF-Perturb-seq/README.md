# Hon WTC11 Cardiomyocyte TF Perturb-seq (Hon CM)

Production TF Perturb-seq from the Hon lab (cardiomyocyte differentiation).
Synapse canonical: [`syn73582673`](https://www.synapse.org/Synapse:syn73582673) — Weizhou's local CRISPR pipeline run, locally mirrored as `weizhou_syn74520421/` (QC re-run only).
Production GCS run: `gs://igvf-pertub-seq-pipeline-data/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/2026_04_15/outs/seqspec_v3/` — mirrored locally as `seqspec_v3/`.

## Layout

```
setup/                          # Shared input generation (Adam's code)
├── scripts/                    # 1_–5_*.sh + generate_per_sample.py + make_patched_metadata.py
├── configs/                    # Base Nextflow .config (seqspec_v3 — the new prod run)
└── samplesheets/               # Canonical sample_metadata.csv
seqspec/                        # Hon lab seqspec yamls (rna/guide/hash) — legacy
synapse_inference_mudata/       # Sara/Weizhou mudata download (gitignored)
HonLabInternal/                 # Hon-internal cNMF results (gitignored entirely)
seqspec_v3/                     # Production CRISPR pipeline run (GCS 2026_04_15)
├── crispr_pipeline/
│   ├── pipeline_info/          # params_*.json + versions (TRACKED — small)
│   └── pipeline_outputs/ / pipeline_dashboard/ / anndata/ (all gitignored)
├── qc/                         # QC re-run locally (results gitignored)
└── cnmf/051126_honcm_torchcnmf_KskillA/   # cNMF Stage 1 (May 2026, KskillA pattern)
    ├── Script/                 # TRACKED — Stage 1 SLURM wrappers
    └── Data/ / Result/         # gitignored (bulk)
weizhou_syn74520421/            # Weizhou's Synapse-imported h5mu run (QC-only on our side)
└── qc/                         # QC outputs (gitignored)
2026_04_19_no_spacer/           # Older energy-distance run (no associated CRISPR pipeline locally)
└── energy_distance/            # ED outputs (config JSONs TRACKED; results gitignored)
```

## Pipeline runs

| Local run | Source | Notes |
|---|---|---|
| `seqspec_v3` | GCS `2026_04_15/outs/seqspec_v3/` | New production run, May 2026 |
| `weizhou_syn74520421` | Synapse [`syn74520421`](https://www.synapse.org/Synapse:syn74520421) → Weizhou's [`syn73582673`](https://www.synapse.org/Synapse:syn73582673) | Canonical h5mu; QC-only locally |
| `2026_04_19_no_spacer` | Older ED run | No CRISPR pipeline (no_spacer modality) |

cNMF Stage 1 (`051126_honcm_torchcnmf_KskillA`) was kicked off against the new `seqspec_v3` h5mu — see issue [#20](https://github.com/adamklie/tf_perturb_seq/issues/20).

## Reproduce

```bash
DS=datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq
bash $DS/setup/scripts/1_generate_per_sample_metadata.sh
bash $DS/setup/scripts/2_upload_to_gcp.sh
bash $DS/setup/scripts/3_patch_gcp_files.sh
bash $DS/setup/scripts/4_run_CRISPR_pipeline.sh
bash $DS/setup/scripts/5_run_energy_distance.sh
```
