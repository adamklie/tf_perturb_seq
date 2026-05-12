# Hon WTC11 Cardiomyocyte TF Perturb-seq (Hon CM)

Production TF Perturb-seq from the Hon lab (cardiomyocyte differentiation).
Synapse canonical: [`syn74520421`](https://www.synapse.org/Synapse:syn74520421) — Weizhou's CRISPR pipeline run, locally mirrored as `2026_04_19_no_spacer/` (QC + energy_distance re-runs).
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
2026_04_19_no_spacer/           # Weizhou's CRISPR pipeline run (= Synapse syn74520421)
├── qc/                         # our QC re-run on Weizhou's MuData (gitignored outputs)
└── energy_distance/            # Adam's ED run on Weizhou's MuData (configs TRACKED; results gitignored)
```

## Pipeline runs

| Local run | Source | Notes |
|---|---|---|
| `seqspec_v3` | GCS `2026_04_15/outs/seqspec_v3/` | New production run, May 2026 |
| `2026_04_19_no_spacer` | Weizhou's CRISPR pipeline run, Synapse [`syn74520421`](https://www.synapse.org/Synapse:syn74520421); MuData [`syn74522725`](https://www.synapse.org/Synapse:syn74522725) | **Canonical** — both ED and cNMF (Alexandra) run on its MuData |

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
