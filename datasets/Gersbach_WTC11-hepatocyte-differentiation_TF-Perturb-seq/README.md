# Gersbach WTC11 Hepatocyte TF Perturb-seq (Hep)

WTC11 hepatocyte-differentiation TF Perturb-seq from the Gersbach lab. Synapse: [`syn70518849`](https://www.synapse.org/Synapse:syn70518849); Sara's processed h5mu mirrored locally as `sara_synapse_syn74842722` (downloaded from [`syn74842722`](https://www.synapse.org/Synapse:syn74842722) on 2026-05-11). See issue [#28](https://github.com/adamklie/tf_perturb_seq/issues/28) for open gaps.

## Layout

```
setup/                              # Shared input generation (Adam's code)
├── scripts/                        # 1_–4_*.sh + create_specfile.R, create_inference_volcano.py
├── configs/                        # Base Nextflow .config (cleanser_initial)
└── samplesheets/                   # sample_metadata.csv + hep_measurement_sets.txt + hep_ms_aux_pairs.tsv + generate.log + seqspec_BL138_S24_L008.yaml
sara_synapse_syn74842722/           # Run from Sara's Synapse-imported h5mu
├── crispr_pipeline/                # (empty — pipeline wasn't re-run; she ran it on her side)
└── qc/                             # QC outputs from local re-run (initial_qc + mapping_gene + mapping_guide + intended_target)
                                    # PNGs/TSVs/PDFs gitignored; on-disk only
synapse_inference_mudata/           # Sara's h5mu (gitignored)
```

## Pipeline runs

Only one run on disk:

| Local run | Source | Notes |
|---|---|---|
| `sara_synapse_syn74842722` | Sara's h5mu from [`syn74842722`](https://www.synapse.org/Synapse:syn74842722) | QC re-run locally; CRISPR pipeline wasn't re-run on our side |

No GCS canonical (the CRISPR pipeline run was done by Sara, not on the IGVF GCS).

## Reproduce

```bash
DS=datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq

# If/when we re-run the CRISPR pipeline ourselves:
bash $DS/setup/scripts/1_generate_per_sample_metadata.sh
bash $DS/setup/scripts/2_upload_to_gcp.sh
bash $DS/setup/scripts/3_patch_gcp_files.sh
bash $DS/setup/scripts/4_run_CRISPR_pipeline.sh

# Or use Sara's h5mu directly (current path):
ls $DS/synapse_inference_mudata/        # Sara's mudata
bash scripts/qc_array.sh                # local QC re-run
```
