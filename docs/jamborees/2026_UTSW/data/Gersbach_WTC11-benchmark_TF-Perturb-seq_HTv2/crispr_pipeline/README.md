# Gersbach HTv2 — CRISPR pipeline (bona fide reference)

This folder is the bona-fide reference CRISPR pipeline output for the 2026 UTSW jamboree. It is the reference against which the per-output schemas in `docs/jamborees/2026_UTSW/schemas/` are verified and the layout that production-dataset bundles should match.

## Run

| | |
|---|---|
| Dataset | `Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2` |
| Run label | `cleanser_800_mito_15pc` |
| Source | HPC: `aklie@nrnb-login.ucsd.edu:/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/runs/cleanser_800_mito_15pc/` |
| Guide assignment | `cleanser` |
| QC: `QC_min_genes_per_cell` | 800 |
| QC: `QC_pct_mito` | 15 |

Naming convention: `<guide-assignment-method>_<min_genes_per_cell>_mito_<pct_mito>pc`.

## What's in this folder

Canonical 3-folder pipeline output layout (`pipeline_dashboard/`, `pipeline_info/`, `pipeline_outputs/`).

| Path | Contents | Source | Source size |
|---|---|---|---|
| `pipeline_dashboard/` | `dashboard.html`, `inference_mudata.h5mu` (1.1 GB), `additional_qc/` (per-gene/guide/intended-target/trans QC tables), `evaluation_output/` (cis/trans bedgraph/bedpe), `figures/`, `guide_seqSpec_plots/`, `svg/`, `benchmark_output/` (TF benchmarking) | HPC `runs/cleanser_800_mito_15pc/pipeline_dashboard/` | ~1.0 GB |
| `pipeline_outputs/` | Final perturbo TSVs: `perturbo_cis_per_element_output.tsv.gz`, `perturbo_cis_per_guide_output.tsv.gz`, `perturbo_trans_per_element_output.tsv.gz`, `perturbo_trans_per_guide_output.tsv.gz` | HPC `runs/cleanser_800_mito_15pc/pipeline_outputs/` | ~150 MB |
| `pipeline_info/` | `nf_core_pipeline_software_versions.yml`, `params_2026-04-02_19-16-29.json` | GCS `gs://igvf-pertub-seq-pipeline-data/scratch/bioinfolucas/Benchmark_cleanser_800_mito_15pc/Gersbach_HTV2/pipeline_info/` | ~3 KB |

## What's not here

- `anndata/`, `calibration/`, `tf/` — additional artifacts that exist in the HPC source run but aren't part of the canonical pipeline output layout. We'll fold the downstream ones in later as needed.
- `pipeline_dashboard/evaluation_output/{cis_perturbo,cis_sceptre,trans_evaluation}.bedgraph` — empty (0-byte) in the source run, skipped during upload (Synapse rejects 0-byte files).

## Synapse

Mirrored 2026-05-09 to `2026_UTSW/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/crispr_pipeline/` — [`syn74885574`](https://www.synapse.org/Synapse:syn74885574).

## Schema

Full machine-readable schema: [`docs/jamborees/2026_UTSW/schemas/crispr_pipeline.json`](../../../schemas/crispr_pipeline.json).
Analysis-level walkthrough: [`docs/analysis/CRISPR_PIPELINE_OUTPUTS.md`](../../../../analysis/CRISPR_PIPELINE_OUTPUTS.md).
