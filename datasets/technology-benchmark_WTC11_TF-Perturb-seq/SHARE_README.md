# Technology Benchmark Analysis — Sharing Guide

## What this is

Analysis code for comparing 5 Perturb-seq technologies on the same WTC11 iPSC benchmark dataset. This was built for the TFP3 pillar project's cross-technology comparison and is fairly specific to our 5-dataset benchmark setup. That said, the individual notebooks and source modules may be useful as starting points for similar analyses on other datasets.

## What's most reusable

- **QC modules** (`src/tf_perturb_seq/qc/`) — standalone CLI scripts for computing gene expression, guide capture, and intended target metrics from any `inference_mudata.h5mu`. These are general-purpose.
- **Inference calibration** (`src/tf_perturb_seq/inference/calibrate.py`) — empirical null calibration of Perturbo results using non-targeting controls. Works on any dataset with NTC elements.
- **Color/config system** (`config/`) — YAML-based color palettes loaded via `config/loader.py`. Easy to adapt.

## What's benchmark-specific

- The **comparison notebooks** (`bin/2_qc/`, `bin/3_inference/`, `bin/0_cross_run_comparison.ipynb`) are hardcoded to read from manifest TSVs that point to our 5 benchmark datasets. They'd need modification for different datasets.
- The **data sync** (`bin/1_get_data/`) assumes our specific GCP bucket structure.
- The **manifests** (`manifests/`) are our dataset paths.

## Directory structure

```
bin/
  0_cross_run_comparison.ipynb   # CLEANSER vs Sceptre comparison
  run_all_notebooks.sh           # Papermill runner for all notebooks
  generate_report.py             # HTML report generation
  1_get_data/                    # GCP sync and manifest generation
  2_qc/                          # Cross-tech QC notebooks (3 notebooks + detailed README)
  3_inference/                   # Calibration + DEG + TF enrichment notebooks
  4_integrate/                   # scVI/scANVI integration
  5_upload/                      # Synapse upload
manifests/                       # TSV manifests pointing to dataset paths
results/cross_tech_comparison/   # All outputs land here
config/colors/                   # Color palettes (YAML)
src/tf_perturb_seq/
  qc/                            # Standalone QC modules (mapping_gene.py, mapping_guide.py, intended_target.py)
  inference/                     # Calibration module (calibrate.py)
```

## How to run

```bash
# 1. Sync data from GCP (edit paths inside first)
bash bin/1_get_data/1_sync_benchmark.sh

# 2. Run all QC + inference notebooks for a pipeline run
bash bin/run_all_notebooks.sh cleanser_unified

# 3. Run cross-run comparison (after running both cleanser_unified and sceptre_v11)
jupyter notebook bin/0_cross_run_comparison.ipynb
```

Each notebook is parameterized with `RUN_LABEL` via papermill. The runner script injects this automatically.

## Key files to read first

1. `bin/2_qc/README.md` — detailed description of every plot and metric
2. `src/tf_perturb_seq/qc/plan.md` — QC module schemas and design
3. `src/tf_perturb_seq/inference/plan.md` — calibration procedure and output schemas
4. `config/colors/technology-benchmark_WTC11_TF-Perturb-seq.yaml` — color scheme

## Environment

Python 3.10, managed with `uv`. Key dependencies: scanpy, anndata, mudata, pandas, numpy, matplotlib, seaborn, scipy. Install with `uv sync` from the project root.
