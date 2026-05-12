# Pipelines

## Overview

The full TF Perturb-seq workflow has 5 stages. Each dataset directory contains numbered shell scripts corresponding to these stages.

```
Raw fastqs ─→ [1] IGVF Portal ─→ [2] CRISPR Pipeline ─→ [3] QC ─→ [4] Energy Distance ─→ [5] Gene Programs
                  (upload)          (GCP/Nextflow)        (local)     (local/SLURM)         (local/SLURM)
```

## Stage 1: Upload to IGVF Portal

Data producers upload raw fastqs and metadata to the [IGVF data portal](https://data.igvf.org/).

**Expected portal structure per dataset:**
- One or more measurement sets for scRNA-seq (split by cellular sub pools)
- One auxiliary set per measurement set for gRNA-seq
- Optionally, a second auxiliary set per measurement set for HTO-seq

**Scripts:**
- `src/tf_perturb_seq/portal/generate_per_sample.py` — Query portal API for metadata
- `src/tf_perturb_seq/gcp/upload_to_gcp.py` — Upload files to GCS
- `src/tf_perturb_seq/gcp/validate_gcp_paths.py` — Verify uploads

## Stage 2: CRISPR Pipeline

Run the [IGVF CRISPR Pipeline](https://github.com/IGVF/CRISPR_Pipeline) (Nextflow) on GCP.

**Prerequisites:**
- GCP access to `igvf-pertub-seq-pipeline` project
- Service account credentials
- Pipeline config file (per-dataset, stored in `datasets/<name>/`)
- Guide metadata file (from `ref/guide_libraries/harmonized/`)
- Seqspec file (technology-specific)

**Per-dataset scripts:**
```bash
# Generate sample metadata from portal
bash datasets/<DATASET>/1_generate_per_sample_metadata.sh

# Upload fastqs to GCP (dry run first)
DRY_RUN=true bash datasets/<DATASET>/2_upload_to_gcp.sh
bash datasets/<DATASET>/2_upload_to_gcp.sh

# Patch gzipped files if needed
bash datasets/<DATASET>/3_patch_gcp_files.sh

# Run pipeline
RUN_IN_BACKGROUND=true bash datasets/<DATASET>/4_run_CRISPR_pipeline.sh
```

**Primary output:** `inference_mudata.h5mu` — MuData object with gene expression, guide assignments, and inference results.

**Inference methods:**
- `sceptre` — Current default for DE calling
- `cleanser` — Alternative assignment/inference
- `perturbo` — Being evaluated (outputs confidence estimates)

See [`crispr_pipeline/CRISPR_PIPELINE.md`](crispr_pipeline/CRISPR_PIPELINE.md) for detailed instructions, and [`crispr_pipeline/CRISPR_PIPELINE_OUTPUTS.md`](crispr_pipeline/CRISPR_PIPELINE_OUTPUTS.md) for the output directory reference.

## Stage 3: QC Pipeline

Quality control on CRISPR pipeline output MuData.

**Script:** `scripts/run_qc_pipeline.sh`

**Per-dataset:** `bash datasets/<DATASET>/4_run_qc_pipeline.sh [--dry-run]`

**QC modules** (in `src/tf_perturb_seq/qc/`):
- `mapping_gene.py` — Gene expression metrics (UMIs, genes detected, mito%)
- `mapping_guide.py` — Guide capture metrics (UMIs, guides/cell, cells/guide)
- `intended_target.py` — Knockdown efficiency (log2FC, auROC)
- Trans effects — DEGs per guide

**Key thresholds** (from team consensus):
| Metric | Threshold |
|--------|-----------|
| Median gene UMI/cell | > 3000 |
| Median gRNA UMI/cell | > 1000 |
| % cells with sgRNA | > 80% |
| Cells per target (post-filter) | ~750 |
| Knockdown threshold | 60% (Hon lab standard) |

## Stage 4: Energy Distance Analysis

Quantifies perturbation effects using energy distance metrics. Identifies significant perturbations and clusters them by phenotype.

**External tool:** [energy_dist_pipeline](https://github.com/Chikara-Takeuchi/energy_dist_pipeline) (git submodule in `external/`)

**Script:** `scripts/run_energy_distance_pipeline.sh`

**Per-dataset:** `bash datasets/<DATASET>/5_run_energy_distance.sh [--use-harmony] [--dry-run]`

**Pipeline steps:**
1. **MuData Adapter** — Convert MuData to pipeline inputs (h5ad + pickle + csv)
   - Normalize, log1p, HVG selection, regress confounders, scale
   - PCA: 50 components (configurable)
2. **gRNA Filtering** — DISCO test for outlier guides, k-means for non-targeting QC
3. **Energy Distance** — Permutation tests (each target vs. non-targeting background)
4. **Visualization** — Diagnostic plots
5. **Phenotype Clustering** — Affinity propagation clustering of significant perturbations

See [`energy_dist/ENERGY_DISTANCE.md`](energy_dist/ENERGY_DISTANCE.md) and [`energy_dist/ENERGY_DISTANCE_OUTPUTS.md`](energy_dist/ENERGY_DISTANCE_OUTPUTS.md).

## Stage 5: Gene Program Discovery

Unsupervised discovery of gene expression programs using **consensus non-negative matrix factorization (cNMF)**, wrapped by the [PerturbNMF](https://github.com/EngreitzLab/PerturbNMF) tool (supersedes the older `cNMF_benchmarking`).

**External tool:** `external/PerturbNMF/` (git submodule)

**Per-dataset run directory:** `datasets/<dataset>/<run>/cnmf/<run_name>/` — per-stage SLURM scripts in `Script/`, outputs in `Result/<run_name>/`.

**Pipeline stages:**
1. **Stage 1 — Inference** (GPU torch-cNMF): 8-K sweep at `dt = 2.0`, ~10–20 iterations per K
2. **Stage 2a — Evaluation**: perturbation association, GO / geneset / trait enrichment, explained variance (per K)
3. **Stage 2b — U-test calibration**: fake-targeting null distribution for FDR calibration (per K)
4. **Stage 3a — K-selection panel**: stability + EV + enrichment + perturbation-recovery vs K (group decision)
5. **Stage 3c — Per-target PDFs**: one volcano + UMAP + per-program effect plot per perturbed TF
6. **Stage 3e — Excel summary**: program-level + target-level + full perturbation-association tables

Stage 3b (per-program PDFs) is currently deferred — see upstream [issue #7](https://github.com/EngreitzLab/PerturbNMF/issues/7).

**Runtime:** Stage 1 takes 3–5 h GPU; full Stage 2/3 sequence runs in ~12–16 h CPU per dataset (see PerturbNMF.md compute-budget table).

See [`cnmf/PerturbNMF.md`](cnmf/PerturbNMF.md) for the how-to-run conventions and [`cnmf/cNMF_OUTPUTS.md`](cnmf/cNMF_OUTPUTS.md) for the output directory reference.

## Standard preprocessing for integrative analysis

After pipeline output, the standard preprocessing for cross-dataset analysis:
1. Depth normalization + log1p
2. Select top 2k-3k highly variable genes
3. PCA on variable genes
4. Clustering (Leiden/Louvain)
5. UMAP for 2D visualization
