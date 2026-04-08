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
- `scripts/generate_per_sample.py` — Query portal API for metadata
- `scripts/upload_to_gcp.py` — Upload files to GCS
- `scripts/validate_gcp_paths.py` — Verify uploads

## Stage 2: CRISPR Pipeline

Run the [IGVF CRISPR Pipeline](https://github.com/IGVF/CRISPR_Pipeline) (Nextflow) on GCP.

**Prerequisites:**
- GCP access to `igvf-pertub-seq-pipeline` project
- Service account credentials
- Pipeline config file (per-dataset, stored in `datasets/<name>/`)
- Guide metadata file (from `ref/finalized_annotation_files/`)
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

See [docs/analysis/CRISPR_PIPELINE.md](analysis/CRISPR_PIPELINE.md) for detailed instructions.

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

## Stage 5: Gene Program Discovery

Unsupervised discovery of gene expression programs using cNMF.

**External tools:** `external/cNMF_benchmarking/` (git submodule)

**Per-dataset:** `bash datasets/<DATASET>/6_run_cnmf.sh`

**Key steps:**
1. **k-selection** — Run cNMF across range of k values (e.g., 50-300), evaluate by:
   - Unique perturbations recovered (FDR cutoff based on safe-targeting genes)
   - Unique pathways recovered
2. **Program inference** — Run cNMF at selected k with many iterations (700k-1M)
3. **GSEA** — Pathway enrichment for each program
4. **TF-program assignment** — Map perturbations to programs

**Runtime:** ~1 business week for k-selection + final run (parallelized, GPU-accelerated)

See [docs/analysis/cNMF.md](analysis/cNMF.md) for details.

## Standard preprocessing for integrative analysis

After pipeline output, the standard preprocessing for cross-dataset analysis:
1. Depth normalization + log1p
2. Select top 2k-3k highly variable genes
3. PCA on variable genes
4. Clustering (Leiden/Louvain)
5. UMAP for 2D visualization
