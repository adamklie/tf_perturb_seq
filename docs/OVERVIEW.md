# TF Perturb-seq Pillar Project (TFP3)

## What is this project?

A multi-lab IGVF consortium effort to perform CRISPRi Perturb-seq targeting ~2000 transcription factors (TFs) and epigenetic modifiers across multiple human cell lineages. The goal is to build an integrated picture of how TFs regulate gene programs, how those programs vary across cell types, and how they connect to human disease.

## Why does it matter?

- **Gene programs**: Define which TF perturbations activate or repress specific transcriptional programs in each cell type
- **Cross-lineage comparison**: Understand which regulatory relationships are conserved vs. lineage-specific
- **Disease connection**: Link TF perturbation effects to GWAS variants and disease-relevant cell types
- **Prediction**: Build models that predict perturbation effects in unseen cell types using sequence models, chromatin, and GRN features

## Paper vision

The project aims to produce an integrative analysis paper organized around 5-6 figures:

| Figure | Focus | Key Analysis |
|--------|-------|-------------|
| **Fig 1** | Overview + descriptive | Experimental design, processing pipeline, # significant TFs per lineage |
| **Fig 2** | Gene expression programs | cNMF program discovery, pathway enrichment, perturbation sensitivity, GWAS enrichment |
| **Fig 3** | Multiome integration + GRNs | ChromBPNet, E2G linking, TF importance comparison (multiome vs. Perturb-seq) |
| **Fig 4** | Deep dive on specific TF(s) | Understudied TF with high effect, pathway map, mechanism proposal, GWAS SNPs |
| **Fig 5** | TF regulatory grammar (ML) | Supervised model predicting perturbation effects from sequence/chromatin/GRN features |
| **Fig 6** | Genomic variation (TBD) | Pleiotropy interpretation |

See the full paper outline: [Google Doc](https://docs.google.com/document/d/1oPDLB3SvD8pKqzwJdIFgQTcDlEK-P_NaC7EJxS7Fsl4/edit) or local copy: `TF Perturb-seq Pillar Project - Integrative Analysis Paper Outline.md`

## Three analysis sub-aims

1. **Basic integrative analysis** (Adam, Sid, Chikara) — Energy distance, DE networks, cross-lineage comparison
2. **Predicting perturbation effects** (Weizhou) — PerturbNet, counterfactual predictions
3. **Integrating sequence models** (Revant, Gary) — ChromBPNet, E2G, cis-GRNs, variant interpretation

## Repository structure

```
tf_perturb_seq/
├── datasets/           # One directory per experiment (benchmark + production)
├── docs/               # Project documentation (you are here)
│   ├── analysis/       # Pipeline-specific docs (CRISPR_PIPELINE.md, cNMF.md)
│   ├── dataset_template/  # Template for setting up new datasets
│   └── jamborees/      # Annual meeting materials
├── external/           # Git submodules (energy_dist_pipeline, cNMF tools)
├── ref/                # Reference files (guide metadata, gene annotations, gene sets)
├── scripts/            # Shared pipeline scripts (QC, energy distance, upload, etc.)
├── src/tf_perturb_seq/ # Python package (QC modules, inference, utilities)
├── tests/              # Test suite
├── scratch/            # Exploratory work (gitignored, not committed)
├── config/             # Color schemes and shared configs
└── pyproject.toml      # Python project metadata and dependencies (managed with uv)
```
