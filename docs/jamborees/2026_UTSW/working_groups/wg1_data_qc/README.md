# WG1 — Data summarization, pipeline, and QC

Topic 1 / Figure 1. The goal is to establish the primary building blocks for every downstream analysis: agree on guide detection and repression criteria, harmonize general statistics across the five production datasets, summarize transcriptome-wide significance from energy distance, and surface cross-lineage shared TFs and trans-target overlaps. End-of-jamboree deliverables include the Figure 1 panels and a written set of QC and harmonization recommendations for downstream WGs.

## Questions

- Which guides are detected and repress their targets in each system?
- How do general statistics (cell counts, gRNA and scRNA MOI, %mito) compare across datasets?
- How many TFs significantly alter the transcriptome per lineage, and how many are shared across lineages?
- Cluster TF perturbations by energy distance in each lineage; flag TFs whose effect differs sharply across lineages.
- What is the distribution of inferred trans targets per TF in each system, and how much do trans-target sets overlap across lineages?
- How should clonal expansion, doublet removal, and DEG calibration be folded into upstream analysis?

## Data

| Dataset | CRISPR pipeline | Energy distance |
|---|:---:|:---:|
| Hon WTC11 Cardiomyocyte | ready | ready |
| Huangfu HUES8 Definitive Endoderm | ready | ready |
| Huangfu HUES8 Embryonic Stem Cell | ready | ready |
| Gersbach WTC11 Hepatocyte | ready | - |
| Engreitz WTC11 Endothelial | - | - |

Per-dataset cards live under [`../../data/`](../../data/).

## Issues

- *[FILL IN issue link]*: ED calibration

## Working flow

Three steps, in order:

1. **Brainstorm** — sketch example figures, summary tables, and pseudocode that answer the questions above. Capture this in a notebook, doc, or notes in this folder.
2. **Execute** — run the analyses; commit code (notebooks, scripts, supporting docs) to this folder on GitHub.
3. **Share** — upload reusable outputs (figures, tables, intermediate data) to WG1's Synapse folder [`syn74954079`](https://www.synapse.org/Synapse:syn74954079) (mirrored `working_groups/wg1_data_qc/` under [`syn64423137/2026_UTSW/`](https://www.synapse.org/Synapse:syn64423137)). Record the syn ID for each upload in the [Outputs](#outputs) table below so the next person can find it.

## Outputs

| Object | Syn ID | Description | Owner |
|---|---|---|---|
| *[FILL IN as outputs land]* | | | |
