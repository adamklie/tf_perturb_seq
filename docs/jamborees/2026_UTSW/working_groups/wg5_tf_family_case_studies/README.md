# WG5 — TF family case studies

Topic 2.4 / Figure 4. The goal is a deep dive into specific TF families (e.g., C2H2 zinc fingers, homeodomains) with newly implicated roles in lineage differentiation: which families show coordinated activity, what biology their regulated genes converge on, and which family members are worth a focused case study. End-of-jamboree deliverables include a family activity scorecard, per-family pathway-enrichment summaries, and one or two written case studies for Figure 4.

## Questions

- Which TF families surface as interesting from WG1–WG3 outputs?
- For those families, what pathways are enriched in their regulated genes, and how do DEG fold-changes overlay?
- Where does ATAC-seq accessibility flag differential binding of family motifs upstream of top regulated genes?
- Which GWAS SNPs sit in or near these TFs and their regulatory elements in the relevant cell type?

## Data

| Dataset | CRISPR pipeline | Energy distance |
|---|:---:|:---:|
| Hon WTC11 Cardiomyocyte | ready | ready |
| Huangfu HUES8 Definitive Endoderm | ready | ready |
| Huangfu HUES8 Embryonic Stem Cell | ready | ready |
| Gersbach WTC11 Hepatocyte | ready | - |
| Engreitz WTC11 Endothelial | - | - |

Family taxonomy uses `jaspar_tf_family` (curated, preferred) with fallback to `lambert_2018_dbd` (prefixed `DBD:`) — defined in [`../../reference/tf_metadata.tsv`](../../reference/tf_metadata.tsv). Per-dataset cards live under [`../../data/`](../../data/).

## Issues

- *[FILL IN issue link]*: ED calibration caveat — `distance_mean > NC max` is the calibration-robust significance proxy used for family scorecards.
- *[FILL IN issue link]*: Per-family ATAC/motif analysis and GWAS SNP overlay depend on the multiome stream and external GWAS catalog, neither staged in this folder.

## Working flow

Three steps, in order:

1. **Brainstorm** — sketch example figures, summary tables, and pseudocode that answer the questions above. Capture this in a notebook, doc, or notes in this folder.
2. **Execute** — run the analyses; commit code (notebooks, scripts, supporting docs) to this folder on GitHub.
3. **Share** — upload reusable outputs (figures, tables, intermediate data) to WG5's Synapse folder [`syn74954085`](https://www.synapse.org/Synapse:syn74954085) (mirrored `working_groups/wg5_tf_family_case_studies/` under [`syn64423137/2026_UTSW/`](https://www.synapse.org/Synapse:syn64423137)). Record the syn ID for each upload in the [Outputs](#outputs) table below so the next person can find it.

## Outputs

| Object | Syn ID | Description | Owner |
|---|---|---|---|
| *[FILL IN as outputs land]* | | | |
