# WG3 — Disease and GWAS

Topic 2.2 / Figure 2. The goal is to connect TF regulatory activity to human disease: which TFs regulate disease genes in each lineage, where convergent vs divergent activity shows up across lineages, and how GWAS variants relate to the regulatory targets of these TFs. End-of-jamboree deliverables include a disease-TF activity panel for Figure 2 and a shortlist of TFs nominated for the case-study working groups.

## Questions

- Which TFs regulate disease and GWAS genes in each lineage?
- For TFs that are disease genes in multiple lineages, is their activity convergent or divergent?
- Which GWAS variants sit near important TFs or near the regulatory elements upstream of their downstream targets, and how should they be annotated?

## Data

| Dataset | CRISPR pipeline | Energy distance |
|---|:---:|:---:|
| Hon WTC11 Cardiomyocyte | caveat | ready |
| Huangfu HUES8 Definitive Endoderm | ready | caveat |
| Huangfu HUES8 Embryonic Stem Cell | ready | caveat |
| Gersbach WTC11 Hepatocyte | caveat | blocked |
| Engreitz WTC11 Endothelial | blocked | blocked |

Gene-disease associations come from HPO's `genes_to_disease.txt` (MONDO + OMIM), cached at [`../../reference/gene_disease_associations.tsv`](../../reference/gene_disease_associations.tsv). Per-dataset cards live under [`../../data/`](../../data/).

## Issues

- *[FILL IN issue link]*: ED calibration caveat — `distance_mean` is the trustworthy signal; `pval_mean` is anti-conservative for the Huangfu runs. Use `distance_mean > NC max` as the calibration-robust significance proxy.
- *[FILL IN issue link]*: GWAS variant overlay relies on an external GWAS catalog not yet staged in the repo.

## Working flow

Three steps, in order:

1. **Brainstorm** — sketch example figures, summary tables, and pseudocode that answer the questions above. Capture this in a notebook, doc, or notes in this folder.
2. **Execute** — run the analyses; commit code (notebooks, scripts, supporting docs) to this folder on GitHub.
3. **Share** — upload reusable outputs (figures, tables, intermediate data) to WG3's Synapse folder [`syn74954083`](https://www.synapse.org/Synapse:syn74954083) (mirrored `working_groups/wg3_disease_gwas/` under [`syn64423137/2026_UTSW/`](https://www.synapse.org/Synapse:syn64423137)). Record the syn ID for each upload in the [Outputs](#outputs) table below so the next person can find it.

## Outputs

| Object | Syn ID | Description | Owner |
|---|---|---|---|
| *[FILL IN as outputs land]* | | | |
