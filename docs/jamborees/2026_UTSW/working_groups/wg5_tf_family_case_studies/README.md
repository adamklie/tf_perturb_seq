# WG5 — TF family case studies

**Topic 2.4 / Figure 4.** Goal: deep-dive analysis of specific TFs / TF families with newly implicated roles in lineage differentiation (e.g., ZNF factors).

## Questions

From [`../../WORKING_GROUPS.md`](../../WORKING_GROUPS.md):

- Identify interesting TF families using outputs from Working Groups 1–3.
- Pathway analysis of top gene programs for these TFs, with DEG fold-change overlays.
- Use ATAC-seq data to assess differential accessibility of binding motifs upstream of top regulated genes.
- Highlight GWAS SNPs in or near these TFs or their regulatory elements linked to disease in the relevant cell type.

## Artifacts in this folder

| ID | File | Status | What it answers | Source data |
|---|---|---|---|---|
| WG5-A | `tf_family_activity_scorecard.tsv` | 🟡 partial | Per-family: n members, n significant in any lineage, mean / max distance, disease fraction, lineage specificity — feeds deep-dive selection | `tf_metadata.tsv` (`jaspar_tf_family`, `lambert_2018_dbd`) + per-dataset `pval_edist_full.csv` + WG3-A disease flags |
| WG5-B | `family_<family>_enrichment.tsv` (one file per family selected) | 🔴 blocked | KEGG / GO pathway enrichment for the family's regulated genes | cNMF programs + target gene lists (cNMF not yet run for all datasets) |

## Out of scope here

- **ATAC-seq accessibility + motif analysis**: multiome stream (same as WG4-C).
- **GWAS SNP overlay**: external resource (same as WG3-C).
