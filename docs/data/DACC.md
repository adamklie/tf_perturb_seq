# DACC File Format Issues & Action Plan

Reference specs: [CPN FG File Formats](https://docs.google.com/document/d/139MbMlf22oi7VOh6UcmvooOgrmggyqlVcl487roNsZo/edit)

## Context

The DACC (Ian Whaling, Ben Hitz, Khine Lin) is ready to accept CPN FG file submissions to the portal/catalog. 
Mpathi has been working on uploading examples
This doc is a mostly LLM generated summary of the state of file on the portal and outstanding issues or questions.

---

## Gene Universe

**Spec:** Two required columns — `gene` (Ensembl ID, GENCODE V43) and `gene_symbol` (GENCODE V43). Represents the genes used as input to program inference (i.e., the top ~2000 HVGs selected before cNMF, not all expressed genes).

**On portal:** [IGVFFI6966LMRS](https://data.igvf.org/tabular-files/IGVFFI6966LMRS/) on analysis set [IGVFDS6332VCTO](https://data.igvf.org/analysis-sets/IGVFDS6332VCTO/). Status: **invalidated**.

**Validation errors:** 286 rows fail a `constraint-error` on the `gene` field (required but empty). These are HVGs that have a symbol but no Ensembl ID in the submitted file.

**Issues:**
1. 286/~2000 genes are missing Ensembl IDs. These should all be mappable to GENCODE V43 — need to investigate why the mapping failed.
2. Ben Hitz raised that the gene universe has 2,000 genes while the gene programs file has 25,312. This is expected — cNMF uses 2k HVGs as input but refits to all expressed genes. The naming is confusing but correct per the spec.

## Gene Programs

**Spec:** Required columns — `program_id` (format `GP_{number}`), `gene` (Ensembl ID, GENCODE V43), `gene_symbol`, `score` (program loading). Optional columns — `program_annotation`, `program_annotation_term`, `program_annotation_id`. Required headers: `#description` and `#sample_summary_short`. One file per cell type/sample.

**On portal:** [IGVFFI9914SKEC](https://data.igvf.org/tabular-files/IGVFFI9914SKEC/) on analysis set [IGVFDS6332VCTO](https://data.igvf.org/analysis-sets/IGVFDS6332VCTO/). Status: **invalidated**.

**Validation errors:**
1. **`incorrect-label` (6)** — Column names don't match the spec. The cNMF output uses its own naming rather than `program_id`, `gene`, `gene_symbol`, `score`.
2. **`missing-label` (1)** — A column expected by the validator is absent from the header.
3. **`type-error` (497)** — The `score` field (field 4) has non-numeric values in ~half the rows.
4. **`missing-cell` (496)** — Rows have fewer columns than the header, suggesting a malformed CSV (delimiter issue or ragged rows from optional columns).

**Issues:**
1. Column renaming needed — map cNMF output columns to spec schema (`program_id`, `gene`, `gene_symbol`, `score`).
2. The CSV itself appears malformed — type errors and missing cells suggest a structural problem (tab vs comma delimiter, or ragged rows). Need to inspect the actual file.
3. The file contains 25,312 genes (refitted loadings across all expressed genes, not just the 2k HVG input). This is correct cNMF behavior.
4. Unclear whether `#description` and `#sample_summary_short` headers are present — the validator may check these separately.

## Gene Program Regulators

**Spec:** Required columns — `program_id`, `gene` (Ensembl ID of perturbed gene), `gene_symbol`, `log2FC` (mean program scores perturbed vs reference), `reference_group`, `test_statistic`, `p_nominal_nlog10`, `fdr_nlog10`, `fdr_method`. Optional — `promoter_coordinates` (hg38). Required headers: `#description` and `#sample_summary_short`. Perturbation type (CRISPRi/CRISPRa) should be specified separately.

**On portal:** [IGVFFI9218RTDZ](https://data.igvf.org/tabular-files/IGVFFI9218RTDZ/) on analysis set [IGVFDS6332VCTO](https://data.igvf.org/analysis-sets/IGVFDS6332VCTO/). Status: **validated**.

**Issues:**
This passed validation, are we happy with contents?

## Global Differential Expression

**Spec:** The current file we are looking at for this is pySpade-specific, defined by the Hon lab ([PDF](https://api.data.igvf.org/documents/88d2ed1f-3c0a-4df4-b14f-c203d028899f/@@download/attachment/pySpade_DE_table_definitions.pdf)). Key columns: `gene_names`, `region` (perturbation), `fc`, `Significance_score`, `num_cell`, `cpm_perturb`, `cpm_bg`.

**On portal:** [IGVFFI5989UAVX](https://data.igvf.org/tabular-files/IGVFFI5989UAVX/) on analysis set [IGVFDS6332VCTO](https://data.igvf.org/analysis-sets/IGVFDS6332VCTO/). Status: **validated, released**. File format spec: `gary-hon:pySpade_DE_table_definitions`.

**Issues:**
1. **Is there a generic spec for global differential expression?** The pySpade format is method-specific. The CRISPR FG pipeline (Sceptre/CLEANSER/Perturbo) produces DE with a different schema (e.g., our calibrated trans results have `element_id`, `tested_gene_id`, `log2fc`, `empirical_pval`, etc.).
2. **Missing Ensembl IDs.** The pySpade format only has `gene_names`. Mpathi is adding Ensembl IDs for Ian to map catalog edges.
3. **GRN edge creation.** Ian wants to derive TF→gene edges from these DE results. The portal currently has no API endpoint for this (Ian: "we are missing that"). This is a DACC dependency.
4. **For benchmark datasets** run through the CRISPR FG pipeline, we produce calibrated trans results with a different schema. How these get submitted is an open question — do they go under `global differential expression` with a new spec, or under `predicted TF to gene regulatory interactions`?

## Other file types (not yet submitted)

The [CPN FG spec](https://docs.google.com/document/d/139MbMlf22oi7VOh6UcmvooOgrmggyqlVcl487roNsZo/edit) defines several additional file types that we have not submitted. Specs exist in the Google Doc but we haven't attempted to produce or upload any of these yet.

| File Type | Spec Defined? | Can We Produce It? | Notes |
|-----------|:---:|:---:|-------|
| Predicted TF-Gene Regulatory Interactions | Yes | Unclear | Need to discuss |
| TF Universe | Yes | Yes (from guide metadata) |  |
| Element Universe (.bed) | Yes | Yes (from guide metadata) |  |
