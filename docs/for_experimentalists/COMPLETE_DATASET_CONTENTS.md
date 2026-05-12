# Complete dataset contents

What a "fully baked" CRISPRi Perturb-seq dataset should contain when handed off to a collaborator. Use this to spot-check whether a bundle is missing something.

This doc cross-references three things:

1. **Our internal CRISPR-pipeline outputs** — what the IGVF pipeline emits.
2. **The IGVF DACC submission spec** — what we have to deliver back to the consortium. See [`docs/data/DACC.md`](../data/DACC.md) for the source-of-truth audit.
3. **The reference IGVF analysis set [IGVFDS6332VCTO](https://data.igvf.org/analysis-sets/IGVFDS6332VCTO/)** — Gary Hon's WTC11 cardiomyocyte dataset; a fully released example showing what a finished submission looks like.

---

## The complete bundle, layer by layer

### Layer 1 — Reference data (shared across all datasets in a project)

| File | Format | What | DACC name | Status |
|---|---|---|---|---|
| Genome annotation (GENCODE V43 + IGVF release) | `.gtf.gz` | GTF used by every dataset | (standard external file) | Pinned per project |
| Guide library | `.csv.gz` / `.tsv.gz` | Every guide + spacer + intended target + genomic location | Used to generate `TF Universe` + `Element Universe` (see below) | Pinned per project |
| TF metadata | `.tsv` | One row per TF target gene; HGNC + Lambert 2018 DBD + JASPAR family | — | Generated per project |
| Experimental metadata | `.tsv` | One row per dataset; lab + cell line + chemistry + measurement-set count | — | Generated per project |

### Layer 2 — Raw / processed count data (per dataset)

| File | Format | What | Pipeline source | DACC content_type |
|---|---|---|---|---|
| Raw feature-barcode matrix | `.hdf5` or `.pkl` | Sparse cell × feature counts pre-filter | CellRanger / kallisto output | `Raw feature barcode matrix` |
| Filtered feature-barcode matrix | `.hdf5` | Sparse cell × feature counts post-filter (whitelist + min counts) | CellRanger / pipeline filter step | `Filtered feature barcode matrix` |
| Sparse gene count matrix | `.h5ad` | AnnData with the QC-passed cells × genes | `preprocessanndata/` | `Sparse gene count matrix` |

In our pipeline, these usually live inside `pipeline_dashboard/` (the figures/dashboard staging area) and `pipeline_outputs/`.

### Layer 3 — Inference MuData + perturbo (per dataset)

| File | Format | What | Pipeline source | DACC content_type |
|---|---|---|---|---|
| `inference_mudata.h5mu` | `.h5mu` | The integrated MuData; cells × genes + cells × guides + cells × HTO | `pipeline_outputs/`, `mergedresults/`, `pipeline_dashboard/` | — (not a DACC-deliverable; downstream tables are) |
| `perturbo_cis_per_element.tsv.gz` | `.tsv.gz` | Per-perturbation cis-DE — did the knockdown work? | `pipeline_outputs/` | maps to `Local differential expression` |
| `perturbo_trans_per_element.tsv.gz` | `.tsv.gz` | Per-perturbation trans-DE — genome-wide consequence | `pipeline_outputs/` | maps to `Global differential expression` |
| `perturbo_cis_per_guide.tsv.gz` | `.tsv.gz` | Same but one row per (guide, gene) | `pipeline_outputs/` | (sub-deliverable of cis) |
| `perturbo_trans_per_guide.tsv.gz` | `.tsv.gz` | Same but one row per (guide, gene) | `pipeline_outputs/` | (sub-deliverable of trans) |
| `dashboard.html` + `figures/` + `additional_qc/` | mixed | QC dashboard + plots | `pipeline_dashboard/` | — |
| `pipeline_info/params_<run>.json` | `.json` | Nextflow params used | `pipeline_info/` | — |
| `pipeline_info/software_versions.yml` | `.yml` | Tool versions | `pipeline_info/` | — |

### Layer 4 — Energy distance (per dataset)

| File | Format | What | Pipeline source | DACC content_type |
|---|---|---|---|---|
| `pval_edist_full.csv` | `.csv` | Per-target energy distance + 20 permutation p-values | Energy-distance pipeline Step 2 | — (no current DACC spec) |
| `targeting_outlier_table.csv` | `.csv` | Outlier-guide flags | Step 1 | — |
| `non_targeting_outlier_table.csv` | `.csv` | Outlier-NTC flags | Step 1 | — |
| `target_by_target_matrix.csv` (optional) | `.csv` | Pairwise TF × TF distance matrix | Step 3 | — |
| `edist_embedding_info.csv` (optional) | `.csv` | 2D t-SNE of TFs | Step 3 | — |
| `figures/*.pdf`, `config_step2.json`, `logs/` | mixed | Provenance | All steps | — |

### Layer 5 — cNMF / gene programs (per dataset)

| File | Format | What | Pipeline source | DACC content_type |
|---|---|---|---|---|
| `Inference.gene_spectra_score.k_<k>.dt_2_0.txt` | `.tsv` | Gene loadings (programs × genes, z-scored) at selected k | Stage 1 inference | maps to `Gene programs` (after reformat — see DACC.md) |
| `Inference.gene_spectra_tpm.k_<k>.dt_2_0.txt` | `.tsv` | Gene loadings (TPM-normalized) | Stage 1 inference | (alternate) `Gene programs` |
| `Inference.usages.k_<k>.dt_2_0.consensus.txt` | `.tsv` | Cell × program usages | Stage 1 inference | — |
| `Inference.overdispersed_genes.txt` | text | HVG list cNMF was fit on | Stage 1 inference | maps to `Gene universe` (after reformat) |
| `<k>_perturbation_association_results_all.txt` | `.tsv` | Per-program × TF perturbation log2FC + q-value | Stage 2 evaluation | maps to `Gene program regulators` (after reformat) |
| `<k>_geneset_enrichment.txt` | `.tsv` | MSigDB enrichments per program | Stage 2 evaluation | (annotation; optional in DACC spec) |
| `<k>_GO_term_enrichment.txt` | `.tsv` | GO enrichments per program | Stage 2 evaluation | (annotation; optional in DACC spec) |
| `k_selection.png` + `k_selection_stats.df.npz` | mixed | k-selection plot + stats | Stage 3a | — |
| `Inference/clustering.k_<K>.dt_2_0.png` for every K | `.png` | Sweep clustering pngs (provenance) | Stage 1 | — |
| `README.txt` | text | k-selection rationale + group decision notes | — | — |

### Layer 6 — DACC-spec deliverables (derived; for the IGVF portal)

These don't come out of the pipelines directly — they're reformatted from upstream outputs to match the DACC schemas. The reference set IGVFDS6332VCTO has all of these.

| File | Format | What | Derived from | Status |
|---|---|---|---|---|
| `TF Universe` | `.tsv` | One row per unique TF whose promoter is in the guide library | Guide library (Layer 1) | ✅ generable now |
| `Element Universe` | `.bed` | BED of every genomic element targeted by the library | Guide library | ✅ generable now |
| `Gene Universe` | `.tsv` | Two cols: `gene` (ENSG), `gene_symbol`. The HVG list cNMF was fit on. | `overdispersed_genes.txt` | Needs cNMF run |
| `Gene Programs` | `.tsv` | Cols: `program_id`, `gene`, `gene_symbol`, `score` (+ optional annotation). One file per cell type. | `gene_spectra_score.k_<sel>.dt_2_0.txt` | Needs cNMF run + reformat |
| `Gene Program Regulators` | `.tsv` | Cols: `program_id`, `gene` (ENSG of perturbed gene), `gene_symbol`, `log2FC`, `reference_group`, `test_statistic`, `p_nominal_nlog10`, `fdr_nlog10`, `fdr_method`. | `<k>_perturbation_association_results_all.txt` | Needs cNMF run + reformat |
| `Global differential expression` | `.tsv` | Per-perturbation × per-gene log2FC + significance. Two specs exist: pySpade-flavor (the Hon lab default) and a generic CRISPR-pipeline-flavor. | `perturbo_trans_per_element.tsv.gz` | Needs reformat per spec |
| `Local differential expression` | `.tsv` | Per-guide × per-target-gene (cis) DE | `perturbo_cis_per_element.tsv.gz` | Needs reformat |
| `Predicted TF-Gene Regulatory Interactions` | `.tsv` (spec defined; we don't produce it yet) | Edges TF → gene with weights | derived from `perturbo_trans_per_element.tsv.gz` (FDR-filtered) | Pending implementation |

---

## Mapping to the reference IGVF analysis set

For comparison, the [IGVFDS6332VCTO](https://data.igvf.org/analysis-sets/IGVFDS6332VCTO/) analysis set (Gary Hon, WTC11 cardiomyocyte) has these released files. Use this as a template for what a finished submission looks like:

| IGVF accession | File format | Content type | Status |
|---|---|---|---|
| [`IGVFFI3617IJOW`](https://data.igvf.org/tabular-files/IGVFFI3617IJOW/) | hdf5 | Filtered feature barcode matrix | Released ✅ |
| [`IGVFFI4735RXPI`](https://data.igvf.org/tabular-files/IGVFFI4735RXPI/) | pkl | Raw feature barcode matrix | Released ✅ |
| [`IGVFFI7637STPX`](https://data.igvf.org/tabular-files/IGVFFI7637STPX/) | h5ad | Sparse gene count matrix | Released ✅ |
| [`IGVFFI5989UAVX`](https://data.igvf.org/tabular-files/IGVFFI5989UAVX/) | csv | Global differential expression (pySpade) | Released ✅ |
| [`IGVFFI7298IERA`](https://data.igvf.org/tabular-files/IGVFFI7298IERA/) | csv | Local differential expression | Released ✅ |
| [`IGVFFI0830FXFI`](https://data.igvf.org/tabular-files/IGVFFI0830FXFI/) | csv | Global differential expression (second) | In progress |
| [`IGVFFI6966LMRS`](https://data.igvf.org/tabular-files/IGVFFI6966LMRS/) | csv | Gene universe | In progress |
| [`IGVFFI9218RTDZ`](https://data.igvf.org/tabular-files/IGVFFI9218RTDZ/) | csv | Gene program regulators | In progress |
| [`IGVFFI9914SKEC`](https://data.igvf.org/tabular-files/IGVFFI9914SKEC/) | csv | Gene programs | In progress |

Workflows attached: **Hon Perturb-seq Workflow** (IGVFWF0190YALU; released) + **Hon cNMF Workflow** (IGVFWF8652HUAV; in progress).

---

## "Is my bundle complete?" — checklist

Mark each box for your dataset:

- [ ] CRISPR-pipeline outputs
  - [ ] `pipeline_dashboard/` (dashboard.html + additional_qc/ + figures/)
  - [ ] `pipeline_outputs/` (inference_mudata.h5mu + perturbo TSVs)
  - [ ] `pipeline_info/` (params + software versions)
- [ ] Energy-distance outputs (if run)
  - [ ] `pval_edist_full.csv`
  - [ ] Outlier tables for targeting + non-targeting
  - [ ] Figures + config + logs
- [ ] cNMF outputs (if run)
  - [ ] Stage-1 inference at selected k
  - [ ] Stage-2 evaluation at selected k
  - [ ] Sweep-as-provenance (all-k clustering pngs + k_selection.png)
  - [ ] `README.txt` with k-selection rationale
- [ ] DACC deliverables (if submitting to the portal)
  - [ ] TF Universe
  - [ ] Element Universe (BED)
  - [ ] Gene Universe
  - [ ] Gene Programs
  - [ ] Gene Program Regulators
  - [ ] Global differential expression (project's chosen spec)
  - [ ] Local differential expression

Any unchecked box should be tracked as a known gap on the dataset's status page.

---

## See also

- [`docs/data/DACC.md`](../data/DACC.md) — DACC file-format audit (what specs exist, what's validated, what's open).
- [`docs/jamborees/<event>/schemas/`](../jamborees/) — per-event JSON schemas for the canonical pipeline outputs (CRISPR pipeline, energy distance, cNMF).
- [`UNDERSTAND_CRISPR_OUTPUTS.md`](UNDERSTAND_CRISPR_OUTPUTS.md), [`UNDERSTAND_ENERGY_DISTANCE.md`](UNDERSTAND_ENERGY_DISTANCE.md), [`UNDERSTAND_CNMF_OUTPUTS.md`](UNDERSTAND_CNMF_OUTPUTS.md) — interpretation guides for each output type.
