# Complete dataset contents

What a fully-baked CRISPRi Perturb-seq dataset on Synapse should contain when handed off to a working group. Use this to spot-check whether your dataset is missing something.

This doc walks through the **output-type subfolders** each dataset has on Synapse, using Huangfu HUES8 Definitive Endoderm ([`syn74834951`](https://www.synapse.org/Synapse:syn74834951)) as the worked example. Every dataset on Synapse follows the same layout — the subfolder names are stable across labs and cell types, even when some output types haven't been produced yet.

```
<dataset>/                          (e.g. Huangfu DE = syn74834951)
├── crispr_pipeline/                CRISPR pipeline outputs (pipeline_outputs, pipeline_dashboard, pipeline_info)
├── calibration/                    DEG-calibrated DE tables (NTC-derived null; per-TF BH FDR)
├── qc/                             Pipeline-output QC (intended-target repression, mapping QC for guides + genes)
├── energy_distance/                Per-target energy distance + outlier tables
└── cnmf/                           Gene programs (cNMF inference + evaluation + interpretation)
```

Reference data shared across all datasets in the project (genome, guide library, TF metadata, experimental metadata) lives at the project root, not inside each dataset — see [Reference data](#reference-data) below.

---

## Reference data (shared across datasets)

These live at the project level, not inside each dataset folder. See [`../README.md`](../README.md) "Synapse IDs at a glance" for syn IDs.

| File | Format | What |
|---|---|---|
| Genome annotation (GENCODE V43 + IGVF release) | `.gtf.gz` | Pinned per project; used by every dataset |
| Guide library | `.csv.gz` / `.tsv.gz` | One row per guide; spacer + intended target + genomic location |
| TF metadata (`tf_metadata.tsv`) | `.tsv` | One row per TF target gene; HGNC + Lambert 2018 DBD + JASPAR family |
| Experimental metadata (`experimental_metadata.tsv`) | `.tsv` | One row per dataset; lab + cell line + chemistry + measurement-set count |

---

## CRISPR pipeline (`crispr_pipeline/`)

Worked example: [`syn74834952`](https://www.synapse.org/Synapse:syn74834952). Three sub-directories:

| Subfolder | What's inside |
|---|---|
| `pipeline_outputs/` | `inference_mudata.h5mu` (cells × genes + cells × guides + cells × HTO) and `perturbo_{cis,trans}_per_{element,guide}_output.tsv.gz` (raw DE tables) |
| `pipeline_dashboard/` | `dashboard.html`, `additional_qc/`, `evaluation_output/`, `figures/`, `guide_seqSpec_plots/`, `svg/`, and a copy of `inference_mudata.h5mu` |
| `pipeline_info/` | `nf_core_pipeline_software_versions.yml` + per-run `params_<timestamp>.json` |

Key files to know:

| File | Format | What |
|---|---|---|
| `inference_mudata.h5mu` | `.h5mu` | Integrated MuData; load with `mudata.read_h5mu()` |
| `perturbo_cis_per_element_output.tsv.gz` | `.tsv.gz` | Per-perturbation cis-DE — did the knockdown work? |
| `perturbo_trans_per_element_output.tsv.gz` | `.tsv.gz` | Per-perturbation trans-DE — genome-wide consequence |
| `perturbo_{cis,trans}_per_guide_output.tsv.gz` | `.tsv.gz` | Same but one row per (guide, gene) |
| `dashboard.html` | `.html` | Click-through QC report |
| `params_<run>.json` | `.json` | Nextflow params used for the run |

Interpretation guide: [`CRISPR.md`](CRISPR.md).

---

## Calibration (`calibration/`)

Worked example: [`syn74920615`](https://www.synapse.org/Synapse:syn74920615). Output of post-processing the CRISPR pipeline's perturbo DE tables to recalibrate p-values against the non-targeting-control null (eCDF or t-fit), then per-TF BH FDR.

| File | Format | What |
|---|---|---|
| `<prefix>_calibrated_all_results.tsv` | `.tsv` | All perturbation × gene tests, calibrated |
| `<prefix>_calibrated_cis_results.tsv` | `.tsv` | Cis tests only (perturbation near its intended target) |
| `<prefix>_calibrated_direct_target_results.tsv` | `.tsv` | Subset restricted to the intended-target gene |
| `<prefix>_calibrated_trans_results.tsv` | `.tsv` | Trans tests only |

`<prefix>` is the dataset + pipeline-run label (e.g. `Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq_muddy_penguin`).

Why calibration matters: raw perturbo p-values can be anti-conservative on some chemistries. The calibrated tables are what downstream WGs should use for FDR-controlled hit calling. Plan + reference implementation: [`../../../../src/tf_perturb_seq/inference/calibrate.py`](../../../../src/tf_perturb_seq/inference/calibrate.py). Drives [Issue #11](https://github.com/adamklie/tf_perturb_seq/issues/11).

---

## QC (`qc/`)

Worked example: [`syn74918479`](https://www.synapse.org/Synapse:syn74918479). Three sub-directories of plots + metric tables for spot-checking a pipeline run.

### `intended_target/`

Does the knockdown work for each TF whose target is in the library?

| File | Format | What |
|---|---|---|
| `<prefix>_intended_target_metrics.tsv` | `.tsv` | Per-TF: target gene, observed log2FC, p-value, FDR, hit flag |
| `<prefix>_intended_target_results.tsv` | `.tsv` | Per-test results table feeding the metrics |
| `<prefix>_intended_target_log2fc_distribution.png` | `.png` | Distribution of intended-target log2FCs (NTC vs targeting) |
| `<prefix>_intended_target_volcano.png` | `.png` | Volcano of TFs (effect size × significance) |
| `<prefix>_intended_target_roc_pr_curves.png` | `.png` | ROC + PR curves for intended-target recovery |

### `mapping_gene/`

Per-cell gene-UMI distributions and batch effects.

| File | Format | What |
|---|---|---|
| `<prefix>_gene_metrics.tsv` | `.tsv` | Per-cell + per-batch gene-mapping QC metrics |
| `<prefix>_gene_histograms.png` | `.png` | UMI distributions across cells |
| `<prefix>_gene_histograms_by_batch.png` | `.png` | Same, split by batch |
| `<prefix>_gene_knee_plot.png` | `.png` | Cell-barcode knee for the gene library |
| `<prefix>_gene_cells_per_batch.png` | `.png` | Cell counts per batch |

### `mapping_guide/`

Same idea for the guide library.

| File | Format | What |
|---|---|---|
| `<prefix>_guide_metrics.tsv` | `.tsv` | Per-cell + per-batch guide-mapping QC metrics |
| `<prefix>_guide_per_guide_capture.tsv` | `.tsv` | Per-guide capture rate across cells |
| `<prefix>_guide_histograms.png` | `.png` | Guide-UMI distributions |
| `<prefix>_guide_histograms_by_batch.png` | `.png` | Same, split by batch |
| `<prefix>_guide_knee_plot.png` | `.png` | Cell-barcode knee for the guide library |

Use the `intended_target_metrics.tsv` to decide whether a dataset's effect-size calls are trustworthy before pulling in its calibration outputs downstream.

---

## Energy distance (`energy_distance/`)

Worked example: [`syn74883327`](https://www.synapse.org/Synapse:syn74883327). One per-target distance + outlier set per dataset.

| File | Format | What |
|---|---|---|
| `pval_edist_full.csv` | `.csv` | Per-target energy distance + permutation p-values |
| `targeting_outlier_table.csv` | `.csv` | Outlier-guide flags |
| `non_targeting_outlier_table.csv` | `.csv` | Outlier-NTC flags |
| `discordance_gRNA.csv` | `.csv` | Per-guide discordance summary |
| `target_by_target_matrix.csv` (optional) | `.csv` | Pairwise TF × TF distance matrix |
| `edist_embedding_info.csv` (optional) | `.csv` | 2D embedding of TFs |
| `image/*.pdf` | `.pdf` | Distribution + cutoff plots |
| `config{1_2,3}.json` + `logs/` | mixed | Provenance |

Interpretation guide: [`ENERGY_DISTANCE.md`](ENERGY_DISTANCE.md).

---

## cNMF (`cnmf/`)

Worked example: [`syn74893844`](https://www.synapse.org/Synapse:syn74893844). Stage-1 inference + Stage-2 evaluation + Stage-3 k-selection + post-hoc interpretation, all under one folder.

| Subfolder | What's inside |
|---|---|
| `Inference/` | `Inference.gene_spectra_score.k_<k>.dt_2_0.txt` (gene loadings), `Inference.usages.k_<k>.dt_2_0.consensus.txt` (cell × program usages), `Inference.overdispersed_genes.txt`, k-selection plot + stats, per-k clustering pngs |
| `Evaluation/` | One subfolder per `<k>_<dt>` with `<k>_perturbation_association_results_all.txt`, `<k>_geneset_enrichment.txt`, `<k>_GO_term_enrichment.txt` |
| `Interpretation/Summary_table/` | Per-program summary tables |
| `Plot/` | `k_selection/` (k-selection figure) + `Perturb_gene_<k>_<dt>/` (per-program × per-TF perturbation plots) |
| `Config/` | YAML configs for every Stage 1 / 2a / 3a / 3c run (provenance) |
| `adata/` | `cNMF_<k>_<dt>.h5mu` — MuData with usages + scores joined to the upstream `inference_mudata.h5mu` |
| `README_<run>.txt` | `.txt` | k-selection rationale + group decision notes |

Interpretation guide: [`CNMF.md`](CNMF.md).

---

## Mapping to the reference IGVF analysis set

For comparison, the [IGVFDS6332VCTO](https://data.igvf.org/analysis-sets/IGVFDS6332VCTO/) analysis set (Hon WTC11 cardiomyocyte) is a fully released submission. Use it as a template for what a finished IGVF-portal-deliverable set looks like:

| IGVF accession | File format | Content type | Status |
|---|---|---|---|
| [`IGVFFI3617IJOW`](https://data.igvf.org/tabular-files/IGVFFI3617IJOW/) | hdf5 | Filtered feature barcode matrix | Released |
| [`IGVFFI4735RXPI`](https://data.igvf.org/tabular-files/IGVFFI4735RXPI/) | pkl | Raw feature barcode matrix | Released |
| [`IGVFFI7637STPX`](https://data.igvf.org/tabular-files/IGVFFI7637STPX/) | h5ad | Sparse gene count matrix | Released |
| [`IGVFFI5989UAVX`](https://data.igvf.org/tabular-files/IGVFFI5989UAVX/) | csv | Global differential expression (pySpade) | Released |
| [`IGVFFI7298IERA`](https://data.igvf.org/tabular-files/IGVFFI7298IERA/) | csv | Local differential expression | Released |
| [`IGVFFI0830FXFI`](https://data.igvf.org/tabular-files/IGVFFI0830FXFI/) | csv | Global differential expression (second) | In progress |
| [`IGVFFI6966LMRS`](https://data.igvf.org/tabular-files/IGVFFI6966LMRS/) | csv | Gene universe | In progress |
| [`IGVFFI9218RTDZ`](https://data.igvf.org/tabular-files/IGVFFI9218RTDZ/) | csv | Gene program regulators | In progress |
| [`IGVFFI9914SKEC`](https://data.igvf.org/tabular-files/IGVFFI9914SKEC/) | csv | Gene programs | In progress |

Workflows attached: **Hon Perturb-seq Workflow** (IGVFWF0190YALU; released) + **Hon cNMF Workflow** (IGVFWF8652HUAV; in progress). The DACC-spec reformatting audit lives at [`../../../data/DACC.md`](../../../data/DACC.md).

---

## "Is my dataset complete?" — checklist

Mark each box for your dataset:

- [ ] `crispr_pipeline/`
  - [ ] `pipeline_outputs/` (`inference_mudata.h5mu` + perturbo TSVs)
  - [ ] `pipeline_dashboard/` (`dashboard.html` + `additional_qc/` + `figures/`)
  - [ ] `pipeline_info/` (params + software versions)
- [ ] `calibration/`
  - [ ] `<prefix>_calibrated_all_results.tsv`
  - [ ] `<prefix>_calibrated_cis_results.tsv`
  - [ ] `<prefix>_calibrated_direct_target_results.tsv`
  - [ ] `<prefix>_calibrated_trans_results.tsv`
- [ ] `qc/`
  - [ ] `intended_target/` (metrics + results + 3 plots)
  - [ ] `mapping_gene/` (metrics + 4 plots)
  - [ ] `mapping_guide/` (metrics + per-guide capture + 3 plots)
- [ ] `energy_distance/` (if run)
  - [ ] `pval_edist_full.csv`
  - [ ] Outlier tables for targeting + non-targeting
  - [ ] `image/` + `logs/` + config JSON(s)
- [ ] `cnmf/` (if run)
  - [ ] `Inference/` at the selected k (plus all-k pngs for provenance)
  - [ ] `Evaluation/<k>_<dt>/` for the selected k
  - [ ] `Plot/k_selection/` + `Plot/Perturb_gene_<k>_<dt>/`
  - [ ] `Interpretation/Summary_table/`
  - [ ] `README_<run>.txt` with k-selection rationale

Any unchecked box should be tracked as a known gap on the dataset's status page ([`../data/<dataset>/README.md`](../data/)).

---

## See also

- [`../../../data/DACC.md`](../../../data/DACC.md) — DACC file-format audit (what specs exist, what's validated, what's open).
- [`../data/schemas/`](../data/schemas/) — per-output JSON schemas (CRISPR pipeline, energy distance, cNMF, experimental metadata, guides).
- [`CRISPR.md`](CRISPR.md), [`ENERGY_DISTANCE.md`](ENERGY_DISTANCE.md), [`CNMF.md`](CNMF.md) — interpretation guides for each output type.
