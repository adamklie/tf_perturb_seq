# CRISPR Pipeline Output Directory Structure

Reference GCP path: `gs://igvf-pertub-seq-pipeline-data/scratch/bioinfolucas/Benchmark_cleanser_unified/<DatasetName>/`

Based on audit of GaryHonDataset CLEANSER run (2026-03-25).

## Pipeline stages and their outputs

### 1. Sequence processing

| Directory | Contents | Useful for |
|-----------|----------|------------|
| `seqspecparser/` | Parsed seqspec files (barcode_file.txt, rna/crispr/hashing parsed seqspecs) | Verifying seqspec configuration |
| `seqspeccheck/` | Position tables + plots for guide/hashing seqspecs | Debugging seqspec issues |
| `downloadreference/` | Transcriptome index (.idx) + t2g mapping | Reference verification |
| `downloadgtf/` | GENCODE GTF | Reference verification |
| `createguideref/` | Guide index, mismatch ref, t2guide mapping | Guide reference verification |
| `createhashingref/` | Hashing index, mismatch ref, t2g mapping | HTO reference verification |

### 2. Mapping (per measurement set)

| Directory | Contents | Useful for |
|-----------|----------|------------|
| `mappingscrna/<MS>_ks_transcripts_out/` | kallisto bus output: **run_info.json** (read counts, alignment rate), matrix.ec, output.bus, counts_unfiltered/ | **Total reads, % aligned** — key for sequencing depth comparison |
| `mappingguide/<MS>_ks_guide_out/` | Same structure for guide reads | Guide sequencing depth |
| `mappinghashing/<MS>_ks_hashing_out/` | Same structure for HTO reads (Hon only) | HTO sequencing depth |

**Key file:** `run_info.json` contains `n_processed` (total reads), `n_pseudoaligned` (mapped reads), `p_pseudoaligned` (% mapped).

### 3. Preprocessing & QC

| Directory | Contents | Useful for |
|-----------|----------|------------|
| `preprocessanndata/` | `filtered_anndata.h5ad`, `rna_concatenated_adata.h5ad`, figures/ | Post-filtering cell counts, QC plots |
| `anndata/` | `concatenated_adata.h5ad` | Pre-QC combined data |
| `filter/` | Per-MS filtered h5ad files (hashing filtered) | Per-lane cell counts post-HTO filtering |

### 4. Guide assignment & HTO demultiplexing

| Directory | Contents | Useful for |
|-----------|----------|------------|
| `prepare/` | Per-MS MuData files (pre-guide-assignment) | Input to guide assignment |
| `guide/` | Per-MS MuData after guide assignment | Per-lane guide assignment results |
| `demultiplex/` | Per-MS hashing filtered/unfiltered demux h5ad | HTO demultiplexing QC |
| `hashing/` | Concatenated hashing demux h5ad | Combined HTO results |

### 5. MuData creation & merging

| Directory | Contents | Useful for |
|-----------|----------|------------|
| `createmudata/` | Combined MuData (mudata.h5mu), guide_concatenated_adata.h5ad | Pre-inference combined object |
| `mudata/` | concat_mudata.h5mu | Alternative combined MuData |

### 6. Inference

| Directory | Contents | Useful for |
|-----------|----------|------------|
| `sceptre/` | inference_mudata.h5mu, sceptre per-element/per-guide results, chunk_manifest.tsv | Sceptre inference results (run even in CLEANSER pipeline) |
| `inference/` | inference_mudata.h5mu, sceptre per-element/per-guide results | Additional inference outputs |
| `mergedresults/` | Final merged inference_mudata.h5mu, per-element and per-guide TSVs | **Primary inference outputs** |

### 7. Evaluation & TF enrichment

| Directory | Contents | Useful for |
|-----------|----------|------------|
| `evaluation/` | evaluation_output/, plots/ | Pipeline built-in evaluation |
| `tf/benchmark_output/` | TF enrichment tables (ChIP-seq validation) | TF target enrichment analysis |

### 8. Dashboard & final outputs

| Directory | Contents | Useful for |
|-----------|----------|------------|
| `pipeline_dashboard/` | **dashboard.html**, inference_mudata.h5mu, additional_qc/ (gene/guide/intended_target/trans metrics) | **Primary QC outputs we use in slides** |
| `pipeline_outputs/` | Final inference_mudata.h5mu, perturbo cis/trans per-element/per-guide TSVs | **Primary inference outputs for downstream** |
| `additional/` | additional_qc/ (duplicate of pipeline_dashboard QC) | Backup QC |
| `pipeline_info/` | params JSON, software versions YAML | **Actual config used in run** |

## Extracting sequencing depth per dataset

From `mappingscrna/*/run_info.json`:

| Dataset | Measurement Sets | Total scRNA Reads | Reads (M) |
|---------|:---:|------------------:|----------:|
| Hon | 3 | 1,946,507,920 | 1,947M |
| Huangfu | 4 | 1,496,571,113 | 1,497M |
| Engreitz | 16 | 6,425,724,497 | 6,426M |
| Gersbach GEM-Xv3 | 4 | 6,404,395,504 | 6,404M |
| Gersbach HTv2 | 4 | 6,382,797,522 | 6,383M |

## Read counts by modality (from `mapping*/*/run_info.json`)

| Dataset | scRNA Reads | scRNA % Aligned | Guide Reads | Guide % Aligned | HTO Reads | HTO % Aligned |
|---------|:-----------:|:---:|:-----------:|:---:|:---------:|:---:|
| Hon (3 MS) | 1,947M | 91.5% | 568M | 62.5% | 161M | 94.9% |
| Huangfu (4 MS) | 1,497M | 90.3% | 381M | **16.1%** | — | — |
| Engreitz (16 MS) | 6,426M | 91.0% | 768M | 90.9% | — | — |
| Gersbach GEM-Xv3 (4 MS) | 6,404M | 85.8% | 477M | 95.8% | — | — |
| Gersbach HTv2 (4 MS) | 6,383M | 91.1% | 692M | 95.9% | — | — |

**Note:** Huangfu guide alignment is 16.1% — this is a known aspect of 10x 3' v3 chemistry that has already been troubleshot. Not a concern.

Full per-MS breakdown: `results/cross_tech_comparison/per_ms_read_counts.tsv`

## Filtered anndata sizes (proxy for cell counts)

| Dataset | `filtered_anndata.h5ad` size |
|---------|:---:|
| Hon | 2.71 GB |
| Huangfu | 2.16 GB |
| Engreitz | 3.79 GB |
| Gersbach GEM-Xv3 | 2.83 GB |
| Gersbach HTv2 | 2.84 GB |

**TODO:** Extract actual cell counts from these h5ad files or from the final `inference_mudata.h5mu`.

## Other pipeline outputs worth downloading

| Directory | Contents | Value |
|-----------|----------|-------|
| `preprocessanndata/figures/` | knee_plot_scRNA.png, scatterplot_scrna.png, violinplot_scrna.png | Upstream QC — shows cell calling, UMI distributions before guide assignment |
| `evaluation/plots/` | trans_perturbo_barplot_direct_vs_control.png, precision_recall_roc.png, volcano_plot.png | Pipeline's built-in evaluation of inference results |
| `evaluation/evaluation_output/` | bedgraph/bedpe files for cis/trans perturbo and sceptre | Could use for genomic visualizations |
| `seqspeccheck/guide_position_table.csv` | Per-sample guide detection HitRatio and configuration scoring | Explains guide alignment differences across technologies |
| `filter/` | Per-MS hashing filtered h5ad files (Hon only) | HTO demultiplexing QC |

## TODO: Build reusable extraction script

This analysis should be automated into a script that can be re-run whenever we get new pipeline outputs. Target location: `scripts/extract_pipeline_metrics.py`

Inputs: GCP base path, list of datasets
Outputs:
- Per-MS read counts TSV (all modalities)
- Pipeline config comparison table
- Upstream QC figures (downloaded)
- Cell counts at each stage

## Pipeline configuration comparison (from GCP `pipeline_info/params_*.json`)

| Parameter | Hon | Huangfu | Engreitz | Gersbach GEM-Xv3 | Gersbach HTv2 |
|-----------|-----|---------|----------|-------------------|---------------|
| `QC_min_genes_per_cell` | 500 | 500 | 500 | **1000** | **1000** |
| `QC_pct_mito` | **15** | 20 | 20 | 20 | 20 |
| `is_10x3v3` | false | **true** | false | false | false |
| `capture_method` | CROP-seq | CROP-seq | CROP-seq | **direct-capture** | CROP-seq |
| `ENABLE_DATA_HASHING` | **true** | false | false | false | false |
| `reverse_complement_guides` | true | **false** | true | true | true |
| `GUIDE_ASSIGNMENT_method` | cleanser* | cleanser* | cleanser* | cleanser* | cleanser* |

*CLEANSER run shown. Sceptre run uses `sceptre` for guide assignment; all other params identical except Gersbach HTv2 Sceptre run has no params (may not have completed).
