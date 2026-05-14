# Pipeline QC — Production Datasets

Cross-dataset pipeline QC for the five TFP3 production datasets, in preparation for WG1 Figure 1 panels at the [2026 UTSW jamboree](../../../README.md). Mirrors the structure of the technology-benchmark cross-comparison ([docs/manuscripts/CRISPRi_tech_benchmark/bin/2_qc/](../../../../../manuscripts/CRISPRi_tech_benchmark/bin/2_qc/)), but applied to differentiated lineages rather than a single WTC11 cell line.

**Goal.** Harmonize general statistics (cell counts, gRNA/scRNA MOI, %mito) and quality metrics across the five production datasets and surface a consistent set of comparison panels (gene/guide metrics, intended-target repression, guide capture) for the jamboree.

## Datasets

| Dataset ID | Canonical run | Notes |
|---|---|---|
| `Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq` | `2026_04_19_no_spacer` | QC TSVs in `2026_04_19_no_spacer/qc/`; h5mu downloaded from Synapse [syn74520421](https://www.synapse.org/Synapse:syn74520421) and lives at `synapse_inference_mudata/inference_mudata.h5mu` (no local CRISPR pipeline run) |
| `Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq` | `muddy_penguin` | GCS `2026_04_09/outs/muddy_penguin/` |
| `Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq` | `sceptre_v1` | GCS `2026_04_13/outs/sceptre_v1/` |
| `Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq` | `sara_synapse_syn74842722` | Sara's h5mu from Synapse [syn74842722](https://www.synapse.org/Synapse:syn74842722); only QC re-run locally, no full pipeline |
| `Engreitz_WTC11-endothelial-cells_TF-Perturb-seq` | *TBD* | Not yet onboarded; placeholder card at [data/Engreitz_WTC11-endothelial-cells_TF-Perturb-seq/](../../../data/Engreitz_WTC11-endothelial-cells_TF-Perturb-seq/) |

Per-dataset cards live under [`docs/jamborees/2026_UTSW/data/`](../../../data/).

## Pipeline parameters

The four onboarded datasets were run with three different `(operator, platform, library-chemistry)` combinations. The `params_*.json` for every run has been pulled to the local repo (Synapse for Hon CM and Gersbach Hep; GCS for the two Huangfu runs) and is referenced from [`manifests/production_pipeline_params.tsv`](manifests/production_pipeline_params.tsv). A rendered version of the table (lineage-colored headers, outlier cells highlighted) lives at [`results/cross_production_qc/pipeline_params_table.pdf`](results/cross_production_qc/pipeline_params_table.pdf) (regenerate with [`scripts/3_pipeline_params_table.py`](scripts/3_pipeline_params_table.py)).

![Pipeline parameters by dataset](results/cross_production_qc/pipeline_params_table.png)

Differing parameters across runs:

| Param | Hon CM (`2026_04_19_no_spacer`) | Huangfu DE (`muddy_penguin`) | Huangfu ESC (`sceptre_v1`) | Gersbach Hep (`sara_synapse_syn74842722`) |
|---|---|---|---|---|
| Operator / platform | Weizhou Qian, UMich HPC | Adam, GCP Batch | Yan Yang, GCP Batch | Sara Geraghty, Duke HPC |
| `ENABLE_DATA_HASHING` | true | false | false | false |
| `use_igvf_reference` | **false** (custom Gencode v46) | true | true | true |
| `is_10x3v3` | false | **true** | **true** | false |
| `reverse_complement_guides` | true | false | false | true |
| `spacer_tag` | `""` (no spacer) | `GAGTACATGGGGG` | `GAGTACATGGGGG` | `TAGCTCTTAAAC` |
| `QC_min_genes_per_cell` | 500 | 500 | 500 | **2500** |
| `QC_pct_mito` | 15 | 20 | 20 | 15 |
| `QC_barcode_filter` | `knee2` | `knee2` | `knee2` | **`none`** |
| `GUIDE_ASSIGNMENT_method` | sceptre | sceptre | sceptre | **cleanser** |
| `GUIDE_ASSIGNMENT_capture_method` | CROP-seq | CROP-seq | CROP-seq | **direct-capture** |
| `INFERENCE_SCEPTRE_control_group` | default | default | default | complement |
| `INFERENCE_SCEPTRE_GENE_CHUNK_SIZE` | 1000 | 1000 | 1000 | 1000 |
| Base container | `pinellolab/crispr_pipeline/conda-docker:latest` | `sjiang9/conda-docker:0.3` | `sjiang9/conda-docker:0.3` | `pinellolab/crispr_pipeline/conda-docker:latest` (sif) |
| `max_cpus` / `max_memory` | 16 / 1500 GB | 128 / 256 GB | 128 / 256 GB | 12 / 600 GB |

Common across all four runs: `ENABLE_SCRUBLET=false`, `DUAL_GUIDE=false`, `QC_min_cells_per_gene=0.05`, `Multiplicity_of_infection=high`, `INFERENCE_max_target_distance_bp=1000000`, `INFERENCE_SCEPTRE_side=both`, `INFERENCE_SCEPTRE_grna_integration_strategy=union`, `INFERENCE_SCEPTRE_resampling_approximation=skew_normal`, sceptre `sjiang9/sceptre-igvf:0.1`, perturbo `ghcr.io/pinellolab/perturbo:sha-f3dc8ca`, cleanser `gersbachlab-bioinformatics/cleanser:1.2.1`.

## Inputs

Each onboarded dataset has its QC outputs at `datasets/<dataset>/<run>/qc/`, with the standard three subdirectories produced by `scripts/qc_array.sh`:

```
datasets/<dataset>/<run>/qc/
  mapping_gene/
    gene_metrics.tsv          # per-batch gene-expression QC
  mapping_guide/
    guide_metrics.tsv         # per-batch guide-capture QC
    per_guide_capture.tsv     # per-guide cells/UMI table
  intended_target/
    intended_target_metrics.tsv
    intended_target_results.tsv
```

The benchmark drives its notebooks from a manifest TSV (`manifests/<run-label>_qc_paths.tsv`, see [docs/manuscripts/CRISPRi_tech_benchmark/manifests/](../../../../../manuscripts/CRISPRi_tech_benchmark/manifests/)). The same pattern will be adopted here: one row per dataset with columns `dataset | qc_dir | gene_metrics | guide_metrics | per_guide_capture | intended_target_results | intended_target_metrics`. The manifest itself is **deferred to a follow-up iteration** — this README defines the contract.

## Planned notebooks

Three notebooks, ported from [bin/2_qc/](../../../../../manuscripts/CRISPRi_tech_benchmark/bin/2_qc/), driven by the manifest above. All outputs land in `results/cross_production_qc/`.

### 1. `1_metrics_comparison.ipynb` — gene + guide metrics

Aggregate-level barplots and per-lane breakdowns. Mirrors `gene_metrics_comparison.pdf`, `guide_metrics_comparison.pdf`, `metrics_by_lane_comparison.pdf` from the benchmark.

Panels:
- **Number of cells** — total cells passing QC per dataset
- **Median Tx UMIs / cell** — sequencing depth
- **Median genes / cell** — library complexity
- **Median Tx %MT UMIs** — cell quality indicator
- **Median sgRNA UMIs / cell** — guide library depth
- **Median guides / cell** — guide MOI
- **Fraction cells with sgRNA UMIs > 0** — guide detection rate
- **Median cells / guide** — coverage uniformity
- **Mean assigned sgRNA / cell** — assignment-level MOI
- **Cells with exactly 1 assigned sgRNA** — singlet fraction
- Per-lane (per-batch) versions of the six core metrics above

## Output directory

```
results/cross_production_qc/
  *.pdf            # plots
  *.tsv / *.tsv.gz # intermediate combined tables
```

## Color scheme

Lineage-themed palette in [config/colors/production_TF-Perturb-seq.yaml](../../../../../../config/colors/production_TF-Perturb-seq.yaml):

| Dataset | Lineage cue | Hex |
|---|---|---|
| Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq | cardiac red | `#D7263D` |
| Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq | endoderm gold | `#F6BD16` |
| Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq | stem-cell teal | `#1B998B` |
| Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq | royal purple | `#7C3F98` |
| Engreitz_WTC11-endothelial-cells_TF-Perturb-seq | vascular pink | `#E36A9E` |

Load via the existing config loader ([config/loader.py](../../../../../../config/loader.py)):

```python
from tf_perturb_seq.config.loader import load_colors
colors = load_colors("production_TF-Perturb-seq", "dataset_colors")
order  = load_colors("production_TF-Perturb-seq", "dataset_order")
```

## Status / outstanding work

- **Engreitz endothelial** is not yet onboarded — no `datasets/Engreitz_WTC11-endothelial-cells_TF-Perturb-seq/` directory exists. Notebooks should skip it gracefully until it lands.
- **Gersbach hepatocyte** only has QC re-run locally over Sara's h5mu; there is no full local CRISPR pipeline run, and `per_guide_capture.tsv` may be absent — confirm before running notebook 3.
- **Manifest TSV** (`manifests/<run-label>_qc_paths.tsv`) is not yet written; the next iteration ports it from the benchmark conventions.
- **Notebooks** (1/2/3 above) are not yet written; they will be ported from [bin/2_qc/](../../../../../manuscripts/CRISPRi_tech_benchmark/bin/2_qc/).

## Synapse outputs

Final figures and combined tables should be mirrored to WG1's Synapse folder [`syn74954079`](https://www.synapse.org/Synapse:syn74954079) under `working_groups/wg1_data_qc/pipeline_qc/`, per the [WG1 README](../README.md) sharing flow.
