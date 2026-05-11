# Energy Distance Pipeline Output Directory Structure

Reference output paths:

- **Production-pattern run (verified)**: [`syn74381167`](https://www.synapse.org/Synapse:syn74381167) (`edist_analysis_result` under `Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2`). Steps 1 + 2 + 2.1 only — Step 3 not yet run.
- Pipeline-author example: `tf_perturb_seq/external/energy_dist_pipeline/pipeline_output/`
- Earlier HPC run: `/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/results/energy_distance/2026_04_08_prelim/`
- Production-dataset runs: not yet generated as of 2026-05-07.

> **Verification status — VERIFIED for Steps 1 + 2 + 2.1 against `syn74381167`.** File names + column structure of `pval_edist_full.csv`, `targeting_outlier_table.csv`, `non_targeting_outlier_table.csv`, and the four `image/*.pdf` figures are confirmed. Step 3 outputs (`target_by_target_matrix.csv`, `edist_embedding_info.csv`) and the optional `discordance_gRNA.csv` / `annotation_file_table.csv` were not present in that run — schemas for those remain sampled.

## Pipeline stages and their outputs

The pipeline is run as four scripts in sequence (see `docs/analysis/ENERGY_DISTANCE.md` for the full guide). Each stage is responsible for a different set of output files inside `OUTPUT_FOLDER`.

### 0. Preprocessing

Either `bin/0_preprocess.py` (legacy) or `preprocess_mudata.py` from the [`Chikara-Takeuchi/energy_dist_TFperturb`](https://github.com/Chikara-Takeuchi/energy_dist_TFperturb) wrapper.

| File | Contents | Useful for |
|------|----------|------------|
| `<OUTPUT_FOLDER>/annotation_file_table.csv` | Per-gRNA annotation following the IGVF DACC `guide_rna_sequences` schema. Columns: `guide_id`, `intended_target_name`, `type`, `spacer`. | Used internally by all downstream steps; same content as the canonical guide library. |
| `<OUTPUT_FOLDER>/pca_dataframe.pickle` | Cell × PC matrix extracted from `adata.obsm["X_pca"]`. | Cached so subsequent steps don't re-derive PCA. |
| `<OUTPUT_FOLDER>/gRNA_dictionary.pickle` | Dict mapping `gRNA_name → list of cell barcodes`. | Cached so subsequent steps don't re-build the gRNA→cell map. |

### 1. gRNA filtering — `bin/1_filtereing_gRNA.py`

| File | Contents | Useful for |
|------|----------|------------|
| `<OUTPUT_FOLDER>/targeting_outlier_table.csv` | Per-targeting-gRNA outlier statistics (DISCO test + hypergeometric ranks against sibling gRNAs targeting the same region). Index = `gRNA_name`. Column: `pval_outlier`. | Identifies outlier gRNAs that disagree with siblings; downstream steps drop these. |
| `<OUTPUT_FOLDER>/non_targeting_outlier_table.csv` | Per-non-targeting-gRNA outlier flag (K-means on pairwise energy-distance matrix between non-targeting gRNAs). Index = `gRNA_name`. Column: `pval_outlier`. | Identifies non-targeting gRNAs that don't behave like the rest; treated as non-baseline. |
| `<OUTPUT_FOLDER>/discordance_gRNA.csv` _(legacy / not exported by the new wrapper)_ | Per-gRNA discordance status — flags guides whose effect disagrees with siblings for the same target. Columns: `gRNA_name`, `status`. | QC; understanding which guides are reliable. |

### 2. Energy distance vs. non-targeting — `bin/2_e_distance_nontargeting.py`

| File | Contents | Useful for |
|------|----------|------------|
| `<OUTPUT_FOLDER>/pval_edist_full.csv` | **Primary results.** One row per target region. Index = target id in format `<ENSG>\|<chr>:<start>-<end>` (built from the gRNA annotation by `preprocess_mudata.py`). Columns: `cell_count`, `type`, `distance_0..distance_19`, `pval_0..pval_19`, `distance_mean`, `pval_mean`, `pval_mean_log`, `distance_mean_log`. In the new wrapper, the `type` column only takes the value `targeting`; positive/negative-control attribution lives in the input annotation, not this column. | **Significance + magnitude per TF.** Working Group 1's main input for "TFs that significantly alter the transcriptome". |
| `<OUTPUT_FOLDER>/figures/gRNA_stat.pdf` | Histogram of gRNAs per target — dropped vs. outlier vs. retained. | QC. |
| `<OUTPUT_FOLDER>/figures/e-dist_distribution.pdf` | Distribution of energy distance values vs. p-values (scatter + density). | QC; cutoff selection. |
| `<OUTPUT_FOLDER>/figures/e-dist_cutoff_value.pdf` | Heatmap of #significant targets across (p-value cutoff × energy-distance cutoff) grids. | Aiding cutoff choice for Step 3. |
| `<OUTPUT_FOLDER>/figures/e-dist_cutoff_value_NEG_CONTROL.pdf` | Same heatmap restricted to negative controls. | Cutoff calibration. |

### 2.1. Diagnostic plots — `bin/2_1_Plot_figure.py`

This step is plot-only; figures land in `<OUTPUT_FOLDER>/figures/` (the four PDFs above are owned by step 2.1, even though they're listed under step 2 because they describe step-2 results).

### 3. Energy distance among regions — `bin/3_e_distance_among_regions.py`

Run separately (`run_step3.sh` / `run_container_step3.sh`) after picking cutoffs in `config_clustering.json`. Operates **only on significant target regions**.

| File | Contents | Useful for |
|------|----------|------------|
| `<OUTPUT_FOLDER>/target_by_target_matrix.csv` | Square pairwise energy-distance matrix between all significant target regions. First column = target id; remaining columns = one per target. | Working Group 1: clustering / UpSet of TFs with shared effects across lineages. |
| `<OUTPUT_FOLDER>/edist_embedding_info.csv` | 2D t-SNE embedding of significant targets (computed from the distance matrix), with cluster labels (Affinity Propagation by default). Columns: `index` (target id), `x`, `y`, `cluster`. | Visualization of TF perturbation similarity. |

## Marker conventions in `pval_edist_full.csv`

- `type=target` — guides target a real region (TF gene promoter)
- `type=positive control` — guides target a positive control gene (e.g., AARS — should always show a strong effect)
- `type=negative control` — guides target a negative control region (should be near baseline)
- Non-targeting guides do **not** get rows in this table; they're the implicit baseline.

## Per-dataset run status (2026-05-08)

| Dataset | Run? | Run path / notes |
|---------|:---:|------|
| Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq | ⏳ | submission ready: `datasets/Hon_WTC11-.../5_run_energy_distance.sh` (sources MuData from Synapse [`syn74522725`](https://www.synapse.org/Synapse:syn74522725) — Hon's GCS run is incomplete as of today). |
| Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq | ⏳ | submission ready: `datasets/Huangfu_..._definitive-endoderm.../5_run_energy_distance.sh` (sources MuData from `gs://.../muddy_penguin/inference_mudata.h5mu`). |
| Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq | ⏳ | submission ready: `datasets/Huangfu_..._embryonic-stemcell.../5_run_energy_distance.sh` (sources MuData from `gs://.../sceptre_v1/inference_mudata.h5mu`). |
| Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq | ☑ | already run by Sara Geraghty (Duke). Output layout we'll match: [`syn74381167`](https://www.synapse.org/Synapse:syn74381167) (Gersbach HTv2 benchmark — same group's tooling). Production hepatocyte output to be supplied by Sara. |
| Engreitz_WTC11-endothelial-cells_TF-Perturb-seq | ☐ | gated on inference MuData (no portal data yet). |
| Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2 (benchmark, reference only) | ☑ | [`syn74381167`](https://www.synapse.org/Synapse:syn74381167) — used to verify the schema. |

## Synapse layout

When production runs land, mirror to:
`syn64423137/2026_UTSW/datasets/<dataset_id>/energy_distance/`

A placeholder folder already exists at [`syn74823521`](https://www.synapse.org/Synapse:syn74823521) (`edist_result`) under the project root.

## TODOs

- [ ] Verify column lists against a real production output (numbers will scale; check for any new columns).
- [ ] Build a cross-dataset summary TSV: per-dataset count of TFs significant at FDR thresholds (rolled up from each dataset's `pval_edist_full.csv`).
- [ ] Build a cross-dataset clustered heatmap of per-target energy distances (joins each dataset's `target_by_target_matrix.csv` on shared targets).
