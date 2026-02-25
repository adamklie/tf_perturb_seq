# TF-Perturb-seq Technology Benchmark

Master tracking document for running the CRISPR Pipeline on WTC11 benchmark datasets.

## Known Issues

- Engreitz seqspecs not yet available on portal
- **Seqspec YAML files from portal are gzipped** - Pipeline expects uncompressed `.yaml` files but portal provides `.yaml.gz`. Must decompress to `patch/` directory along with barcode_onlist and guide_design files.
- **Huangfu barcode_onlist mismatch** - The onlist associated with the Huangfu benchmark datasets (`IGVFFI9487JPEN`) doesn't match what's in the seqspec. Based on Gary's previous pipeline runs, the proper onlist is `IGVFFI4695IKAL` (https://data.igvf.org/tabular-files/IGVFFI4695IKAL/). Fixed in patched_v2.
- **Pipeline `seqSpecParser` gzip bug** - The pipeline copies barcode files as-is (`cp "${barcode_file}" barcode_file.txt`) without decompressing, so `kb count -w barcode_file.txt` fails on gzipped content. Workaround: pre-decompress to `patch/` directory.
- **Missing `num_expressed_genes` in 3/5 Synapse MuDatas** — see Issue 1 below.
- **Engreitz batch ID mismatch between gene and guide modalities** — see Issue 2 below.

### Issue 1: Missing `num_expressed_genes` column in Synapse-downloaded MuDatas

**Affected datasets:** Gersbach GEM-Xv3 (`jade_dolphin`), Gersbach HTv2 (`amber_sparrow`), Engreitz (`cobalt_heron`)
**Unaffected datasets:** Hon (`golden_falcon`), Huangfu (`silver_otter`)

**Symptom:** `gene_metrics_comparison.pdf` and `metrics_by_lane_comparison.pdf` are missing the "Median Genes" panel for the 3 affected datasets. The `genes_median` through `genes_q75` columns are empty in `combined_gene_metrics.tsv`.

**Root cause:** The `num_expressed_genes` column is absent from `gene.obs` in the inference MuDatas downloaded from Synapse for Gersbach GEM-Xv3, Gersbach HTv2, and Engreitz. The QC script (`src/tf_perturb_seq/qc/mapping_gene.py`, line 84) checks for this column and silently skips all `genes_*` metrics when it's missing.

**Columns present in ALL 5 datasets' `gene.obs`:**
`batch`, `cov1`, `batch_number`, `n_counts`, `log1p_n_genes_by_counts`, `total_gene_umis`, `log1p_total_counts`, `pct_counts_in_top_50_genes`, `pct_counts_in_top_100_genes`, `pct_counts_in_top_200_genes`, `pct_counts_in_top_500_genes`, `total_counts_mt`, `log1p_total_counts_mt`, `percent_mito`, `total_counts_ribo`, `log1p_total_counts_ribo`, `pct_counts_ribo`

**Column present ONLY in Hon + Huangfu:** `num_expressed_genes`

**Why the difference:** Hon and Huangfu MuDatas were produced by pipeline runs that computed `num_expressed_genes` (likely a newer pipeline version or different config). The Gersbach and Engreitz MuDatas were uploaded to Synapse by other groups (Alex Barrera, Engreitz lab) using runs that did not include this column.

**Potential fix:** The `log1p_n_genes_by_counts` column IS present in all 5 datasets. This is `log(1 + n_genes_by_counts)` from scanpy's `sc.pp.calculate_qc_metrics()`, where `n_genes_by_counts` is the number of genes with at least 1 count — equivalent to `num_expressed_genes`. The value can be recovered via `np.expm1(log1p_n_genes_by_counts)`. Options:
1. Update `mapping_gene.py` to fall back to deriving from `log1p_n_genes_by_counts` when `num_expressed_genes` is missing
2. Materialize the column into `.obs` on the affected MuDatas before re-running QC
3. Add `num_expressed_genes` to the plotting functions as well (they currently skip histograms entirely if the column is absent)

After fixing, re-run the QC pipeline for the 3 affected datasets and regenerate the cross-tech comparison notebook.

### Issue 2: Engreitz batch ID mismatch between gene and guide modalities

**Affected dataset:** Engreitz (`cobalt_heron`) only

**Symptom:** `metrics_by_lane_comparison.pdf` is missing `guide_umi_median` and `guides_per_cell_mean` for all 15 Engreitz lanes, even though the per-lane guide metrics TSV (`Engreitz_WTC11-benchmark_TF-Perturb-seq_guide_metrics.tsv`) is fully populated.

**Root cause:** The gene and guide modalities in the Engreitz MuData use completely different analysis set accessions as batch labels. There is zero overlap:
- `gene.obs['batch']`: `IGVFDS0208BWWR_1`, `IGVFDS0208BWWR_2`, `IGVFDS1283GTHY_1`, ... (15 batches)
- `guide.obs['batch']`: `IGVFDS2013SMCO_1`, `IGVFDS2013SMCO_2`, `IGVFDS2018NMHT_1`, ... (15 batches)
- Intersection: empty set

The cross-tech comparison notebook (`1_metrics_comparison.ipynb`) joins gene and guide per-lane metrics on the `batch` column. Since no Engreitz batch labels match across modalities, the guide columns come back as NaN for all Engreitz lanes.

**Why the difference:** The other 4 datasets use the same analysis set accession for both gene and guide modalities (i.e., the RNA and guide libraries share a batch label). Engreitz appears to use separate analysis sets for RNA vs guide libraries, resulting in different `batch` values in each modality. This is likely an artifact of how the Engreitz pipeline was configured (CC Perturb-seq uses a different library structure than 10x-based assays).

**Potential fix:** Options:
1. Build a mapping between Engreitz gene batch IDs and guide batch IDs (they should correspond 1:1 based on the same underlying sample) and harmonize before running QC
2. Modify the comparison notebook to handle mismatched batch IDs by falling back to overall (non-per-lane) guide metrics for datasets where per-lane join fails
3. Ask the Engreitz group to confirm the correct batch pairing

## 2026-02-15: New seqspec update

Replaced all seqspec paths across the 4 benchmark datasets with new shared seqspecs from Lucas:
`gs://igvf-pertub-seq-pipeline-data/scratch/bioinfolucas/new_03_seqspec/`

| Dataset | Source samplesheet | scRNA seqspec | gRNA seqspec | hash seqspec |
|---------|-------------------|---------------|--------------|--------------|
| Gersbach GEM-Xv3 | `sample_metadata_gcp_2026_02_04_patched.csv` | `alex_IGVFDS6237URFJ_rna_10x_5_gex.yaml` | `alex_IGVFDS6237URFJ_guide_10x_5_gex.yaml` | - |
| Gersbach HTv2 | `sample_metadata_gcp_2026_02_04_patched.csv` | `alex_IGVFDS6673ZFFG_rna_chomium5.yaml` | `alex_IGVFDS6673ZFFG_guide_HTV2_chomium5.yaml` | - |
| Hon | `sample_metadata_gcp_2026_01_18_patched.csv` | `gary_IGVFDS4761PYUO_rna.yaml` | `gary_IGVFDS4761PYUO_guide.yaml` | `gary_IGVFDS4761PYUO_hash.yaml` |
| Huangfu | `sample_metadata_gcp_2026_01_30_patched_v2.csv` | `huangfu_IGVFDS1889TBEY_rna.yaml` | `huangfu_IGVFDS1889TBEY_guide.yaml` | - |

New samplesheets saved as `sample_metadata_gcp_2026_02_15.csv` in each dataset directory (local and GCS).

## Previous patching notes (Hon sample sheet, 2026-01-18)

The following summarizes what was done to the Hon_WTC11-benchmark_TF-Perturb-seq sample sheet to get the pipeline to run successfully on GCP on 2026-01-18.

| Change | Original Value | Patched Value | Reason |
|--------|----------------|---------------|--------|
| **Hash file path** (`barcode_hashtag_map`) | `gs://.../2024/12/20/73fc75aa-.../IGVFFI3028BJFI.tsv.gz` | `gs://.../patch/IGVFFI3028BJFI_with_header.tsv.gz` | Original file missing required `hash_id\tsequence` header |
| **Barcode onlist path** (`barcode_onlist`) | `gs://.../2024/12/05/996cf65e-.../IGVFFI9487JPEN.tsv.gz` | `gs://.../patch/IGVFFI9487JPEN.tsv` | Pipeline renames file to `.txt` without decompressing; `kb count -w` fails on gzipped content |
| **Guide design path** (`guide_design`) | `gs://.../2025/08/04/cc3d00eb-.../IGVFFI5765HMZH.tsv.gz` | `gs://.../patch/IGVFFI5765HMZH.tsv` | Precautionary unzip (see note below) |
| **Lane column** | `3`, `5`, `7`, `8` (various per row) | `` (empty for all rows) | Match working sample sheet structure |
| **Sequencing run column** | `1` | `1` (unchanged) | Already matched working sample sheet |

## Dataset Overview

| Dataset | Lab | Technology | Status |
|---------|-----|------------|---------|
| [Engreitz](../Engreitz_WTC11-benchmark_TF-Perturb-seq/README.md) | Engreitz | CC Perturb-seq | inference_mudata downloaded from Synapse (2026-02-17); QC pipeline run |
| [Gersbach GEM-Xv3](../Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/README.md) | Gersbach | 10x GEM-X 3' | New seqspecs uploaded (2026-02-15), ready for pipeline run |
| [Gersbach HTv2](../Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/README.md) | Gersbach | 10x HT v2 | New seqspecs uploaded (2026-02-15), ready for pipeline run |
| [Hon](../Hon_WTC11-benchmark_TF-Perturb-seq/README.md) | Hon | 10x 5' HT v2 w/ HTO | Pipeline run completed; new seqspecs uploaded (2026-02-15) |
| [Huangfu](../Huangfu_WTC11-benchmark_TF-Perturb-seq/README.md) | Huangfu | 10x 3' v3 | New seqspecs uploaded (2026-02-15), ready for pipeline run |
---

## Engreitz_WTC11-benchmark_TF-Perturb-seq

- **Analysis Set Accession**: `IGVFDS6673ZFFG`
- **Sample Metadata CSV**: Not yet generated
- **Guide metadata**: `IGVFFI5765HMZH`
- **Seqspecs**: Validating on portal
- **Barcode Onlist**: TODO
- **Status**: Portal missing seqspecs, need to follow up with Engreitz lab
- **Key pipeline parameters:**
   - `ENABLE_DATA_HASHING = false`
   - `Multiplicity_of_infection = 'high'`

## Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3

- **Analysis Set Accession**: `IGVFDS6237URFJ`
- **Sample Metadata CSV**: `sample_metadata_gcp_2026_02_15.csv` (based on patched version)
- **Guide Metadata**: `IGVFFI5765HMZH` (decompressed in `patch/`)
- **Barcode Onlist**: `IGVFFI1697XAXI` (decompressed in `patch/`)
- **Seqspecs**: `alex_IGVFDS6237URFJ_rna_10x_5_gex.yaml` / `alex_IGVFDS6237URFJ_guide_10x_5_gex.yaml`
- **Status**: New seqspecs uploaded (2026-02-15), ready for pipeline run
- **Key Parameters:**
   - `ENABLE_DATA_HASHING = false`

## Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2

- **Analysis Set Accession**: `IGVFDS6673ZFFG`
- **Sample Metadata CSV**: `sample_metadata_gcp_2026_02_15.csv` (based on patched version)
- **Guide Metadata**: `IGVFFI5765HMZH` (decompressed in `patch/`)
- **Barcode Onlist**: `IGVFFI9487JPEN` (decompressed in `patch/`)
- **Seqspecs**: `alex_IGVFDS6673ZFFG_rna_chomium5.yaml` / `alex_IGVFDS6673ZFFG_guide_HTV2_chomium5.yaml`
- **Status**: New seqspecs uploaded (2026-02-15), ready for pipeline run
- **Key Parameters:**
   - `ENABLE_DATA_HASHING = false`

## Hon_WTC11-benchmark_TF-Perturb-seq

- **Analysis Set Accession**: `IGVFDS4761PYUO`
- **Sample Metadata CSV**: `sample_metadata_gcp_2026_02_15.csv` (based on patched version)
- **Guide Metadata**: `IGVFFI5765HMZH` (decompressed in `patch/`)
- **Seqspecs**: `gary_IGVFDS4761PYUO_rna.yaml` / `gary_IGVFDS4761PYUO_guide.yaml` / `gary_IGVFDS4761PYUO_hash.yaml`
- **Barcode Onlist**: `IGVFFI9487JPEN` (decompressed in `patch/`)
- **Hashing Map**: Patched `IGVFFI3028BJFI` with header
- **Status**: Pipeline successfully run to completion; new seqspecs uploaded (2026-02-15)
- **Key Parameters:**
   - `ENABLE_DATA_HASHING = true`
   - `Multiplicity_of_infection = 'high'`
- **Successful Runs (see configs)**
   - harmless_teacher: initial run with hashing
   - beautiful_tornado: harmless teacher but with cleanser instead of sceptre
   - introverted_frankenstein: disabled hashing, back to sceptre

## Huangfu_WTC11-benchmark_TF-Perturb-seq

- **Analysis Set Accession**: `IGVFDS1889TBEY`
- **Sample Metadata CSV**: `sample_metadata_gcp_2026_02_15.csv` (based on patched_v2 version)
- **Guide Metadata**: `IGVFFI5765HMZH` (decompressed in `patch/`)
- **Barcode Onlist**: `IGVFFI4695IKAL` (decompressed in `patch/`)
- **Seqspecs**: `huangfu_IGVFDS1889TBEY_rna.yaml` / `huangfu_IGVFDS1889TBEY_guide.yaml`
- **Modalities**: scRNA, gRNA (no cell hashing)
- **Status**: New seqspecs uploaded (2026-02-15), ready for pipeline run
- **Key Parameters:**
   - `ENABLE_DATA_HASHING = false`
   - `DUAL_GUIDE = false`
   - `is_10x3v3 = true` (**unique to this dataset!**)
   - `Multiplicity_of_infection = 'high'`
   - `GUIDE_ASSIGNMENT_method = 'sceptre'`
   - `GUIDE_ASSIGNMENT_capture_method = 'CROP-seq'`
   - `use_igvf_reference = true`


Hon_WTC11-benchmark_TF-Perturb-seq
Benchmark
syn73743227

Huangfu_WTC11-benchmark_TF-Perturb-seq
Benchmark
syn72386406

Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2
Benchmark
syn73712262

Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3
Benchmark
syn73712373

Engreitz_WTC11-benchmark_TF-Perturb-seq
Benchmark
syn73615563
