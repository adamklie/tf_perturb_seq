# Data

## Dataset naming convention

`<Lab>_<CellLine>-<condition>_TF-Perturb-seq`

All datasets live under `datasets/`. Each follows a numbered script convention:
1. `1_generate_per_sample_metadata.sh` — Query IGVF portal for metadata
2. `2_upload_to_gcp.sh` — Upload fastqs to GCS
3. `3_run_CRISPR_pipeline.sh` (or `4_*`) — Run the CRISPR analysis pipeline
4. `5_run_energy_distance.sh` — Energy distance analysis (**TODO**)
5. `6_run_cnmf.sh` — Gene program discovery (**TODO**)

## IGVF Portal Structure

The CRISPR FG pipeline requires an **analysis set** on the portal to generate the sample metadata CSV. An analysis set groups measurement sets (scRNA) and auxiliary sets (gRNA, HTO) together.

**What a complete portal setup looks like** (see Huangfu WTC11 benchmark as reference):
- **Analysis Set** (e.g., `IGVFDS8556LYRW`) — groups all file sets for the dataset
  - **Measurement Sets** — one per 10x lane/pool (scRNA-seq data)
    - Must have: `strand_specificity`, `onlist_files`, `onlist_method`
    - Contains R1/R2 sequence files with seqspecs
  - **Auxiliary Sets** — one per measurement set (gRNA sequencing)
    - Must link back to its measurement set via `measurement_sets`
    - Type: `gRNA sequencing`
  - **Construct Library Set** — exactly one, with `guide_design` file
    - Must have one `integrated_content_files` entry with `content_type: "guide RNA sequences"`

**Portal search for TF Perturb-seq datasets:**
https://data.igvf.org/search/?type=MeasurementSet&preferred_assay_titles=Perturb-seq&collections=TF+Perturb-seq+Project

## Technology benchmark datasets

Cross-lab comparison using the same WTC11 iPSC line and shared ~50-gene guide library.

| Dataset | Lab | Technology | Cell Line | Portal Accession | Pipeline Runs | Status |
|---------|-----|------------|-----------|------------------|---------------|--------|
| [Hon_WTC11-benchmark](../datasets/Hon_WTC11-benchmark_TF-Perturb-seq/) | Hon | 10x 5' HT v2 w/ HTO | WTC11 | [IGVFDS4761PYUO](https://data.igvf.org/analysis-sets/IGVFDS4761PYUO) | cleanser_unified, sceptre_v11 | Complete (stages 1-7) |
| [Huangfu_WTC11-benchmark](../datasets/Huangfu_WTC11-benchmark_TF-Perturb-seq/) | Huangfu | 10x 3' v3 | WTC11 | [IGVFDS8556LYRW](https://data.igvf.org/analysis-sets/IGVFDS8556LYRW) | cleanser_unified, sceptre_v11 | Complete (stages 1-6) |
| [Engreitz_WTC11-benchmark](../datasets/Engreitz_WTC11-benchmark_TF-Perturb-seq/) | Engreitz | CC Perturb-seq | WTC11 | [IGVFDS5057HJKP](https://data.igvf.org/analysis-sets/IGVFDS5057HJKP) | cleanser_unified, sceptre_v11 | Complete (stages 1-4) |
| [Gersbach_WTC11-benchmark_GEM-Xv3](../datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/) | Gersbach | 10x GEM-X 3' | WTC11 | [IGVFDS6673ZFFG](https://data.igvf.org/analysis-sets/IGVFDS6673ZFFG) | cleanser_unified, sceptre_v11 | Complete (stages 1-6) |
| [Gersbach_WTC11-benchmark_HTv2](../datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/) | Gersbach | 10x HT v2 | WTC11 | [IGVFDS6237URFJ](https://data.igvf.org/analysis-sets/IGVFDS6237URFJ) | cleanser_unified, sceptre_v11 | Complete (stages 1-6) |

Cross-dataset benchmark analyses live in: [technology-benchmark_WTC11_TF-Perturb-seq/](../datasets/technology-benchmark_WTC11_TF-Perturb-seq/)

## Production datasets

Full TF library (~2000 targets) applied to differentiated lineages.

| Dataset | Lab | Cell Line | Lineage | Portal Status | Pipeline Status |
|---------|-----|-----------|---------|---------------|-----------------|
| [Hon_WTC11-cardiomyocyte-differentiation](../datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/) | Hon | WTC11 | Cardiomyocyte (12-day) | ? | CRISPR pipeline troubleshooting (Weizhou) |
| [Hon_WTC11-cardiomyocyte-differentiation_Pilot](../datasets/Hon_WTC11-cardiomyocyte-differentiation_Pilot-TF-Perturb-seq/) | Hon | WTC11 | Cardiomyocyte (pilot) | [IGVFDS4389OUWU](https://data.igvf.org/analysis-sets/IGVFDS4389OUWU) | Early stage |
| [Huangfu_HUES8-definitive-endoderm](../datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/) | Huangfu | HUES8 | Definitive endoderm | **BLOCKED** (see below) | Scripts scaffolded, waiting on portal |
| [Huangfu_HUES8-embryonic-stemcell](../datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/) | Huangfu | HUES8 | Embryonic stem cell | **BLOCKED** (see below) | Scripts scaffolded, waiting on portal |
| [Gersbach_WTC11-hepatocyte](../datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/) | Gersbach | WTC11 | Hepatocyte | ? | CRISPR pipeline troubleshooting (Sara) |

### Huangfu HUES8 portal status (audited 2026-03-25)

Both datasets have measurement sets, auxiliary sets, and a construct library set on the portal. The **only missing piece** is the analysis set that groups them together.

**What exists:**

| | Measurement Sets | Auxiliary Sets (gRNA) | Construct Library Set | Analysis Set | Seqspecs |
|-|:---:|:---:|:---:|:---:|:---:|
| HUES8 Stemcell | 8 | 8 | [IGVFDS3299AXST](https://data.igvf.org/construct-library-sets/IGVFDS3299AXST) (guide: [IGVFFI8270UPKB](https://data.igvf.org/tabular-files/IGVFFI8270UPKB)) | **None** | **Missing** on files |
| HUES8 Endoderm | 8 | 8 | same as above | **None** | **Missing** on files |
| WTC11 Benchmark (for reference) | 4 | 4 | separate CLS | [IGVFDS8556LYRW](https://data.igvf.org/analysis-sets/IGVFDS8556LYRW) | Present |

**Measurement set metadata:** Looks good — `strand_specificity`, `onlist_files`, `onlist_method` all present on measurement sets.

**HUES8 Stemcell measurement sets -> auxiliary sets:**
| Measurement Set | Auxiliary Set (gRNA) |
|---|---|
| IGVFDS0746MYRH | IGVFDS8263ODBC |
| IGVFDS0956UUQN | IGVFDS2137GIIX |
| IGVFDS2908DIHX | IGVFDS2178ZGIR |
| IGVFDS2940LYGK | IGVFDS3927XOZY |
| IGVFDS3348KRMQ | IGVFDS0173JYTZ |
| IGVFDS5934EZKT | IGVFDS9550TKIZ |
| IGVFDS6247FKDR | IGVFDS6962YCLR |
| IGVFDS8623ONYE | IGVFDS4608EUUL |

**HUES8 Endoderm measurement sets -> auxiliary sets:**
| Measurement Set | Auxiliary Set (gRNA) |
|---|---|
| IGVFDS0788TUIL | IGVFDS7673CCYF |
| IGVFDS1313VGHL | IGVFDS5105RAVI |
| IGVFDS1403KZYV | IGVFDS0367BQIW |
| IGVFDS1437XFJK | IGVFDS6182MVVF |
| IGVFDS2129VHBD | IGVFDS2661AFQU |
| IGVFDS2520ZEYH | IGVFDS7765TYSM |
| IGVFDS3434HGEI | IGVFDS5564SHUB |
| IGVFDS4079VXPX | IGVFDS1454WEUZ |

**ACTION NEEDED** (coordinate with Huangfu lab / Denis Torre / DACC):
1. **Create analysis sets** — one for stemcell (grouping 8 MS + 8 aux), one for endoderm (grouping 8 MS + 8 aux). This is the primary blocker for running `generate_per_sample.py`.
2. **Add seqspecs** to sequence files — R1 files currently have no seqspecs linked. Either upload seqspecs to the portal or provide fallback YAML paths at runtime (like the benchmark did).
3. Note: auxiliary sets are not tagged with `collections: TF Perturb-seq Project` — may want to fix for discoverability.

Legacy internal processing outputs have been archived to `scratch/2026_03_25_legacy_archives/`.

## Public datasets

| Dataset | Description |
|---------|-------------|
| [ENCODE_WTC11_ChIP-seq](../datasets/ENCODE_WTC11_ChIP-seq/) | 90 ENCODE TF ChIP-seq experiments in WTC11; used for validation |

## Shared data locations

- **Synapse**: https://www.synapse.org/Synapse:syn63675917
- **GCP**: `igvf-pertub-seq-pipeline` project on Google Cloud
- **IGVF Portal**: https://data.igvf.org/
- **Guide metadata on portal**: https://data.igvf.org/tabular-files/IGVFFI5765HMZH/

## Reference files (ref/)

Located in `ref/` at the repo root:
- `finalized_annotation_files/` — Guide-to-target mapping files used by pipelines
- GTF genome annotation files
- Control guide lists
- Gene sets and target lists

## Pipeline run outputs

For benchmark datasets, pipeline outputs are stored in `datasets/<name>/runs/` with subdirectories per run configuration (e.g., `cleanser_unified/`, `sceptre_v11/`). The primary output is an `inference_mudata.h5mu` MuData object containing gene expression, guide assignments, and inference results.
