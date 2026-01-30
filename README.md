# IGVF TF Perturb-seq Pillar Project

Code for processing and analysis of CRISPR screening data targeting a large set of known transcription factors and epigenetic modifiers across multiple cell states.

## Quick Start

```bash
# Clone the repository
git clone https://github.com/adamklie/tf_perturb_seq.git
cd tf_perturb_seq

# Install dependencies (using uv)
uv sync
```

## Repository Structure

```
tf_perturb_seq/
├── datasets/                       # Dataset directories (one per experiment)
│   ├── <dataset_name>/             # Individual dataset
│   │   ├── README.md
│   │   ├── 1_generate_per_sample_metadata.sh
│   │   ├── 2_upload_to_gcp.sh
│   │   ├── 3_run_CRISPR_pipeline.sh
│   │   ├── patch_files.sh          # Decompress gzipped files
│   │   ├── sample_metadata.csv
│   │   └── sample_metadata_gcp_*.csv
│   └── technology-benchmark_WTC11_TF-Perturb-seq/
│       ├── config/                 # Analysis configuration
│       ├── bin/                    # Cross-dataset analysis scripts
│       └── results/                # Analysis outputs
├── docs/                           # Documentation
│   ├── CRISPR_PIPELINE.md          # Pipeline workflow guide
│   └── dataset_template/           # Template for new datasets
├── scripts/                        # Shared Python scripts
│   ├── generate_per_sample.py      # Query IGVF portal for metadata
│   ├── upload_to_gcp.py            # Upload files to GCS
│   ├── setup_dataset.py            # Bootstrap new dataset
│   └── validate_gcp_paths.py       # Verify GCS uploads
├── scratch/                        # Exploratory work (not committed)
└── ref/                            # Reference files (guide metadata, etc.)
```

## Pipeline Overview

The full workflow consists of five main stages:

1. **Upload data to IGVF Portal** - Data producers upload raw fastq and associated metadata
2. **Run CRISPR Pipeline** - Process raw data (locally or on GCP)
3. **Downstream QC** - Quality control analysis
4. **Energy Distance Analysis** - Perturbation effect quantification
5. **Gene Program Discovery** - cNMF and related analyses

### Stage 1: Upload Data to IGVF Portal
Fastq files should be uploaded to the [IGVF data portal](https://data.igvf.org/) by data producing groups.

Expected structure (see [example](https://data.igvf.org/analysis-sets/IGVFDS4389OUWU/)):
1. One or more measurement sets for scRNA-seq (split by cellular sub pools)
2. One auxiliary set per measurement set for gRNA-seq
3. Optionally, a second auxiliary set per measurement set for HTO-seq

TODO

### Stage 2: CRISPR Pipeline (GCP)

See [docs/CRISPR_PIPELINE.md](docs/CRISPR_PIPELINE.md) for detailed instructions.

1. **IGVF Analysis Set** - Accession ID for your dataset on the IGVF portal
2. **GCP Access** - Authenticated with `igvf-pertub-seq-pipeline` project
3. **CRISPR Pipeline** - Local clone of [IGVF CRISPR Pipeline](https://github.com/IGVF/CRISPR_Pipeline)

```bash
cd /path/to/tf_perturb_seq

# 2a. Generate sample metadata from IGVF portal
bash datasets/<DATASET_NAME>/1_generate_per_sample_metadata.sh

# 2b. Upload files to GCP (dry run first)
DRY_RUN=true bash datasets/<DATASET_NAME>/2_upload_to_gcp.sh
bash datasets/<DATASET_NAME>/2_upload_to_gcp.sh

# 2c. Patch gzipped files (if needed)
bash datasets/<DATASET_NAME>/patch_files.sh

# 2d. Run the pipeline
RUN_IN_BACKGROUND=true bash datasets/<DATASET_NAME>/3_run_CRISPR_pipeline.sh
```

### Stage 3: QC Pipeline
TODO

### Stage 4: Energy Distance Analysis
TODO

### Stage 5: Gene Program Discovery
TODO

## Datasets

### Technology Benchmark (WTC11)

Cross-lab comparison using the same WTC11 cell line:

| Dataset | Lab | Technology | Status |
|---------|-----|------------|--------|
| [Hon](datasets/Hon_WTC11-benchmark_TF-Perturb-seq/) | Hon | 10x 5' HT v2 w/ HTO | Complete |
| [Huangfu](datasets/Huangfu_WTC11-benchmark_TF-Perturb-seq/) | Huangfu | 10x 3' v3 | In progress |
| [Engreitz](datasets/Engreitz_WTC11-benchmark_TF-Perturb-seq/) | Engreitz | CC Perturb-seq | Waiting on seqspecs |
| [Gersbach GEM-Xv3](datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/) | Gersbach | 10x GEM-X 3' | Waiting on seqspecs |
| [Gersbach HTv2](datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/) | Gersbach | 10x HT v2 | Waiting on seqspecs |

See [technology-benchmark_WTC11_TF-Perturb-seq/](datasets/technology-benchmark_WTC11_TF-Perturb-seq/) for cross-dataset analyses.

### Production Datasets:

| Dataset | Description |
|---------|-------------|
| Hon_WTC11-cardiomyocyte-differentiation | Pools ABCD → cardiomyocytes (12 days) |
| Huangfu_HUES8-definitive-endoderm-differentiation | Pools ABCD → definitive endoderm |
| Huangfu_HUES8-embryonic-stemcell | Pools ABCD in HUES8 hESCs |
| Gersbach_WTC11-hepatocyte-differentiation | WTC11 → hepatocytes |

Full tracking: [Google Sheets](https://docs.google.com/spreadsheets/d/1ThGCNsxa2KUfboRweK1p_55GlgwaX1D19O5MTB6MbSs/edit?gid=0#gid=0)


