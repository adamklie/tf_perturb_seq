# IGVF TF Perturb-seq Pillar Project

Code for processing and analysis of CRISPR screening data targeting a large set of known transcription factors and epigenetic modifiers across multiple cell states.

## Quick Start

```bash
# Clone the repository
git clone https://github.com/adamklie/tf_perturb_seq.git
cd tf_perturb_seq

# Install dependencies (using uv)
uv sync

# Set up a new dataset
python scripts/setup_dataset.py --name "Lab_CellLine-state_TF-Perturb-seq" --accession IGVFDSXXXXXX

# Follow the pipeline steps in docs/CRISPR_PIPELINE.md
```

## Repository Structure

```
tf_perturb_seq/
├── datasets/                    # Dataset directories (one per experiment)
│   └── <dataset_name>/
│       ├── bin/                 # Analysis scripts
│       │   ├── 1_CRISPR_pipeline/
│       │   ├── 2_qc/
│       │   └── 3_gene_program_discovery/
│       ├── dataset_config.yaml  # Dataset-specific parameters
│       └── sample_metadata.csv  # Generated sample sheet
├── docs/                        # Documentation
│   └── CRISPR_PIPELINE.md      # Pipeline workflow guide
├── scripts/                     # Shared Python scripts
│   ├── generate_per_sample.py  # Query IGVF portal for metadata
│   ├── upload_to_gcp.py        # Upload files to GCS
│   ├── setup_dataset.py        # Bootstrap new dataset
│   └── validate_gcp_paths.py   # Verify GCS uploads
├── src/                         # Python modules
├── ref/                         # Reference files (guide metadata, etc.)
└── templates/                   # Dataset templates
```

## Pipeline Overview

The processing pipeline consists of three main steps:

1. **Generate Sample Metadata** - Query IGVF portal for dataset files and create sample sheet
2. **Upload to GCP** - Transfer files from IGVF/S3 to Google Cloud Storage
3. **Run CRISPR Pipeline** - Execute Nextflow pipeline on Google Cloud Batch

See [docs/CRISPR_PIPELINE.md](docs/CRISPR_PIPELINE.md) for detailed instructions.

## TODO

### Run processing pipeline on all datasets
- [ ] [Run CRISPR FG pipeline on all datasets](https://github.com/adamklie/tf_perturb_seq/issues/6)
- [ ] [Harmonized pipeline config](https://github.com/adamklie/tf_perturb_seq/issues/4)

### Downstream analysis
- [ ] [Run energy distance analysis](https://github.com/adamklie/tf_perturb_seq/issues/7)
- [ ] [Run cNMF and gene program evaluations on each dataset](https://github.com/adamklie/tf_perturb_seq/issues/8)

## Datasets

This repository contains multiple datasets in `datasets/`, each with its own folder:

1. `Hon_TF_Perturb_Seq_Pilot` -- A pilot set of 55 TFs introduced transfected into WTC11 iPSCs and differentiated to cardiomyocytes (12 days)
2. [`Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq`](https://docs.google.com/presentation/d/1FDr7KE873Er1CezXQNqvqPKecoNEfOrBd8ixoTUsOSE/edit#slide=id.p) -- Pools ABCD transfected into WTC11 iPSCs and differentiated to cardiomyocytes (12 days)
3. [`Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq`](https://docs.google.com/presentation/d/1anFVk5SsNikMiVymgnp2283khwRSh6V5dMlTs4m_V0M/edit#slide=id.p) -- Pools ABCD transfected into HUES8 hESCs and differentiated to defenitive endoderm (?? days)
4. [`Huangfu_HUES8-embryonic-stemcell_TF-Perturb-seq`](https://docs.google.com/presentation/d/13mw1IA4diHEHsFRSyK25IoGjpMXCQ12bvPJ2-tt73Ek/edit#slide=id.p) -- Pools ABCD transfected into HUES8 hESCs
5. `Hon_H9-neuron-differentiation_TF-Perturb-seq` --
6. `Engreitz_WTC11-endothelial-differentiation_TF-Perturb-seq` --
7. `Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq` --
8. `Huangfu_HUES8-pancreatic-lineages-differentiation_TF-Perturb-seq`--

Dataset details and status is tracked here: https://docs.google.com/spreadsheets/d/1ThGCNsxa2KUfboRweK1p_55GlgwaX1D19O5MTB6MbSs/edit?gid=0#gid=0

## Data upload
Fastq files should be uploaded to the [IGVF data portal](https://data.igvf.org/) by relevant data producing groups. This analysis set for `Hon_TF_Perturb_Seq_Pilot` is a good example of the expected structure: https://data.igvf.org/analysis-sets/IGVFDS4389OUWU/.

1. One or more measurement sets associated with the single-cell RNA-sequencing (likely split by cellular sub pools)
2. One associated auxiliary set for each of measurement set for gRNA-sequencing
3. Optionally a second auxiliary set for each measurement set for HTO-sequencing

Each measurement and auxiliary set needs a set of fastqs in it, as well as seqspecs for each type of sequencing.

## Processing

Each dataset is processed through the [IGVF CRISPR Pipeline](https://github.com/pinellolab/CRISPR_Pipeline) running on Google Cloud.

**For detailed instructions, see [docs/CRISPR_PIPELINE.md](docs/CRISPR_PIPELINE.md).**

### Quick Reference

```bash
# 1. Generate sample metadata from IGVF portal
./bin/1_CRISPR_pipeline/1_generate_per_sample_metadata.sh

# 2. Upload files to GCP
./bin/1_CRISPR_pipeline/2_upload_to_gcp.sh

# 3. Run the pipeline
./bin/1_CRISPR_pipeline/3_run_CRISPR_pipeline.sh
```

### Required Inputs

1. **IGVF Analysis Set** - Accession ID for your dataset on the IGVF portal
2. **Seqspec YAML files** - Read structure definitions for RNA, guide, and hash (if applicable)
3. **Guide metadata** - See [Issue #3](https://github.com/adamklie/tf_perturb_seq/issues/3)
4. **Pipeline config** - Dataset-specific parameters in `dataset_config.yaml`

## Downstream analysis

### Basic transcriptomic analysis

### E-distance analysis

### cNMF gene program analysis

## Preliminary TF Perturb-seq processing and analysis
A preliminary version of 3 processed datasets is available in the following Synapse folder: https://www.synapse.org/Synapse:syn64423137

If you do not have access to Synapse, please request access from the a member of the IGVF Stanford DACC.

## Processing steps
All 3 datasets were first processed with CellRanger to get gene UMI count matrices
 - Cell Ranger 8.0.1
 - refdata-gex-GRCh38-2024-A

Each lab did there own guide assignment
 - Hon Lab: FBA
 - Huangfu Lab: SCEPTRE

 Guide metadata file - since each lab did its own thing with the guide assignment, they ended up using different guide annotation metadata files
 - I generated a merged cleaned up version for more uniform downstream analysis: [file](https://github.com/adamklie/tf_perturb_seq/blob/main/ref/sgRNA_id_master.tsv) |  [methodology](https://docs.google.com/presentation/d/1EIflIKvOoeS7UC_fGktfU-Wvg-k22DvO5ChGfUJtY_Y/edit#slide=id.g341d90379bb_0_32) | [code](https://github.com/adamklie/tf_perturb_seq/blob/main/scripts/guide_mapping.ipynb)

## Current outputs
* `Cell_Ranger_Output`
    * `barcodes.tsv.gz` – barcodes from CellRanger output
    * `filtered_feature_bc_matrix.h5` – cell x gene matrix from CellRanger in h5 format
    * `web_summary.html` -- HTML summary output from CellRanger `aggr` command

* `Perturbation_Information` -- Two formats for cell x guide matrix are generated here depending on the guide assignment tool run
    * `cell_x_sgrna_matrix.pkl`– guide_id x cell_id assignment matrix from PySpade (Hon lab datasets)
    * `sgrna_design_matrix_filtered_combined_control_final.csv` -- guide_id x gene_id assignment matrix from SCEPTRE (Huangfu lab datasets)
    * I generated a more standardized version of the guide matrix for each: [example output](https://www.synapse.org/Synapse:syn64003380) | [code](https://github.com/adamklie/tf_perturb_seq/blob/main/datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/bin/1_create_guide_anndata.ipynb)
        * `guide_assignment_matrix/` - 10x style matrix markdown directory
        * `guide_assignment_matrix.csv` - comma separated file for input into tools like PySpade
        * `guide_assignment_matrix.h5ad` - AnnData format
          
* `Transcriptome_Analysis` -- output from quick and dirty RNA analysis: [pipeline script](https://github.com/adamklie/tf_perturb_seq/blob/main/datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/bin/2_gene_UMI_pipeline.sh) | [example output](https://www.synapse.org/Synapse:syn65491895) | [methodolgy](https://docs.google.com/presentation/d/1MBH1zJmJLu4Dlzb4oRwpqojQXPPtw8hu9Q-lhyrQEfc/edit#slide=id.g341dd80d490_0_5)
      
 * `mdata_filtered.h5mu` -- a preliminary MuData object that is a rough estimation of the structure expected from the CRISPR FG pipeline

## More details on the data

### Guide library
Summary of sgRNA designs: https://docs.google.com/document/d/1Tup8lX42O_sVPxU-kzAFtcM0GAZMJqer7vU3zbM0uFA/edit?tab=t.0#bookmark=id.qk5wiev070

Pool A-C (target genes + positive controls)
Targeting 1983 TFs and epigenetic genes (2210 promoters) x 6 sgRNAs = 13,260 sgRNAs
Pool A: 4420 sgRNAs
Pool B: 4420 sgRNAs
Pool C: 4420 sgRNAs
Synthesized as 3 separate pools of sgRNAs (1+2, 3+4, and 5+6).

Pool D (Positive and Negative controls):
Targeting negative controls: 100 olfactory genes x 6 sgRNAs = 600 sgRNAs (5% of screen) [sgRNA design from Table S4 of Replogle et al Elife 2023]
Non-targeting negative controls: 600 sgRNAs (5% of screen)
Each pool will also have 19 positive control sgRNAs (plus any system-specific ones), each represented 6-10 times to ensure detection: 127 sgRNAs
Total: 1200 + 127 sgRNAs = 1327 sgRNAs

Spreadsheet https://docs.google.com/spreadsheets/d/1WcVgLllWrtO_-h5Ry07WS_ovomRcyznQdZrx2KOntEU/edit?gid=1430289032#gid=1430289032
