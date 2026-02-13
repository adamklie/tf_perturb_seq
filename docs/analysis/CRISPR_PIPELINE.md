# Running the CRISPR Pipeline on GCP

This guide walks through running the IGVF CRISPR Pipeline on Google Cloud Platform from scratch.

## Overview

The workflow consists of four main steps:

1. **Generate sample metadata** - Query IGVF portal for dataset files and create a sample sheet
2. **Upload to GCP** - Transfer files from IGVF/S3 to Google Cloud Storage
3. **Patch files (optional)** - Decompress `.tsv.gz` files if needed
4. **Run pipeline** - Execute the Nextflow CRISPR pipeline on Google Cloud Batch

## Prerequisites

### 1. Software Requirements

| Software | Version | Installation |
|----------|---------|--------------|
| Python | 3.8+ | [python.org](https://python.org) |
| Google Cloud SDK | Latest | [cloud.google.com/sdk](https://cloud.google.com/sdk/docs/install) |
| Nextflow | 23.04+ | [nextflow.io](https://www.nextflow.io/) |
| Java | 17+ | [adoptium.net](https://adoptium.net/) |

### 2. Clone Required Repositories

```bash
# Clone this repository
git clone https://github.com/adamklie/tf_perturb_seq.git
cd tf_perturb_seq

# Clone the CRISPR Pipeline
git clone https://github.com/IGVF/CRISPR_Pipeline.git ../CRISPR_Pipeline
```

### 3. Install Python Dependencies

```bash
# Using uv (recommended)
uv sync

# Or with pip
pip install requests
```

### 4. Configure GCP Access

```bash
# Authenticate with GCP
gcloud auth login
gcloud auth application-default login

# Set project
gcloud config set project igvf-pertub-seq-pipeline

# Verify access to the bucket
gsutil ls gs://igvf-pertub-seq-pipeline-data/
```

### 5. Set IGVF API Credentials

```bash
# Add to your shell profile (~/.bashrc or ~/.zshrc)
export IGVF_API_KEY="your-key"
export IGVF_SECRET_KEY="your-secret"
```

Get credentials from [data.igvf.org](https://data.igvf.org) under your profile settings.

### 6. Configure Nextflow Tower (Optional)

```bash
# For pipeline monitoring
export TOWER_ACCESS_TOKEN="your-tower-token"
```

Get a token from [tower.nf](https://tower.nf).

## Dataset Structure

Each dataset follows this structure:

```
datasets/{DATASET_NAME}/
├── README.md                              # Dataset documentation
├── 1_generate_per_sample_metadata.sh      # Step 1 script
├── 2_upload_to_gcp.sh                     # Step 2 script
├── 3_run_CRISPR_pipeline.sh               # Step 4 script
├── patch_files.sh                         # Step 3 script (optional)
├── sample_metadata.csv                    # Output of step 1
├── sample_metadata_gcp_YYYY_MM_DD.csv     # Output of step 2
└── sample_metadata_gcp_*_patched.csv      # Output of step 3 (optional)
```

## Step-by-Step Workflow

### Step 1: Generate Sample Metadata

This step queries the IGVF data portal to create a CSV file with all sequencing files for your dataset.

#### Configure the Script

Edit `datasets/{DATASET_NAME}/1_generate_per_sample_metadata.sh`:

```bash
# Update these values:
BASE_DIR=/path/to/tf_perturb_seq
DATASET_NAME=Your_Dataset_Name
ACCESSION=IGVFDSXXXXXXXX  # Analysis set ID from IGVF portal
```

Find your Analysis Set accession at [data.igvf.org](https://data.igvf.org).

#### Run

```bash
cd /path/to/tf_perturb_seq
bash datasets/{DATASET_NAME}/1_generate_per_sample_metadata.sh
```

#### Output

Creates `sample_metadata.csv` with columns:
- `R1_path`, `R2_path` - IGVF file accessions
- `file_modality` - `scRNA`, `gRNA`, or `hash`
- `measurement_sets` - Measurement set accession
- `sequencing_run`, `lane` - Sequencing metadata
- `seqspec` - Path to seqspec YAML file
- `barcode_onlist` - Cell barcode whitelist file
- `guide_design` - Guide sequences file
- `barcode_hashtag_map` - HTO mapping (if applicable)

### Step 2: Upload to GCP

This step transfers files from IGVF S3 storage to Google Cloud Storage.

#### Configure the Script

Edit `datasets/{DATASET_NAME}/2_upload_to_gcp.sh`:

```bash
# Update these values:
BASE_DIR=/path/to/tf_perturb_seq
DATASET_NAME=Your_Dataset_Name
PROJECT=igvf-pertub-seq-pipeline
GCS_BUCKET=igvf-pertub-seq-pipeline-data
```

#### Run (Dry Run First)

```bash
# Preview what will be uploaded
DRY_RUN=true bash datasets/{DATASET_NAME}/2_upload_to_gcp.sh

# Perform actual upload
bash datasets/{DATASET_NAME}/2_upload_to_gcp.sh
```

#### Verify Upload

```bash
gsutil ls gs://igvf-pertub-seq-pipeline-data/{DATASET_NAME}/
```

#### Output

Creates `sample_metadata_gcp_YYYY_MM_DD.csv` with GCS paths replacing IGVF accessions.

### Step 3: Patch Files (Optional)

Some pipeline components require uncompressed TSV files. If your `barcode_onlist` or `guide_design` files are gzipped (`.tsv.gz`), decompress them.

#### Configure the Script

Edit `datasets/{DATASET_NAME}/patch_files.sh`:

```bash
BUCKET="gs://igvf-pertub-seq-pipeline-data"
DATASET="Your_Dataset_Name"
DATE="YYYY_MM_DD"  # Match your upload date

# Get these paths from sample_metadata_gcp_*.csv
BARCODE_ONLIST_SRC="${BUCKET}/${DATASET}/${DATE}/.../IGVFFIXXXXXX.tsv.gz"
GUIDE_DESIGN_SRC="${BUCKET}/${DATASET}/${DATE}/.../IGVFFIXXXXXX.tsv.gz"
```

#### Run

```bash
bash datasets/{DATASET_NAME}/patch_files.sh
```

#### Create Patched Sample Metadata

After running the patch script, create a patched CSV:

```bash
# Using Python to update paths
python3 << 'EOF'
import csv

input_file = 'datasets/{DATASET_NAME}/sample_metadata_gcp_YYYY_MM_DD.csv'
output_file = 'datasets/{DATASET_NAME}/sample_metadata_gcp_YYYY_MM_DD_patched.csv'

PATCH_DIR = "gs://igvf-pertub-seq-pipeline-data/{DATASET_NAME}/YYYY_MM_DD/patch"
barcode_patched = f"{PATCH_DIR}/IGVFFIXXXXXX.tsv"
guide_patched = f"{PATCH_DIR}/IGVFFIXXXXXX.tsv"

with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
    reader = csv.DictReader(infile)
    writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames)
    writer.writeheader()
    for row in reader:
        row['barcode_onlist'] = barcode_patched
        row['guide_design'] = guide_patched
        writer.writerow(row)
EOF
```

### Step 4: Run CRISPR Pipeline

Execute the Nextflow pipeline on Google Cloud Batch.

#### Configure the Script

Edit `datasets/{DATASET_NAME}/3_run_CRISPR_pipeline.sh`:

```bash
DATASET_NAME=Your_Dataset_Name
BASE_DIR=/path/to/tf_perturb_seq/datasets/${DATASET_NAME}

# Use patched CSV if you ran step 3, otherwise use the regular one
SAMPLE_METADATA=$BASE_DIR/sample_metadata_gcp_YYYY_MM_DD_patched.csv

PIPELINE_PATH=/path/to/CRISPR_Pipeline
OUTDIR=gs://igvf-pertub-seq-pipeline-data/${DATASET_NAME}/YYYY_MM_DD/outs
```

#### Run

```bash
# Run in foreground (blocks terminal)
bash datasets/{DATASET_NAME}/3_run_CRISPR_pipeline.sh

# Run in background (recommended for long runs)
RUN_IN_BACKGROUND=true bash datasets/{DATASET_NAME}/3_run_CRISPR_pipeline.sh
```

#### Monitor Progress

```bash
# Check local log file
tail -f datasets/{DATASET_NAME}/logs/{DATASET_NAME}_crispr_pipeline_*.log

# Monitor on Nextflow Tower (if configured)
# https://tower.nf
```

#### Output

Pipeline outputs are written to `$OUTDIR` on GCS:
- Gene expression matrices
- Guide assignment results
- QC reports
- MuData objects

## Quick Reference

```bash
# Full workflow for a new dataset
cd /path/to/tf_perturb_seq

# Step 1: Generate metadata
bash datasets/{DATASET_NAME}/1_generate_per_sample_metadata.sh

# Step 2: Upload to GCP
DRY_RUN=true bash datasets/{DATASET_NAME}/2_upload_to_gcp.sh
bash datasets/{DATASET_NAME}/2_upload_to_gcp.sh

# Step 3: Patch files (if needed)
bash datasets/{DATASET_NAME}/patch_files.sh
# Then create patched CSV manually

# Step 4: Run pipeline
RUN_IN_BACKGROUND=true bash datasets/{DATASET_NAME}/3_run_CRISPR_pipeline.sh
```

## Troubleshooting

### IGVF API 403 Forbidden

Ensure credentials are set:
```bash
export IGVF_API_KEY="your-key"
export IGVF_SECRET_KEY="your-secret"
```

### GCS Permission Denied

Re-authenticate:
```bash
gcloud auth login
gcloud auth application-default login
gcloud config set project igvf-pertub-seq-pipeline
```

### "Cannot invoke method toLowerCase() on null object"

The sample metadata CSV has incorrect format. Check:
- File is comma-separated (not TSV)
- `file_modality` contains `scRNA`, `gRNA`, or `hash`

### GCS Paths Don't Exist

Validate paths before running:
```bash
python scripts/validate_gcp_paths.py --input sample_metadata_gcp_*.csv
```

### Nextflow Java Version Error

Install Java 17+:
```bash
# macOS
brew install openjdk@17

# Verify
java -version
```

### Pipeline Fails on Compressed Files

Run the patch script to decompress `.tsv.gz` files:
```bash
bash datasets/{DATASET_NAME}/patch_files.sh
```

## Creating a New Dataset

1. Copy the template:
   ```bash
   cp -r templates/dataset_template datasets/{NEW_DATASET_NAME}
   ```

2. Update all scripts with your dataset name and IGVF accession

3. Add seqspec YAML files for each modality

4. Follow the workflow above

## Related Resources

- [IGVF Data Portal](https://data.igvf.org)
- [CRISPR Pipeline Repository](https://github.com/IGVF/CRISPR_Pipeline)
- [Nextflow Documentation](https://www.nextflow.io/docs/latest/)
- [Nextflow Tower](https://tower.nf)
