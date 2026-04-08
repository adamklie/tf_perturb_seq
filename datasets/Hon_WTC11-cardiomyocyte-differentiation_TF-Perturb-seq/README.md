# Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq

## Overview

| Property | Value |
|----------|-------|
| **Lab** | Hon (UTSW) |
| **Cell Line** | WTC11 (GM25256) |
| **Differentiation** | Cardiomyocyte (cardiac muscle) |
| **IGVF Analysis Set** | [IGVFDS6332VCTO](https://data.igvf.org/analysis-sets/IGVFDS6332VCTO) (released, 78 inputs) |
| **Guide Library** | [IGVFDS3299AXST](https://data.igvf.org/construct-library-sets/IGVFDS3299AXST) (guide file: [IGVFFI8270UPKB](https://data.igvf.org/tabular-files/IGVFFI8270UPKB)) |
| **Technology** | 10x (HTO multiplexed) |
| **Scale** | 32 measurement sets, 26 scRNA lanes, 52 gRNA lanes, 26 HTO lanes |
| **Status** | Prepping for CRISPR pipeline run (Lucas) |

## Pipeline Status

- [x] Step 0: Portal has analysis set, measurement sets, auxiliary sets
- [x] Step 1: Sample metadata generated (exists on GCP)
- [x] Step 2: Fastqs already on GCP
- [ ] Step 3: Patch compressed files (seqspec, barcode_onlist, guide_design are .gz)
- [ ] Step 4: Run CRISPR pipeline (Lucas to run on GCP)
- [ ] Step 5: QC analysis
- [ ] Step 6: Energy distance
- [ ] Step 7: cNMF

## GCP Data Location

```
gs://igvf-pertub-seq-pipeline-data/WTC11_CM_TF_PerturbSeq/
├── WTC11_CM_TF_PerturbSeq_per_sample_accessions.tsv  # Portal accession IDs
├── WTC11_CM_TF_PerturbSeq_per_sample_paths.tsv        # GCS paths (104 rows + header)
├── 2024/  # Fastqs and reference files
└── 2025/  # Seqspecs
```

**Sample metadata TSV** (104 rows): 26 scRNA + 52 gRNA + 26 HTO lanes
- Already has GCS paths for R1/R2, seqspec, barcode_onlist, guide_design
- Columns: R1_path, R1_md5sum, R2_path, R2_md5sum, file_modality, file_set, measurement_sets, sequencing_run, lane, flowcell_id, index, seqspec, barcode_onlist, onlist_method, strand_specificity, guide_design, barcode_hashtag_map

## What needs to happen

1. **Regenerate sample metadata from portal** — The existing TSV on GCP (`WTC11_CM_TF_PerturbSeq_per_sample_paths.tsv`) is old (2024-2025) and may be stale. Run `1_generate_per_sample_metadata.sh` to generate fresh metadata from analysis set IGVFDS6332VCTO.
2. **Check existing GCP fastqs** — `gs://igvf-pertub-seq-pipeline-data/WTC11_CM_TF_PerturbSeq/` has fastqs from a prior upload. Need to verify these are still current or if they should be cleaned up (check with Ian/Hon lab). Re-upload from portal if needed via `2_upload_to_gcp.sh`.
3. **Patch .gz files** — seqspec (.yaml.gz), barcode_onlist (.tsv.gz), and guide_design (.csv.gz) need decompression for the pipeline.
4. **Create pipeline config** — Adapt from Hon benchmark config. Key differences from benchmark: full TF library (~2000 targets vs ~50), HTO multiplexing, more lanes.
5. **Hand off to Lucas** for GCP pipeline execution.

## Legacy Data

`HonLabInternal/` contains outputs from the Hon lab's internal processing pipeline (Cell_Ranger_Output/, Perturbation_information/, Transcriptome_Analysis/). These are not from the standardized CRISPR FG pipeline and should be archived once the new run completes.

## Portal Details (audited 2026-03-25)

| Component | Count | Status |
|-----------|:---:|--------|
| Measurement sets (scRNA) | 26 | On portal, in TFP3 collection |
| Auxiliary sets (gRNA) | 26 × 2 = 52 | On portal |
| Auxiliary sets (HTO) | 26 | On portal |
| Construct library set | 1 (IGVFDS3299AXST) | Released, guide file present |
| Analysis set | 1 (IGVFDS6332VCTO) | Released, 78 inputs |
| Seqspecs | Present (.yaml.gz on GCP) | Need decompression |
