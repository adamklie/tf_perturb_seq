# Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq

## Overview

| Property | Value |
|----------|-------|
| **Lab** | Huangfu |
| **Cell Line** | HUES8 |
| **Differentiation** | Definitive endoderm |
| **IGVF Analysis Set** | TODO: Find accession on portal |
| **Status** | Setting up pipeline scripts |

## Data Description

Production TF Perturb-seq dataset from the Huangfu lab. HUES8 hESCs differentiated to definitive endoderm, with CRISPRi perturbation of the full TF library (~2000 targets). Pools ABCD.

## Pipeline Status

- [ ] Step 0: **BLOCKED** — Need analysis set created on IGVF portal (see below)
- [ ] Step 1: Generate per-sample metadata (need portal accession ID)
- [ ] Step 2: Upload to GCP
- [ ] Step 3: Patch compressed files
- [ ] Step 4: Run CRISPR pipeline
- [ ] Step 5: QC analysis
- [ ] Step 6: Energy distance
- [ ] Step 7: cNMF

## Legacy Data (to be archived)

The following were from an earlier internal processing run and are **not** from the standardized CRISPR FG pipeline. They should be archived to `scratch/` once the new pipeline run is complete:

- `mdata_filtered.h5mu` (3.8 GB, Mar 2025)
- `Cell_Ranger_Output/`
- `Perturbation_information/`
- `Transcriptome_Analysis/`

## Portal Status (audited 2026-03-25)

Measurement sets (8) and auxiliary sets (8 gRNA) exist on the portal. Construct library set is IGVFDS3299AXST with guide file IGVFFI8270UPKB (released).

**Blockers:**
1. **No analysis set exists** — someone needs to create one grouping 8 MS + 8 aux sets
2. **No seqspecs** on R1 files — need upload or fallback YAMLs

Measurement sets: IGVFDS0788TUIL, IGVFDS1313VGHL, IGVFDS1403KZYV, IGVFDS1437XFJK, IGVFDS2129VHBD, IGVFDS2520ZEYH, IGVFDS3434HGEI, IGVFDS4079VXPX

See `docs/DATA.md` for full MS-to-aux mapping.

**ACTION:** Coordinate with Denis Torre / DACC to create analysis set and upload seqspecs.

## Notes

- Data was previously shared via Dropbox from the Huangfu lab
- Raw fastqs are in `fastq/`
- These datasets use 10x 3' v3 chemistry (same as benchmark)
- Config can likely be adapted from `datasets/Huangfu_WTC11-benchmark_TF-Perturb-seq/Huangfu_WTC11-benchmark_TF-Perturb-seq_2026_03_11.config` but will need changes for the larger guide library and potentially different MOI
