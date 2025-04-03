# IGVF TF Perturb-seq Pillar Project
Code for processing and analysis of CRISPR screening data targeting a large set of known transcription factors and epigenetic modifiers across multiple cell states

## TODO
- [ ] [Organize and facilitate participating labs to upload data to the portal](https://github.com/adamklie/tf_perturb_seq/issues/1)
- [ ] [Run CRISPR FG pipeline on pilot TF data](https://github.com/adamklie/tf_perturb_seq/issues/2)
- [ ] [Generate guide metadata files and determine whether one single file can be used for all of them](https://github.com/adamklie/tf_perturb_seq/issues/3)
- [ ] [Decide how we want to run inference with the CRISPR pipeline (local, global, pairwise, other?)](https://github.com/adamklie/tf_perturb_seq/issues/4)
- [ ] [Run CRISPR FG pipeline on all datasets](https://github.com/adamklie/tf_perturb_seq/issues/6)
- [ ] [Implement PySpade for global inference at target level](https://github.com/adamklie/tf_perturb_seq/issues/5)
- [ ] [Run energy distance analysis](https://github.com/adamklie/tf_perturb_seq/issues/7)
- [ ] [Run cNMF](https://github.com/adamklie/tf_perturb_seq/issues/8)

## Datasets
This repository contains multiple (`datasets/`) each with its own folder:

1. `Hon_TF_Perturb_Seq_Pilot` -- A pilot set of 55 TFs introduced transfected into WTC11 iPSCs and differentiated to cardiomyocytes (12 days)
2. [`Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq`](https://docs.google.com/presentation/d/1FDr7KE873Er1CezXQNqvqPKecoNEfOrBd8ixoTUsOSE/edit#slide=id.p) -- Pools ABCD transfected into WTC11 iPSCs and differentiated to cardiomyocytes (11 days)
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
Each dataset will be processed through the [IGVF CRISPR Pipeline](https://github.com/pinellolab/CRISPR_Pipeline)

Follow the instructions at the link above to run the pipeline. Briefly, the required inputs are:

1. Fastq files
First create an `analysis_sets.tsv` file that looks like this:
```
analysis_set_id	input_file_sets	measurement_sets	associated_auxiliary_sets
IGVFDS4389OUWU	IGVFDS3231RKUL	IGVFDS4852MAUM	/auxiliary-sets/IGVFDS6861GRJI/,/auxiliary-sets/IGVFDS5228AHGR/
IGVFDS4389OUWU	IGVFDS2378SIVI	IGVFDS7733OTCA	/auxiliary-sets/IGVFDS0997MTQA/,/auxiliary-sets/IGVFDS7027VRZT/
IGVFDS4389OUWU	IGVFDS7715GHKL	IGVFDS2628PYOH	/auxiliary-sets/IGVFDS8139DLEM/,/auxiliary-sets/IGVFDS0109SCQO/
IGVFDS4389OUWU	IGVFDS8280DBSV	IGVFDS7834BXBN	/auxiliary-sets/IGVFDS8280AMDA/,/auxiliary-sets/IGVFDS6231TYIS/
```
You can then use the `download_fastq.py` script to download all the fastqs to a `fastq_files` directory.

2. YAML Configuration Files:
 - `rna_seqspec.yml`: Defines RNA sequencing structure and parameters
 - `guide_seqspec.yml`: Specifies guide RNA detection parameters
 - `hash_seqspec.yml`: Defines cell hashing structure (required if using cell hashing)
 - `whitelist.txt`: List of valid cell barcodes

3. Metadata files
 - `guide_metadata.tsv` -- See https://github.com/adamklie/tf_perturb_seq/issues/3
 - `hash_metadata.tsv` -- include only if your protocol used HTOs

4. Pipeline configuration file: see https://github.com/adamklie/tf_perturb_seq/issues/4

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
