# IGVF TF Perturb-seq Pillar Project
Code for processing and analysis of CRISPR screening data targeting a large set of known transcription factors and epigenetic modifiers across multiple cell states

## TODO
- [ ] [Organize and facilitate participating labs to upload data to the portal](https://github.com/adamklie/tf_perturb_seq/issues/1)
- [ ] [Run CRISPR FG pipeline on pilot TF data](https://github.com/adamklie/tf_perturb_seq/issues/2)
- [ ] [Generate guide metadata files and determine whether one single file can be used for all of them](https://github.com/adamklie/tf_perturb_seq/issues/3)
- [ ] [Decide how we want to run inference with the CRISPR pipeline (local, global, pairwise, other?)](https://github.com/adamklie/tf_perturb_seq/issues/4)
- [ ] [Run CRISPR FG pipeline on all datasets]()
- [ ] [Implement PySpade for global inference at target level](https://github.com/adamklie/tf_perturb_seq/issues/5)
- [ ] [Run energy distance analysis](https://github.com/adamklie/tf_perturb_seq/issues/7)
- [ ] [Run cNMF](https://github.com/adamklie/tf_perturb_seq/issues/8)

## Datasets
This repository contains multiple (`datasets/`). Each dataset is a separate folder containing all the necessary files to run the initial CRISPR pipeline. The datasets are:

1. `Hon_TF_Perturb_Seq_Pilot` -- A pilot set of 55 TFs introduced transfected into WTC11 iPSCs and differentiated to cardiomyocytes (12 days)
2. `Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq` -- Pools ABCD transfected into WTC11 iPSCs and differentiated to cardiomyocytes (11 days)
3. `Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq` -- Pools ABCD transfected into HUES8 hESCs and differentiated to defenitive endoderm (?? days)
4. `Huangfu_HUES8-embryonic-stemcell_TF-Perturb-seq` -- Pools ABCD transfected into HUES8 hESCs
5. `Hon_H9-neuron-differentiation_TF-Perturb-seq` --
6. `Engreitz_WTC11-endothelial-differentiation_TF-Perturb-seq` --
7. `Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq` --
8. `Huangfu_HUES8-pancreatic-lineages-differentiation_TF-Perturb-seq`--

Dataset details and status is tracked here: https://docs.google.com/spreadsheets/d/1ThGCNsxa2KUfboRweK1p_55GlgwaX1D19O5MTB6MbSs/edit?gid=0#gid=0

## Data upload
Fastq files should be uploaded to the IGVF data portal by relevant data producing groups. This needs to be done systematically to ensure that data can be properly processed. TODO

## Processing
Each dataset will be run through the [IGVF CRISPR Pipeline](https://github.com/pinellolab/CRISPR_Pipeline)

Follow the instructions at the link above to run the pipeline. Briefly, the required inputs are:

1. . Fastq files need to be downloaded to the machine where the pipeline will be run. If downloading you can use the `download_fastq.py` file. You will need to create a TSV file matching `per-sample_file.tsv`. Once downloaded, they should be in a folder titled `fastq_files`.

2. `guide_metadata.tsv` -- THIS SHOULD BE THE SAME FOR EVERY DATASET.

3. `hash_metadata.tsv` -- include only if your protocol used HTOs

4. An updated configuration file `pipeline_input.config`

## Downstream analysis

### PySpade differential expression inference

### Energy distance analysis

### cNMF

## Preliminary TF Perturb-seq Synapse Upload
A preliminary version of 3 datasets is available in the following Synapse folder: https://www.synapse.org/Synapse:syn64423137

If you do not have access to Synapse, please request access from the a member of the IGVF Stanford DACC.

## Processing steps
ALl 3 datasets were first processed with CellRanger to get gene UMI count matrices
 - Cell Ranger 8.0.1
 - refdata-gex-GRCh38-2024-A
Each lab did there own guide assignment
 - Hon Lab: FBA
 - Huangfu Lab: SCEPTRE

## Current outputs
* `Cell_Ranger_Output_(CR8.0.1)`
    * `filtered_feature_bc_matrix.h5` – output from CellRanger
    * `barcodes.tsv.gz` – barcodes from CellRanger output

* `Perturbation_information`
    * Two types of output matrix are generated here depending on the tool (should be sparse matrices in the future)
        * `cell_x_sgrna_matrix.pkl`– guide_id x cell_id assignment matrix from PySpade. This is a pickled pandas dataframe that takes a long time to load in
        * `sgrna_design_matrix_filtered_combined_control_final.csv` -- guide_id x gene_id assignment matrix from SCEPTRE. This is a csv file that also takes a long time to load in


## More details on the pipeline

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
