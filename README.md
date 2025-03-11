# IGVF TF Perturb-seq Pillar Project
Code for processing and analysis of CRISPR screening data targeting a large set of known transcription factors and epigenetic modifiers across multiple cell states

## Datasets
This repository is organized by datasets

1. `Hon_TF_Perturb_Seq_Pilot` -- A pilot set of TFs introduced transfected into WTC11 iPSCs and differentiated to cardiomyocytes (11 days)
2. `Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq` -- Pools ABCD transfected into WTC11 iPSCs and differentiated to cardiomyocytes (11 days)
3. `Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq` -- Pools ABCD transfected into HUES8 hESCs and differentiated to defenitive endoderm (?? days)
4. `Huangfu_HUES8-embryonic-stemcell_TF-Perturb-seq` -- Pools ABCD transfected into HUES8 hESCs
5. `Hon_H9-neuron-differentiation_TF-Perturb-seq` --
6. `Engreitz_WTC11-endothelial-differentiation_TF-Perturb-seq` --
7. `Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq` --
8. `Huangfu_HUES8-pancreatic-lineages-differentiation_TF-Perturb-seq`--

Dataset progress tracked here: https://docs.google.com/spreadsheets/d/1ThGCNsxa2KUfboRweK1p_55GlgwaX1D19O5MTB6MbSs/edit?gid=0#gid=0

## Processing
Each dataset was run through the [IGVF CRISPR Pipeline](https://github.com/pinellolab/CRISPR_Pipeline)

We follow the directory and configuration structure of that pipeline

## Downstream analysis
- [] aaa

## Guide library


## Preliminary TF Perturb-seq Synapse Upload
There is a tf_perturb_seq Synapse project: https://www.synapse.org/Synapse:syn63660123/files/

## Data (re)processing
Preprocessing with CellRanger:
 - Cell Ranger 8.0.1 w/ refdata-gex-GRCh38-2024-A (Huangfu)
Harmonized guide assignments
Harmonized differential expression calling

## Current outputs

```
tf_perturb_seq
├── Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq
│   ├── Cell_Ranger_Output
│   │   ├── filtered_feature_bc_matrix.h5
│   └── Perturbation_information
│       └── cell_x_sgrna_matrix.pkl
│       └── ref_feature.csv
│   ├── Differential_Expression_(PySpade0.1.5)
│   │   ├── filtered_cell_x_gene_matrix.h5ad
│   │   ├── global_DE.csv
│   │   └── local_DE.csv
```

* `Cell_Ranger_Output_(CR8.0.1)`
    * `filtered_feature_bc_matrix.h5` – output from CellRanger

* `Differential_Expression_(PySpade0.1.5)`
    * `filtered_cell_x_gene_matrix.h5ad` – input to differential expression

* `Perturbation_information`
    * `cell_x_sgrna_matrix.pkl`– guide_id x cell_id assignment matrix
    ${image?fileName=cell%5Fx%5Fsgrna%2Epng&align=None&scale=100&responsive=true&altText=}
    * `ref_feature.csv` - guide metadata
    ${image?fileName=ref%5Ffeature%2Epng&align=None&scale=100&responsive=true&altText=}
