## Data (re)processing
Preprocessing with CellRanger:
 - Cell Ranger 8.0.1 w/ refdata-gex-GRCh38-2024-A (Huangfu)
Harmonized guide assignments
Harmonized differential expression calling

## Preliminary TF Perturb-seq Synapse Upload
There is a tf_perturb_seq Synapse project: https://www.synapse.org/Synapse:syn63660123/files/

### Create a machine readable, descriptive name for your dataset. Examples:
- `Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq`
- `Hon_H9-neuron-differentiation_TF-Perturb-seq`
- `Engreitz_WTC11-endothelial-differentiation_TF-Perturb-seq`
- `Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq`
- `Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq`
- `Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq`

### Create the following directory structure locally

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

### With the following file descriptions
* `Cell_Ranger_Output_(CR8.0.1)`
    * `filtered_feature_bc_matrix.h5` – output from CellRanger

* `Differential_Expression_(PySpade0.1.5)`
    * `filtered_cell_x_gene_matrix.h5ad` – input to differential expression

* `Perturbation_information`
    * `cell_x_sgrna_matrix.pkl`– guide_id x cell_id assignment matrix
    ${image?fileName=cell%5Fx%5Fsgrna%2Epng&align=None&scale=100&responsive=true&altText=}
    * `ref_feature.csv` - guide metadata
    ${image?fileName=ref%5Ffeature%2Epng&align=None&scale=100&responsive=true&altText=}

## Upload to Synapse: 
Use the following template notebook: https://gist.github.com/adamklie/beddb4de448f20740843cd22c35bd181# tf_perturb_seq
