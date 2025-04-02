# Files provided by labs
`Hon_sgRNA_index_dacc_annot_reference.csv` -- 14,358 total guides, latest file from Hon lab

`Huangfu_ref_feature.csv` -- 14,364 total guides, latest file from Huangfu lab

# Files extracted from [here](https://docs.google.com/spreadsheets/d/1WcVgLllWrtO_-h5Ry07WS_ovomRcyznQdZrx2KOntEU/edit?gid=1430289032#gid=1430289032)
`negative_controls.tsv` (Neg control sgRNAs (targeting) sheet) -- 600 targeting guides, 100 targets, 6 guides per target

`non_targeting.tsv` (Neg control sgRNAs (non-targeting) sheet) -- 600 non-targeting guides

`positive_controls.tsv` (Pos control sgRNAs sheet) -- 19 targeting guides, 19 targets, 1 guide per target

`target_genes.tsv` (Target genes + sgRNAs (with alt prom) sheet) -- 13,260 targeting guides, 1983 targets (some with multiple promoters), 6 guides per target

# Master generated from merging files
`sgRNA_id_master.tsv` -- 14,364 total guides, merged guides from Hon and Huangfu labs (see `guide_mapping.ipynb`)
