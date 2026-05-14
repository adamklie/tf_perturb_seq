# Data

## Datasets

| Dataset | Folder |
|---|---|
| Hon WTC11 cardiomyocyte | [`Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/`](Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/) |
| Huangfu HUES8 definitive endoderm | [`Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/`](Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/) |
| Huangfu HUES8 embryonic stem cell | [`Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/`](Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/) |
| Gersbach WTC11 hepatocyte | [`Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/`](Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/) |
| Engreitz WTC11 endothelial | [`Engreitz_WTC11-endothelial-cells_TF-Perturb-seq/`](Engreitz_WTC11-endothelial-cells_TF-Perturb-seq/) |

Each per-dataset folder contains only a `README.md` linking to the dataset's Synapse folder. Bulky pipeline outputs (CRISPR pipeline, cNMF, energy distance) live on Synapse under `syn64423137/2026_UTSW/`, not in this directory.

## Cross-cutting files

| Path | What it is |
|---|---|
| [`schemas/`](schemas/) | JSON schemas for every output table and metadata file. See [`schemas/README.md`](schemas/README.md). |
| [`scripts/`](scripts/) | Python scripts that produce or mirror the contents of this directory. See [`scripts/README.md`](scripts/README.md). |
