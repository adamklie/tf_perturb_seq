# Production TF Perturb-seq — cross-dataset manifests and drivers

This is a **cross-dataset** folder (sibling to `datasets/technology-benchmark_WTC11_TF-Perturb-seq/`), not a per-dataset folder. It holds the manifests + scripts used to run analyses across multiple production datasets in a uniform way.

## Layout

```
production_TF-Perturb-seq/
├── README.md
├── bin/
│   └── qc_array.sh                            # SLURM array driver for QC across datasets
├── manifests/
│   ├── production_h5mu.tsv                    # canonical list of per-dataset inference_mudata.h5mu paths
│   ├── qc.tsv                                 # QC parameter manifest (per-dataset)
│   └── 2026_02_05_pipeline_comparison.tsv     # historical pipeline-version comparison
└── results/                                   # cross-dataset analysis outputs (empty for now)
```

## Production datasets in scope

See `manifests/production_h5mu.tsv` for the full list. As of 2026-05-11:
- `Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq`
- `Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq`
- `Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq`
- `Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq`
- `Engreitz_WTC11-endothelial-cells_TF-Perturb-seq` (blocked; no data)
