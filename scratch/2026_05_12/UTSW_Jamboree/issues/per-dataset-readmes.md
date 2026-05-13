# Issue 6 — Per-dataset READMEs

**Status**: ✅ mostly complete (2026-05-09). Top-level + per-analysis READMEs landed for all 5 production datasets + the HTv2 testbed. Will be refreshed in place as each output's status changes.

**Owner of fix**: us.

## TL;DR

The HTv2 testbed already has [`datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/crispr_pipeline/README.md`](../datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/crispr_pipeline/README.md). We need the equivalent for each production dataset, mostly auto-fillable from `synapse_paths.tsv` + `experimental_metadata.tsv`.

## Required files

Each production dataset gets:

```
datasets/<dataset>/
├── README.md                       # dataset-level pointer + status
├── crispr_pipeline/README.md       # what's on Synapse for CRISPR pipeline
├── cnmf/README.md                  # what's on Synapse for cNMF (TBD until runs land)
└── energy_distance/README.md       # what's on Synapse for energy distance
```

Datasets:

- [x] `Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq` — top-level + crispr_pipeline + cnmf + energy_distance READMEs all present
- [x] `Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq` — full set
- [x] `Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq` — full set
- [x] `Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq` — full set (per-analysis files note "awaiting Sara")
- [x] `Engreitz_WTC11-endothelial-cells_TF-Perturb-seq` — top-level + cnmf + energy_distance READMEs (no crispr_pipeline placeholder; the top-level + Issue #4 cover it)
- [x] `Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2` (testbed) — top-level + crispr_pipeline + cnmf READMEs

## Template

The HTv2 CRISPR pipeline README is the template. Each per-dataset / per-analysis README should include:

- 1-line description of what this dataset's analysis is
- The canonical run name + run label (from `experimental_metadata.tsv` `canonical_run_label`)
- Quick-fact table: source (HPC / GCS / local), Synapse path, schema link, bundle size
- "What's in this folder" — the canonical layout for that analysis
- "What's not here" — intentional exclusions, missing pieces, source-run quirks
- Pointers to schemas + analysis walkthroughs

## Acceptance criteria

- [ ] Each production dataset has a top-level `README.md` summarizing what's on Synapse for it.
- [ ] Each per-analysis subdirectory has its own README pointing at the Synapse folder + schema + the analysis walkthrough in `docs/analysis/`.
- [ ] No README is more than ~50 lines — these are pointers, not duplicates of the schemas.
- [ ] [`docs/jamborees/2026_UTSW/README.md`](../README.md) and [`TODO.md`](../TODO.md) Step 5 checkbox marked complete.

## Pointers

| Object | Path |
|---|---|
| Template | [`datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/crispr_pipeline/README.md`](../datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/crispr_pipeline/README.md) |
| Per-dataset metadata | [`reference/experimental_metadata.tsv`](../reference/experimental_metadata.tsv) |
| Synapse path map | [`synapse_paths.tsv`](../synapse_paths.tsv) |
