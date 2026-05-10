# Issue 6 — Per-dataset READMEs

**Status**: ⏳ Small task. One README per production dataset under `datasets/<dataset>/` so collaborators landing in any folder know what's on Synapse without reading the master README.

**Owner of fix**: us. ~30 min once the production CRISPR bundles are settled.

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

- [ ] `Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq` (gated on Issue #2)
- [ ] `Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq` (CRISPR + ED ready to write; cNMF when runs land)
- [ ] `Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq` (CRISPR + ED ready to write; cNMF when runs land)
- [ ] `Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq` (gated on Issue #3)
- [ ] `Engreitz_WTC11-endothelial-cells_TF-Perturb-seq` (gated on Issue #4)

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
