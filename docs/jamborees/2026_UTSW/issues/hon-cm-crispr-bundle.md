# Issue 2 — Hon CM CRISPR pipeline bundle is partial

**Status**: ⚠ Synapse mirror has `dashboard/` + `pipeline_outputs/` but no `pipeline_info/`. The downstream cNMF + energy distance for Hon CM are gated on closing this gap.

**Owner of fix**: **Weizhou** (Hon team) — needs to deliver the rest of the canonical CRISPR pipeline outputs.

## TL;DR

The schema we agreed on for CRISPR pipeline bundles is the IGVF FG pipeline's three terminal directories:

- `pipeline_dashboard/` — dashboard, MuData, additional QC
- `pipeline_info/` — Nextflow params + software versions
- `pipeline_outputs/` — final perturbo cis/trans TSVs

For Hon WTC11 Cardiomyocyte:

| Path | Source | Status |
|---|---|---|
| Synapse [`syn74520421`](https://www.synapse.org/Synapse:syn74520421) | Hon team mirror, run label `2026_04_19_no_spacer` | Has `dashboard/` + `pipeline_outputs/`, **no `pipeline_info/`** |
| GCS `gs://igvf-pertub-seq-pipeline-data/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/2026_04_15/outs/initial_run/` (canonical path per `experimental_metadata.tsv` `gcs_output_path`) | Adam's GCS run | Only has stages 1-2 outputs (`seqspecparser/`, `mappingscrna/`, `mappingguide/`, `prepare/`, `pipeline_info/`) — pipeline didn't reach dashboard/outputs |
| GCS `gs://igvf-pertub-seq-pipeline-data/scratch/bioinfolucas/Benchmark_cleanser_800_mito_15pc/GaryHonDataset/` | Lucas's benchmarking run (cleanser inference, not perturbo) | Has dashboard + info, **but uses cleanser TSV naming (`cis_per_element_results.tsv.gz`), not perturbo (`perturbo_cis_per_element_output.tsv.gz`)** — different inference pipeline |
| HPC `/cellar/users/aklie/projects/tf_perturb_seq/datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/` | Adam's local | No `runs/` dir — no canonical pipeline output staged locally |

We can't trivially assemble the canonical 3-folder bundle from sources we control. Need the Hon team's canonical run.

## Current Synapse layout at `syn74520421`

```
syn74520421                           Folder  "2026_04_19_no_spacer"
├── dashboard/                        syn74526301  ← (named "dashboard", not "pipeline_dashboard")
└── pipeline_outputs/                 syn74520424
```

That's a 2-of-3 canonical bundle. The folder naming (`dashboard/` vs `pipeline_dashboard/`) is a small layout deviation we can live with or rename later.

## Evidence the gap blocks downstream

- [`docs/jamborees/2026_UTSW/schemas/crispr_pipeline.json`](../schemas/crispr_pipeline.json) requires the 3-folder layout.
- [`docs/jamborees/2026_UTSW/schemas/energy_distance.json`](../schemas/energy_distance.json) → `production_runs_status.Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq` says `not run yet (awaiting full crispr_pipeline bundle from Hon team; source MuData syn74522725)`.
- [`docs/analysis/CRISPR_PIPELINE_OUTPUTS.md`](../../../analysis/CRISPR_PIPELINE_OUTPUTS.md) → `pipeline_info/` (`params_*.json` + `nf_core_pipeline_software_versions.yml`) is the only place we capture the actual run config; without it we can't reproduce or audit the run.

## What we want from Weizhou

Concretely, the three canonical folders for the Hon team's `2026_04_19_no_spacer` run (or whichever they consider canonical), uploaded to Synapse — ideally to `2026_UTSW/datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/crispr_pipeline/` to match the production layout, but adding the missing pieces under the existing `syn74520421` is also fine.

Required:
- [ ] `pipeline_info/` — `params_<timestamp>.json` and `nf_core_pipeline_software_versions.yml` from the canonical Nextflow run.

Optional but helpful:
- [ ] Confirmation that the existing `syn74520421/dashboard/` and `pipeline_outputs/` came from the *same* run as the `pipeline_info/` they upload (i.e., they're consistent, not from different launches).
- [ ] Rename `dashboard/` → `pipeline_dashboard/` to match the schema (small ergonomic fix).

## Reference example — what "canonical bundle" looks like

For a complete example of the layout we want, see [`syn74834952`](https://www.synapse.org/Synapse:syn74834952) (Huangfu DE). It has `pipeline_dashboard/` + `pipeline_info/` + `pipeline_outputs/`, all from the same `muddy_penguin` run.

The mirror script we used: [`scripts/mirror_pipeline_outputs.py`](../scripts/mirror_pipeline_outputs.py) (GCS source) — runs from HPC with `--workdir /cellar/users/aklie/scratch/...` to handle the ~63 GB bundle size.

## Acceptance criteria

- [ ] `pipeline_info/` present alongside `dashboard/` (or `pipeline_dashboard/`) and `pipeline_outputs/` on Synapse for Hon CM.
- [ ] [`synapse_paths.tsv`](../synapse_paths.tsv) `crispr_pipeline` column for Hon CM is updated if the path changes.
- [ ] `experimental_metadata.tsv` `canonical_run_label` (currently `initial_run`) reconciled with whatever Weizhou names the canonical run.
- [ ] [`docs/jamborees/2026_UTSW/README.md`](../README.md) at-a-glance row for Hon CM CRISPR moves from ⚠ to ✅.

## Pointers

| Object | Path |
|---|---|
| Existing partial mirror | Synapse [`syn74520421`](https://www.synapse.org/Synapse:syn74520421) |
| Canonical 3-folder example | Synapse [`syn74834952`](https://www.synapse.org/Synapse:syn74834952) (Huangfu DE — what the layout should look like) |
| Schema | [`schemas/crispr_pipeline.json`](../schemas/crispr_pipeline.json) |
| Mirror script (GCS source) | [`scripts/mirror_pipeline_outputs.py`](../scripts/mirror_pipeline_outputs.py) |
| Mirror script (HPC source) | [`scripts/mirror_pipeline_outputs_hpc.py`](../scripts/mirror_pipeline_outputs_hpc.py) |
| Analysis-level walkthrough | [`docs/analysis/CRISPR_PIPELINE_OUTPUTS.md`](../../../analysis/CRISPR_PIPELINE_OUTPUTS.md) |
| Source MuData | Synapse [`syn74522725`](https://www.synapse.org/Synapse:syn74522725) (Hon's cleaned inference MuData) |
