# Issue 3 — Gersbach Hep deliverables (CRISPR pipeline + cNMF + energy distance)

**Status**: ⚠ Synapse mirror exists at [`syn70518849`](https://www.synapse.org/Synapse:syn70518849) but in a non-canonical layout that doesn't match our schemas. No production cNMF or energy distance for Gersbach Hep yet — Sara likely has both runs already and we need them shaped to our formats.

**Owner of fix**: **Sara Geraghty** (Gersbach team / Duke).

## TL;DR

Sara has been running Gersbach Hep through the Gersbach lab's tooling and almost certainly has all three outputs (CRISPR pipeline, cNMF, energy distance) in some form. What we need is for those outputs to:

1. Land in our canonical Synapse paths under `2026_UTSW/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/`
2. Match our three schemas in [`schemas/`](../schemas/)
3. Be from a single canonical run (not a mix of variants)

The HTv2 testbed at [`syn74885574`](https://www.synapse.org/Synapse:syn74885574) is the structural example for what each bundle should look like; the Huangfu DE/ESC mirrors at [`syn74834952`](https://www.synapse.org/Synapse:syn74834952) / [`syn74835010`](https://www.synapse.org/Synapse:syn74835010) (CRISPR) and [`syn74883327`](https://www.synapse.org/Synapse:syn74883327) / [`syn74883475`](https://www.synapse.org/Synapse:syn74883475) (energy distance) are the production-equivalent precedents.

## What's currently on Synapse for Gersbach Hep

```
syn70518849                Folder  "Gersbach hepatocyte" (existing pre-jamboree mirror)
├── Perturbo_outputs/      ← perturbo TSVs in some layout
├── cNMF_inputs/           ← intermediate; not what we want
└── (multiple MuData files in non-canonical positions)
```

**This is not the canonical 3-folder CRISPR bundle.** Mixing variants + intermediates makes it ambiguous which "canonical" run any downstream analysis was based on.

We do **not** want to re-derive from this — too messy. Cleaner to ask Sara for a single canonical run.

## What we want — three deliverables

### 3a. CRISPR pipeline bundle

| | |
|---|---|
| Schema | [`schemas/crispr_pipeline.json`](../schemas/crispr_pipeline.json) |
| Required layout | `pipeline_dashboard/` + `pipeline_info/` + `pipeline_outputs/` (the 3 terminal Nextflow output dirs) |
| Synapse target | `2026_UTSW/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/crispr_pipeline/` |
| Reference example | Huangfu DE [`syn74834952`](https://www.synapse.org/Synapse:syn74834952) (what we want it to look like) |

Source can be HPC, Sara's local, or GCS — Adam can mirror from any of these (we have GCS-source and HPC-source variants of the mirror script).

### 3b. cNMF outputs (selected-k bundle + sweep-as-provenance)

| | |
|---|---|
| Schema | [`schemas/cnmf.json`](../schemas/cnmf.json) |
| Curation rule | Selected-k full data (MuData + all loadings + cell usages + full `Eval/<sel>_<dt>/` + Plot/Annotation) **plus** sweep-as-provenance (`k_selection.png` + stats, all-k clustering pngs, all-k `gene_spectra_score`, all-k `Eval/` TXT bundles, k-selection figure folder, `README.txt` with the selection rationale). Drop intermediate `cnmf_tmp/`, `prog_data/`, per-k MuDatas. |
| Synapse target | `2026_UTSW/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/cnmf/<run_name>/` |
| Bundle size | ~5–7 GB (vs ~72 GB raw run dir — the curation rule drops the bulk) |
| Reference example | Hon WTC11 benchmark `030726_20iter_5KHVG_torch_halsvar_batch_e7` run on HPC at `/cellar/users/aklie/projects/tf_perturb_seq/datasets/Hon_WTC11-benchmark_TF-Perturb-seq/PerturbNMF/` (used as the structural reference; not on Synapse but documented in [`docs/analysis/cNMF_OUTPUTS.md`](../../../analysis/cNMF_OUTPUTS.md)) |

The `README.txt` per run is critical — it captures the selected k + the rationale (one-time clinical-review-board decision per [`docs/analysis/cNMF.md`](../../../analysis/cNMF.md)) so a future reader can revisit the k decision without re-running.

If Sara's existing run has different parameters (e.g., different k sweep, different HVG count) than the schema's reference values, that's fine — we'd just like the parameters captured in the run's `README.txt` so we know what we're getting.

### 3c. Energy distance outputs

| | |
|---|---|
| Schema | [`schemas/energy_distance.json`](../schemas/energy_distance.json) |
| Required files | `pval_edist_full.csv`, `targeting_outlier_table.csv`, `non_targeting_outlier_table.csv`, `config1_2.json`, `image/*.pdf`. Optional: `target_by_target_matrix.csv`, `edist_embedding_info.csv`, `config3.json` (step 3 — only if cutoffs were chosen and clustering was run). |
| Synapse target | `2026_UTSW/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/energy_distance/` |
| Reference example | Huangfu DE [`syn74883327`](https://www.synapse.org/Synapse:syn74883327) — has steps 1 + 2 + 2.1 only (no step 3 yet) |

**Calibration concern (see [Issue #1](edistance-calibration.md))**: Sara should be aware that our Huangfu DE/ESC runs had anti-conservative p-values (all NCs `pval_mean=0`). If Sara's Gersbach Hep run was set up similarly (all-genes PCA, large NC pool), it may have the same issue. Worth a sanity check on her end before delivering.

## Practical paths

Sara can either:

- **Path A (lightest)** — Upload to Synapse herself directly to the target paths above. Adam updates `synapse_paths.tsv` to point at her IDs.
- **Path B (Adam's mirror scripts)** — Sara puts the runs at a GCS path or HPC location Adam can access; Adam runs:
  - [`scripts/mirror_pipeline_outputs.py`](../scripts/mirror_pipeline_outputs.py) (GCS source) or [`scripts/mirror_pipeline_outputs_hpc.py`](../scripts/mirror_pipeline_outputs_hpc.py) (HPC source) for the CRISPR bundle.
  - [`scripts/mirror_edistance_outputs.py`](../scripts/mirror_edistance_outputs.py) for energy distance.
  - cNMF mirror script [`scripts/mirror_cnmf_outputs.py`](../scripts/) doesn't exist yet — Adam will write it once Sara confirms file naming, or Sara can upload directly. The schema in [`schemas/cnmf.json`](../schemas/cnmf.json) → `bundle_inclusion_rule` enumerates the full include/exclude list.

Path A is simplest if Sara is comfortable with synapseclient. Path B keeps the bundle assembly automated on our side.

## Acceptance criteria

For each of the three outputs:
- [ ] Synapse folder exists at the canonical path under `2026_UTSW/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/<analysis>/`.
- [ ] Folder contents validate against the schema in [`schemas/<analysis>.json`](../schemas/) (file presence + filename pattern).
- [ ] [`synapse_paths.tsv`](../synapse_paths.tsv) Gersbach Hep row updated with the new Synapse IDs.
- [ ] [`docs/jamborees/2026_UTSW/README.md`](../README.md) at-a-glance row for Gersbach Hep moves from ⚠ to ✅ for that output.

## Pointers

| Object | Path |
|---|---|
| Existing non-canonical mirror | Synapse [`syn70518849`](https://www.synapse.org/Synapse:syn70518849) |
| CRISPR schema | [`schemas/crispr_pipeline.json`](../schemas/crispr_pipeline.json) |
| cNMF schema | [`schemas/cnmf.json`](../schemas/cnmf.json) |
| Energy distance schema | [`schemas/energy_distance.json`](../schemas/energy_distance.json) |
| CRISPR analysis walkthrough | [`docs/analysis/CRISPR_PIPELINE_OUTPUTS.md`](../../../analysis/CRISPR_PIPELINE_OUTPUTS.md) |
| cNMF analysis walkthrough | [`docs/analysis/cNMF_OUTPUTS.md`](../../../analysis/cNMF_OUTPUTS.md) |
| Energy distance analysis walkthrough | [`docs/analysis/ENERGY_DISTANCE_OUTPUTS.md`](../../../analysis/ENERGY_DISTANCE_OUTPUTS.md) |
| Mirror scripts | [`scripts/mirror_pipeline_outputs.py`](../scripts/mirror_pipeline_outputs.py), [`scripts/mirror_pipeline_outputs_hpc.py`](../scripts/mirror_pipeline_outputs_hpc.py), [`scripts/mirror_edistance_outputs.py`](../scripts/mirror_edistance_outputs.py) |
