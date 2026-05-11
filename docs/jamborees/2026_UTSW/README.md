# 2026 UTSW Jamboree

Working folder for preparing data and documentation ahead of the 2026 UTSW jamboree. This README is the single entry point — it tells you what we're packaging, where everything lives, and how the pieces fit together.

For scientific scope (topics + working groups), see [`TOPICS.md`](TOPICS.md) and [`WORKING_GROUPS.md`](WORKING_GROUPS.md). For open problems and conversation-ready reports per blocker, see [`issues/`](issues/). New participants: start with [`GETTING_STARTED.md`](GETTING_STARTED.md). Today's plan: latest `AGENDA_<date>.md` in this folder.

## What we're packaging

For each of the 5 production datasets:
- **CRISPR pipeline** outputs (dashboard, MuData, perturbo TSVs)
- **cNMF** gene programs (selected k + sweep-as-provenance)
- **Energy distance** perturbation effects (per-target distances + p-values)

Plus cross-dataset **reference data**: IGVF GTF, TF metadata, experimental metadata, guide library.

Everything lives on Synapse under [`syn64423137`](https://www.synapse.org/Synapse:syn64423137) → `2026_UTSW/`. Each output has a JSON schema in [`schemas/`](schemas/) describing its files and columns. Working-group analyses pull from Synapse + the schemas to generate figures and run cross-dataset comparisons.

## Production datasets

| Dataset | Lab | Cell line | Differentiation |
|---------|-----|-----------|-----------------|
| `Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq` | Hon | WTC11 | Cardiomyocyte |
| `Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq` | Huangfu | HUES8 | Definitive endoderm |
| `Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq` | Huangfu | HUES8 | Embryonic stem cell |
| `Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq` | Gersbach | WTC11 | Hepatocyte |
| `Engreitz_WTC11-endothelial-cells_TF-Perturb-seq` | Engreitz | WTC11 | Endothelial |

> _A 6th dataset, Gersbach WTC11 benchmark HTv2 (`Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2`), is mirrored under `2026_UTSW/datasets/` purely as a small testbed for verifying pipeline structure end-to-end. It is **not** part of the production roster and is not enshrined in the schemas._

## Status at a glance

Status as of 2026-05-09. ✅ = on Synapse, canonical layout. ⚠ = on Synapse, partial / non-canonical. ⏳ = pending (data exists somewhere, not yet packaged). ☐ = blocked.

| Dataset | CRISPR pipeline | cNMF | Energy distance |
|---|---|---|---|
| Hon WTC11 Cardiomyocyte | ⚠ [`syn74520421`](https://www.synapse.org/Synapse:syn74520421) — has `dashboard/` + `pipeline_outputs/` but no `pipeline_info/`. Awaiting the rest of the CRISPR outputs from **Weizhou** (Hon team). | ⏳ blocked on full CRISPR bundle | ⏳ blocked on full CRISPR bundle |
| Huangfu HUES8 Definitive Endoderm | ✅ [`syn74834952`](https://www.synapse.org/Synapse:syn74834952) | ⏳ run pending | ✅ [`syn74883327`](https://www.synapse.org/Synapse:syn74883327) ⚠ p-value calibration |
| Huangfu HUES8 Embryonic Stem Cell | ✅ [`syn74835010`](https://www.synapse.org/Synapse:syn74835010) | ⏳ run pending | ✅ [`syn74883475`](https://www.synapse.org/Synapse:syn74883475) ⚠ p-value calibration |
| Gersbach WTC11 Hepatocyte | ⚠ [`syn70518849`](https://www.synapse.org/Synapse:syn70518849) — non-canonical layout. **Sara** to deliver in our format. | ⏳ **Sara** to deliver in our format | ⏳ **Sara** to deliver in our format |
| Engreitz WTC11 Endothelial | ☐ no data on portal | ☐ | ☐ |

**Reference data** (cross-dataset, all on Synapse): TF metadata [`syn74834227`](https://www.synapse.org/Synapse:syn74834227) ✅ • Experimental metadata [`syn74834309`](https://www.synapse.org/Synapse:syn74834309) ✅ • IGVF GTF [`syn74834518`](https://www.synapse.org/Synapse:syn74834518) ✅ • Guide library [`syn74834519`](https://www.synapse.org/Synapse:syn74834519) ✅.

The full mapping (one Synapse path per dataset × output) lives in [`synapse_paths.tsv`](synapse_paths.tsv). Latest dataset-checkpoint snapshot: [`2026_05_09_state.tsv`](2026_05_09_state.tsv) (compare against the original [`2026_05_07_state.png`](2026_05_07_state.png) / [`2026_05_07_state.tsv`](2026_05_07_state.tsv) to see what's moved).

## Where the data lives

| Tier | Location | What's there |
|---|---|---|
| **Local repo** | `tf_perturb_seq/datasets/<dataset>/` | Per-dataset configs, run scripts, simplified summaries, READMEs that point at Synapse. |
| **HPC (UCSD nrnb)** | `aklie@nrnb-login.ucsd.edu:/cellar/users/aklie/projects/tf_perturb_seq` | Most processed outputs and intermediate files. Source of truth for cNMF + energy distance runs. |
| **GCS** | `gs://igvf-pertub-seq-pipeline-data/<dataset>/<YYYY_MM_DD>/outs/<run>/` | CRISPR pipeline outputs (Nextflow target). |
| **Synapse** | [`syn64423137`](https://www.synapse.org/Synapse:syn64423137) → `2026_UTSW/` | **Canonical store for jamboree-distributed artifacts.** |

## Organizing principle: Synapse-as-we-go

We upload each artifact to Synapse **directly from where it lives** (HPC / GCS / IGVF portal) the moment it's ready, rather than staging everything locally first. Then we log the Synapse path in [`synapse_paths.tsv`](synapse_paths.tsv).

Practical rules:
- **This local folder** holds only docs + small simplified outputs (e.g. `*_simplified.tsv` in [`reference/`](reference/)). No bulky artifacts.
- **The dataset analysis folders** (`datasets/<name>/<analysis>/`) are mostly READMEs that describe what's on Synapse and link to it.
- **Filename convention**: `<name>.tsv` is the comprehensive (machine-readable) form on Synapse; `<name>_simplified.tsv` is the human-readable summary that lives in the repo.

### Folder layout

```
2026_UTSW/                            (mirrored on Synapse + this local folder)
├── README.md                         (this file)
├── TOPICS.md, WORKING_GROUPS.md      (jamboree scope)
├── TODO.md                           (step-by-step prep plan)
├── 2026_05_07_state.{png,tsv}        (initial state snapshot, kept as historical reference)
├── 2026_05_09_state.tsv              (latest snapshot)
├── synapse_paths.tsv                 (dataset × output → Synapse ID)
├── schemas/                          (one JSON schema per output table)
├── scripts/                          (generation + mirror scripts)
├── reference/                        (cross-dataset simplified TSVs)
├── portal_snapshots/                 (timestamped IGVF portal snapshots)
└── datasets/
    └── <dataset_name>/
        ├── README.md                 (what's on Synapse for this dataset)
        ├── crispr_pipeline/
        ├── cnmf/
        └── energy_distance/
```

## Outputs

Per-output, column-level documentation lives in [`schemas/`](schemas/) (see [`schemas/README.md`](schemas/README.md) for the index). Each subsection below points at its schema and gives the operational details.

### Reference data

#### IGVF GTF

Reference gene annotation used across all 5 production datasets. Single canonical version, no subsetting.

| | |
|---|---|
| Filename | `IGVFFI9573KOZR.gtf.gz` (54 MB) |
| Local | `reference/IGVFFI9573KOZR.gtf.gz` |
| HPC | `/cellar/users/aklie/projects/tf_perturb_seq/ref/genome/IGVFFI9573KOZR.gtf.gz` |
| IGVF portal | https://data.igvf.org/reference-files/IGVFFI9573KOZR/ |
| Synapse | [`syn74834518`](https://www.synapse.org/Synapse:syn74834518) |
| Schema | _GTF (external standard); not in `schemas/`_ |

#### TF metadata

One row per unique TF target gene. Joins `target_genes.tsv` with the IGVF GTF, HGNC complete set, Lambert et al. 2018, and JASPAR CORE (human).

| | |
|---|---|
| Generation | `scripts/generate_tf_metadata.py` |
| Comprehensive | `reference/tf_metadata.tsv` (1,983 × 16) — Synapse [`syn74834227`](https://www.synapse.org/Synapse:syn74834227) |
| Simplified | `reference/tf_metadata_simplified.tsv` (1,983 × 8) |
| Schema | [`schemas/tf_metadata.json`](schemas/tf_metadata.json), [`schemas/tf_metadata_simplified.json`](schemas/tf_metadata_simplified.json) |

Resolution priority (per `gene_symbol`): GTF direct → HGNC alias → harmonized guide file fallback. Latest run: 1,983 unique target genes, all resolve to an Ensembl ID; 1,868 in Lambert 2018 (94%); 774 with a JASPAR human entry (39%).

#### Experimental metadata

One row per production dataset. Lab / cell line / differentiation / IGVF accessions / pipeline parameters. Hand-curated; chemistry-related fields auto-extracted from each dataset's `*.config`.

| | |
|---|---|
| Generation | `scripts/generate_experimental_metadata.py` |
| Comprehensive | `reference/experimental_metadata.tsv` (5 × 26) — Synapse [`syn74834309`](https://www.synapse.org/Synapse:syn74834309) |
| Simplified | `reference/experimental_metadata_simplified.tsv` (5 × 12) |
| Schema | [`schemas/experimental_metadata.json`](schemas/experimental_metadata.json), [`schemas/experimental_metadata_simplified.json`](schemas/experimental_metadata_simplified.json) |

Marker conventions: `?` = needs to be filled in; `-` = not applicable; numeric `0` = unknown.

#### Guide library

The IGVF-released TF guide library (pools A-D). Used **as-is** from the portal — no local generation, no simplified version.

| | |
|---|---|
| Filename | `IGVFFI8270UPKB.csv.gz` (340 KB; TSV gzipped despite the `.csv.gz` extension) |
| Local | `reference/IGVFFI8270UPKB.csv.gz` |
| IGVF portal | https://data.igvf.org/tabular-files/IGVFFI8270UPKB/ |
| Synapse | [`syn74834519`](https://www.synapse.org/Synapse:syn74834519) |
| Schema | [`schemas/guide_metadata.json`](schemas/guide_metadata.json) |
| Scope | **Pool A-D only.** Hon CM and Gersbach hepatocyte additionally use a pool F set in the actual experiment; tracked in `experimental_metadata.tsv` `guide_pools` column, not mirrored separately. |

### CRISPR pipeline (`datasets/<dataset>/crispr_pipeline/`)

The IGVF CRISPR FG pipeline's three terminal directories per dataset, mirrored as-is.

| | |
|---|---|
| Schema | [`schemas/crispr_pipeline.json`](schemas/crispr_pipeline.json) |
| Bundle size | ~63 GB / dataset (`pipeline_dashboard/` ~40 GB + `pipeline_outputs/` ~23 GB + `pipeline_info/` ~10 KB) |
| Source | GCS (per-dataset `gcs_output_path` in `experimental_metadata.tsv`); for runs done on HPC, the local-on-HPC variant. |
| Synapse target | `2026_UTSW/datasets/<dataset>/crispr_pipeline/` |
| Mirror script (GCS source) | [`scripts/mirror_pipeline_outputs.py`](scripts/mirror_pipeline_outputs.py) |
| Mirror script (HPC source) | [`scripts/mirror_pipeline_outputs_hpc.py`](scripts/mirror_pipeline_outputs_hpc.py) |

> **Mirroring constraint**: the bundle is too large for the laptop and for HPC `/tmp` (only 20 GB). Run the mirror from the HPC and point `--workdir` at `/cellar/users/aklie/scratch/...` (47 TB free). See [`docs/analysis/CRISPR_PIPELINE_OUTPUTS.md`](../../analysis/CRISPR_PIPELINE_OUTPUTS.md) for the recipe.

### cNMF (`datasets/<dataset>/cnmf/`)

Per-dataset cNMF gene-program-discovery outputs from the torch-cNMF pipeline. Used by Working Groups 1 + 2 (gene programs as primary building blocks; cross-lineage program comparisons; regulators per program).

| | |
|---|---|
| Schema | [`schemas/cnmf.json`](schemas/cnmf.json) |
| Bundle size | ~5–7 GB / dataset |
| Source | HPC: `/cellar/users/aklie/projects/tf_perturb_seq/datasets/<dataset>/PerturbNMF/Result/<run_name>/` |
| Synapse target | `2026_UTSW/datasets/<dataset>/cnmf/<run_name>/` |
| Runner | per-dataset `6_run_cnmf.sh` — see [`docs/analysis/cNMF.md`](../../analysis/cNMF.md) for the runbook |
| Mirror script | [`scripts/mirror_cnmf_outputs.py`](scripts/mirror_cnmf_outputs.py) (HPC → Synapse; implements the `schemas/cnmf.json` curation rule). Takes `--selected-k`. |

**Curation rule** (full schema in [`schemas/cnmf.json`](schemas/cnmf.json)):

1. **Selected-k full data** for downstream analysis — integrated MuData, all loading variants, cell usages, full `Eval/<sel>_2_0/`, selected-k `Plot/`, `Annotation/`, `Interpretation/`.
2. **Sweep-as-provenance** so the k decision is auditable without re-running cNMF — `k_selection.png` + stats, all-k clustering pngs, all-k `gene_spectra_score`, all-k `Eval/` TXT bundles, k-selection figure folder, `README.txt` with the selection rationale.

The selected k is decided as a group ("clinical review board" style per [`docs/analysis/cNMF.md`](../../analysis/cNMF.md)). Per-k MuDatas + per-k usages + intermediate caches stay on HPC.

### Energy distance (`datasets/<dataset>/energy_distance/`)

Per-target energy distances + permutation p-values + diagnostic plots, from steps 1, 2, and 2.1 of [`Chikara-Takeuchi/energy_dist_pipeline`](https://github.com/Chikara-Takeuchi/energy_dist_pipeline). Step 3 (target-by-target matrix + clustering) runs separately after picking cutoffs in `config3.json`.

| | |
|---|---|
| Schema | [`schemas/energy_distance.json`](schemas/energy_distance.json) |
| Bundle size | a few hundred MB / dataset (tables + configs + small PDFs; `inference_mudata.h5mu` is **not** included — already mirrored under `crispr_pipeline/`) |
| Source | HPC: `/cellar/users/aklie/projects/tf_perturb_seq/datasets/<dataset>/results/energy_distance/<run_label>/` |
| Synapse target | `2026_UTSW/datasets/<dataset>/energy_distance/` |
| Runner | per-dataset `5_run_energy_distance.sh`, wrapping `scripts/run_energy_distance_pipeline.sh` |
| Mirror script | [`scripts/mirror_edistance_outputs.py`](scripts/mirror_edistance_outputs.py) (HPC → Synapse) |

> **⚠ Known issue (both Huangfu runs, 2026-05-09)**: p-values are anti-conservative — all 100 negative-control targets have `pval_mean=0` despite their distance distribution overlapping the targeting distribution. The permutation null is too tight relative to our distances vs the HTv2 reference (likely cause: PCA computed on all genes vs HVG subset). Use raw `distance_mean` as an effect-size proxy until calibration is fixed; do not threshold on the p-value alone.

## Operations

### Portal snapshots

[`scripts/query_igvf_portal.py`](scripts/query_igvf_portal.py) snapshots the TF Perturb-seq Project on the IGVF portal (excluding Benchmark) into `portal_snapshots/<utc-iso>/`. Each snapshot writes per-type JSON + TSV (MeasurementSet, AuxiliarySet, AnalysisSet, ConstructLibrarySet) and a `manifest.tsv`. A `portal_snapshots/latest` symlink always points to the newest.

```bash
.venv/bin/python docs/jamborees/2026_UTSW/scripts/query_igvf_portal.py
```

For a daily cron at 06:00 UTC:
```cron
0 6 * * * cd /Users/adamklie/Desktop/tfp3/tf_perturb_seq && .venv/bin/python docs/jamborees/2026_UTSW/scripts/query_igvf_portal.py >> docs/jamborees/2026_UTSW/portal_snapshots/cron.log 2>&1
```

> As of the first snapshot (2026-05-07), only Hon measurement sets appear under the `TF Perturb-seq Project` collection on the portal. Gersbach hepatocyte data is on the portal under a different collection; Engreitz endothelial isn't on the portal yet.

### Mirror scripts

| Script | Direction | What it does |
|---|---|---|
| [`scripts/mirror_pipeline_outputs.py`](scripts/mirror_pipeline_outputs.py) | GCS → Synapse | Downloads `pipeline_dashboard/` + `pipeline_info/` + `pipeline_outputs/` from GCS, uploads to `2026_UTSW/datasets/<dataset>/crispr_pipeline/`, records folder ID in `synapse_paths.tsv`. |
| [`scripts/mirror_pipeline_outputs_hpc.py`](scripts/mirror_pipeline_outputs_hpc.py) | HPC → Synapse | Same target layout but reads from a local HPC run dir (used for runs that didn't land on GCS). |
| [`scripts/mirror_edistance_outputs.py`](scripts/mirror_edistance_outputs.py) | HPC → Synapse | Uploads only deliverables (CSVs, configs, `image/`, `logs/`) — skips intermediates and the input MuData. |
| [`scripts/mirror_cnmf_outputs.py`](scripts/mirror_cnmf_outputs.py) | HPC → Synapse | Uploads the curated cNMF bundle (selected-k full data + sweep-as-provenance) per the `schemas/cnmf.json` rule. Records folder ID in `synapse_paths.tsv`'s `cnmf` column. |

### Generation scripts

| Script | What it builds |
|---|---|
| [`scripts/generate_tf_metadata.py`](scripts/generate_tf_metadata.py) | Comprehensive + simplified TF metadata tables. |
| [`scripts/generate_experimental_metadata.py`](scripts/generate_experimental_metadata.py) | Comprehensive + simplified experimental metadata tables. |
| [`scripts/upload_to_synapse.py`](scripts/upload_to_synapse.py) | Idempotent file uploader to the Synapse mirror. |
| [`../../../src/tf_perturb_seq/crispr_pipeline/cross_dataset_pipeline_summary.py`](../../../src/tf_perturb_seq/crispr_pipeline/cross_dataset_pipeline_summary.py) | Aggregates per-dataset CRISPR pipeline metrics from Synapse → [`reference/cross_dataset_pipeline_summary.tsv`](reference/cross_dataset_pipeline_summary.tsv) (cell counts, UMI medians, knockdown stats, perturbo-significant counts). For WG1 data summarization. |
| [`../../../src/tf_perturb_seq/edistance/cross_dataset_edistance_summary.py`](../../../src/tf_perturb_seq/edistance/cross_dataset_edistance_summary.py) | Aggregates per-dataset energy-distance results from Synapse → [`reference/cross_dataset_edistance_summary.tsv`](reference/cross_dataset_edistance_summary.tsv) (n_targets per type, distance medians, calibration-robust "targeting > NC max" counts, p-value diagnostics). For WG1 transcriptome-wide-significance summary. |

### Validators

| Script | Validates |
|---|---|
| [`../../../src/tf_perturb_seq/edistance/validate_edistance_outputs.py`](../../../src/tf_perturb_seq/edistance/validate_edistance_outputs.py) | Energy distance run dir against `schemas/energy_distance.json` (4 layers: presence, schema, value sanity, cross-ref vs HTv2 reference). |
| [`../../../src/tf_perturb_seq/cnmf/validate_cnmf_outputs.py`](../../../src/tf_perturb_seq/cnmf/validate_cnmf_outputs.py) | cNMF run dir against `schemas/cnmf.json` (4 layers: presence, table shapes, value sanity, cross-ref vs Hon benchmark). Takes `--selected-k`. |

## Working groups

The scientific scope is in [`TOPICS.md`](TOPICS.md) (4 topics from primary building blocks → biological questions → modeling → catalog viz) and [`WORKING_GROUPS.md`](WORKING_GROUPS.md) (5 + 1 working groups, each tied to a topic + figure). Every output here is shaped to answer questions in those docs.

## Tracking

Work is tracked on the [TFP3 GitHub Project](https://github.com/users/adamklie/projects/4) — milestones as draft items, tasks as issues. The [`TODO.md`](TODO.md) in this folder is the step-by-step prep plan.
