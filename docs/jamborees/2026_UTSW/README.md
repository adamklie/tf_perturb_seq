# 2026 UTSW Jamboree

Working folder for preparing data and documentation ahead of the 2026 UTSW jamboree.

## Scope

Five **production** datasets are the focus of this jamboree. Bridge / benchmark datasets are tracked elsewhere and are out of scope here.

| Dataset | Lab | Cell line | Differentiation |
|---------|-----|-----------|-----------------|
| `Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq` | Hon | WTC11 | Cardiomyocyte |
| `Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq` | Huangfu | HUES8 | Definitive endoderm |
| `Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq` | Huangfu | HUES8 | Embryonic stem cell |
| `Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq` | Gersbach | WTC11 | Hepatocyte |
| `Engreitz_WTC11-endothelial-cells_TF-Perturb-seq` | Engreitz | WTC11 | Endothelial |

See `2026_05_07_state.tsv` for per-dataset processing status as of 2026-05-07. The original snapshot lives in `2026_05_07_state.png` (which also shows bridge/benchmark datasets — those are not part of this jamboree's scope).

## State columns

The TSV captures four pipeline checkpoints per dataset:

| Column | Meaning |
|--------|---------|
| `igvf_uploaded` | Raw data uploaded to the IGVF portal |
| `crispr_pipeline` | CRISPR perturb-seq pipeline run completed |
| `cnmf_programs` | cNMF gene programs calculated |
| `energy_distance` | Energy distance calculations completed |

Values: `yes` / `no` / `running` / `?` (unknown).

## Files in this folder

| File | Purpose |
|------|---------|
| `TODO.md` | Step-by-step plan for jamboree prep |
| `WORKING_GROUPs.md` | Hand-written outline of the 5 (+1 optional) working groups + their goals |
| `2026_05_07_state.png` | Source-of-truth snapshot of dataset processing state |
| `2026_05_07_state.tsv` | Machine-readable version of the snapshot |
| `synapse_paths.tsv` | Map of (dataset × output type) → Synapse ID, populated as uploads happen |
| `schemas/` | JSON schemas (one per output table) describing columns, types, sources, and notes |
| `scripts/generate_tf_metadata.py` | Builds the comprehensive + simplified TF metadata tables |
| `scripts/generate_experimental_metadata.py` | Builds the comprehensive + simplified experimental metadata tables |
| `scripts/query_igvf_portal.py` | Snapshots IGVF portal state (production-project filesets) into `portal_snapshots/<utc-iso>/` |
| `scripts/upload_to_synapse.py` | Idempotent file uploader to the Synapse mirror |
| `portal_snapshots/` | History of IGVF portal state — one timestamped subdir per snapshot, with raw JSON + per-type TSV summaries |
| `reference/tf_metadata_simplified.tsv` | Simplified human-readable TF metadata |
| `README.md` | This file — extended as decisions are made |

## Portal snapshots

`scripts/query_igvf_portal.py` snapshots the current state of the TF Perturb-seq Project on the IGVF portal (excluding Benchmark) into `portal_snapshots/<utc-iso>/`. Each snapshot writes per-type JSON + TSV (MeasurementSet, AuxiliarySet, AnalysisSet, ConstructLibrarySet) and a `manifest.tsv` summary. A `portal_snapshots/latest` symlink always points to the newest snapshot.

To take a fresh snapshot manually:

```bash
.venv/bin/python docs/jamborees/2026_UTSW/scripts/query_igvf_portal.py
```

To run on a recurring schedule, add a cron entry (e.g., daily at 06:00 UTC):

```cron
0 6 * * * cd /Users/adamklie/Desktop/tfp3/tf_perturb_seq && .venv/bin/python docs/jamborees/2026_UTSW/scripts/query_igvf_portal.py >> docs/jamborees/2026_UTSW/portal_snapshots/cron.log 2>&1
```

**Note** (as of the first snapshot, 2026-05-07): only **Hon** measurement sets currently appear under the `TF Perturb-seq Project` collection on the portal (Gersbach hepatocyte data is on the portal under a different collection; Engreitz endothelial isn't on the portal at all yet).

## Synapse mirror

Top-level Synapse parent: [`syn64423137`](https://www.synapse.org/Synapse:syn64423137) (the `tf_perturb_seq` project folder).

Created for this jamboree:

| Synapse | Type | Path |
|---|---|---|
| [`syn74834225`](https://www.synapse.org/Synapse:syn74834225) | Folder | `2026_UTSW/` |
| [`syn74834226`](https://www.synapse.org/Synapse:syn74834226) | Folder | `2026_UTSW/reference/` |
| [`syn74834227`](https://www.synapse.org/Synapse:syn74834227) | File | `2026_UTSW/reference/tf_metadata.tsv` |
| [`syn74834309`](https://www.synapse.org/Synapse:syn74834309) | File | `2026_UTSW/reference/experimental_metadata.tsv` |
| [`syn74834518`](https://www.synapse.org/Synapse:syn74834518) | File | `2026_UTSW/reference/IGVFFI9573KOZR.gtf.gz` |
| [`syn74834519`](https://www.synapse.org/Synapse:syn74834519) | File | `2026_UTSW/reference/IGVFFI8270UPKB.csv.gz` |
| [`syn74834952`](https://www.synapse.org/Synapse:syn74834952) | Folder | `2026_UTSW/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/crispr_pipeline/` (full bundle: pipeline_dashboard + pipeline_info + pipeline_outputs) |
| [`syn74835010`](https://www.synapse.org/Synapse:syn74835010) | Folder | `2026_UTSW/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/crispr_pipeline/` (full bundle) |

## Where the data lives

- **Local repo:** `tf_perturb_seq/datasets/<dataset_name>/` holds per-dataset configs and metadata.
- **HPC (UCSD nrnb):** `aklie@nrnb-login.ucsd.edu:/cellar/users/aklie/projects/tf_perturb_seq` — most processed outputs and intermediate files.
- **GCS:** `gs://igvf-pertub-seq-pipeline-data/<dataset_name>/YYYY_MM_DD/outs/` — CRISPR pipeline outputs.
- **Synapse:** **canonical store for jamboree-distributed artifacts.** We upload **directly from where the data lives** (HPC / GCS / IGVF portal) to Synapse as each artifact is ready — no full local mirror.

## Strategy: Synapse-as-we-go

Rather than staging everything locally and then mirroring to Synapse at the end, we upload each artifact to Synapse as soon as it's ready, transferring directly from its source (HPC, GCS, IGVF portal). The final deliverable is `synapse_paths.tsv` plus the simplified / human-readable outputs that are small enough to live in this repo.

Practical implications:
- This local folder holds **only**: docs, the state TSV, the synapse path TSV, and small simplified outputs (e.g., `*_simplified.tsv` files in `reference/`).
- Bulky artifacts (MuData, GTFs, cNMF outputs, kallisto indexes) are **not** copied locally — they go straight to Synapse.
- The dataset-level analysis folders (`datasets/<name>/<analysis>/`) primarily hold READMEs that describe what's on Synapse and link to it.

## Folder layout

```
2026_UTSW/
├── README.md
├── TODO.md
├── 2026_05_07_state.{png,tsv}
├── synapse_paths.tsv               # populated as artifacts are uploaded
├── reference/                      # cross-dataset reference data
│   ├── tf_metadata.tsv (or simplified-only — the full table goes to Synapse)
│   ├── tf_metadata_simplified.tsv
│   └── ...
└── datasets/
    └── <dataset_name>/
        ├── README.md               # what's on Synapse for this dataset, with links
        ├── crispr_pipeline/        # README + any small simplified summaries
        ├── cnmf/
        └── energy_distance/
```

Each analysis directory (`crispr_pipeline/`, `cnmf/`, `energy_distance/`) may grow its own subdirectories as needed (e.g., per-run, per-version, per-replicate) — defined per-analysis as we go.

Within each subdir, **simplified** (human-readable) artifacts can live in this repo; the corresponding **detailed** (machine-readable) artifacts are on Synapse. Filename convention: `<name>_simplified.tsv` for the simplified version.

## Outputs to package for the jamboree

Per-table column-level documentation lives in `schemas/` (one JSON file per output, plus a `schemas/README.md` index). The README below only summarizes each output and points to its schema file. For files distributed on the IGVF portal, we prefer linking to the portal URL over storing a local copy.

### Reference

#### IGVF GTF

| | |
|---|---|
| Filename | `IGVFFI9573KOZR.gtf.gz` |
| Local | `reference/IGVFFI9573KOZR.gtf.gz` (54 MB) |
| HPC | `/cellar/users/aklie/projects/tf_perturb_seq/ref/IGVFFI9573KOZR.gtf.gz` |
| IGVF portal | https://data.igvf.org/reference-files/IGVFFI9573KOZR/ |
| Synapse | [`syn74834518`](https://www.synapse.org/Synapse:syn74834518) |
| Schema | _GTF format (standard); not in `schemas/` since it's an external standard._ |
| Description | Reference gene annotation used across all 5 production datasets. Single canonical version, no subsetting. |

#### TF metadata

One row per unique TF target gene. Joins `target_genes.tsv` with the IGVF GTF, HGNC complete set, Lambert et al. 2018, and JASPAR CORE (human).

| | |
|---|---|
| Generation script | `scripts/generate_tf_metadata.py` |
| Comprehensive output | `reference/tf_metadata.tsv` (1,983 × 16) — also on Synapse [`syn74834227`](https://www.synapse.org/Synapse:syn74834227) |
| Simplified output | `reference/tf_metadata_simplified.tsv` (1,983 × 8) |
| Schema | [`schemas/tf_metadata.json`](schemas/tf_metadata.json), [`schemas/tf_metadata_simplified.json`](schemas/tf_metadata_simplified.json) |

**Resolution priority (per `gene_symbol`)**: GTF direct → HGNC alias → harmonized guide file fallback.

**Coverage (latest run)**: 1,983 unique target genes; all resolve to an Ensembl ID (1,946 GTF direct, 37 HGNC alias, 0 harmonized fallback). 1,868 in Lambert 2018 (94%); 774 with a JASPAR human entry (39%).

#### Experimental metadata

One row per production dataset. Lab / cell line / differentiation / IGVF accessions / pipeline parameters. Hand-curated from per-dataset documentation; chemistry-related fields are auto-extracted from each dataset's `*.config` file.

| | |
|---|---|
| Generation script | `scripts/generate_experimental_metadata.py` |
| Comprehensive output | `reference/experimental_metadata.tsv` (5 × 26) — also on Synapse [`syn74834309`](https://www.synapse.org/Synapse:syn74834309) |
| Simplified output | `reference/experimental_metadata_simplified.tsv` (5 × 12) |
| Schema | [`schemas/experimental_metadata.json`](schemas/experimental_metadata.json), [`schemas/experimental_metadata_simplified.json`](schemas/experimental_metadata_simplified.json) |

Marker conventions in the table: `?` = needs to be filled in; `-` = not applicable; numeric `0` = unknown.

#### Guide metadata

The IGVF-released TF guide library (pools A-D), used **as-is** from the portal — no local generation, no simplified version. This file is what gets fed into all downstream guide-assignment / inference steps.

| | |
|---|---|
| Filename | `IGVFFI8270UPKB.csv.gz` |
| Local | `reference/IGVFFI8270UPKB.csv.gz` (340 KB) |
| IGVF portal | https://data.igvf.org/tabular-files/IGVFFI8270UPKB/ |
| Synapse | [`syn74834519`](https://www.synapse.org/Synapse:syn74834519) |
| Format | TSV (gzipped, despite the `.csv.gz` extension), 14,150 rows × 18 cols |
| Schema | [`schemas/guide_metadata.json`](schemas/guide_metadata.json) |
| Scope | **Pool A-D only.** Hon CM and Gersbach hepatocyte additionally use a pool F set in the actual experiment; that's tracked in the `guide_pools` column of `experimental_metadata.tsv` but not mirrored as a separate file. |

### CRISPR pipeline (`datasets/<dataset_id>/crispr_pipeline/`)

The IGVF CRISPR FG pipeline's three terminal directories per dataset, mirrored as-is:

- `pipeline_dashboard/` — dashboard.html, inference_mudata.h5mu, additional_qc/, evaluation_output/, figures/ (~40 GB / dataset)
- `pipeline_info/` — params JSON + software-versions YAML (tiny)
- `pipeline_outputs/` — final inference_mudata.h5mu + perturbo cis/trans per-element/per-guide TSVs (~23 GB / dataset)

| | |
|---|---|
| Schema | [`schemas/crispr_pipeline.json`](schemas/crispr_pipeline.json) |
| Mirror script | `scripts/mirror_pipeline_outputs.py` (GCS → Synapse, uses `gcloud storage rsync` + synapseclient) |
| Bundle size | ~63 GB per dataset |
| Source | GCS (per-dataset `gcs_output_path` in `experimental_metadata.tsv`) |
| Synapse target | `2026_UTSW/datasets/<dataset_id>/crispr_pipeline/` |

**Per-dataset status (2026-05-07):**

| Dataset | Canonical run | Status |
|---|---|---|
| Hon WTC11 Cardiomyocyte | `2026_04_19_no_spacer` (Synapse name; differs from `initial_run` in metadata) | Already on Synapse [`syn74520421`](https://www.synapse.org/Synapse:syn74520421) — has `dashboard/` + `pipeline_outputs/` but **no `pipeline_info/`**. Awaiting full bundle from Hon team. |
| Huangfu HUES8 Definitive Endoderm | `muddy_penguin` | ✅ Mirrored to [`syn74834952`](https://www.synapse.org/Synapse:syn74834952) (full 3-folder bundle, 2026-05-07) |
| Huangfu HUES8 Embryonic Stem Cell | `sceptre_v1` | ✅ Mirrored to [`syn74835010`](https://www.synapse.org/Synapse:syn74835010) (full 3-folder bundle, 2026-05-07) |
| Gersbach WTC11 Hepatocyte | _various_ | Already on Synapse [`syn70518849`](https://www.synapse.org/Synapse:syn70518849) but **not in canonical 3-folder layout** — has `Perturbo_outputs/`, `cNMF_inputs/`, multiple MuData files. Awaiting full bundle from Gersbach team. |
| Engreitz WTC11 Endothelial | — | No data yet |

**Mirroring constraint**: per-bundle size (~63 GB) exceeds local free disk on the laptop, **and HPC `/tmp` is only 20 GB**. The recommended path is to run `scripts/mirror_pipeline_outputs.py` from the HPC and explicitly point `--workdir` at the `/cellar` filesystem (47 TB free):

```bash
# from the HPC, in the project root, for each Huangfu dataset:
gcloud config set account adamklie@igvf-pertub-seq-pipeline.iam.gserviceaccount.com
.venv/bin/python docs/jamborees/2026_UTSW/scripts/mirror_pipeline_outputs.py \
  --dataset Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq \
  --workdir /cellar/users/aklie/scratch/jamboree_pipeline_staging/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq
```

Run the long upload via `nohup ... &` so it survives the SSH session. The script downloads the three pipeline directories from GCS, uploads them under `2026_UTSW/datasets/<dataset_id>/crispr_pipeline/` on Synapse, and records the resulting Synapse folder ID in `synapse_paths.tsv`'s `crispr_pipeline` column. After both runs complete, rsync the updated `synapse_paths.tsv` back to the laptop.

### cNMF (`datasets/<name>/cnmf/`)

_TBD — see step 3._

### Energy distance (`datasets/<dataset_id>/energy_distance/`)

Per-target energy distances + permutation p-values, plus diagnostic plots, from running steps 1, 2, and 2.1 of [`Chikara-Takeuchi/energy_dist_pipeline`](https://github.com/Chikara-Takeuchi/energy_dist_pipeline) on each dataset's inference MuData. Step 3 (target-by-target matrix + clustering) is run separately after picking cutoffs in `config3.json`.

| | |
|---|---|
| Schema | [`schemas/energy_distance.json`](schemas/energy_distance.json) |
| Mirror script | `scripts/mirror_edistance_outputs.py` (HPC → Synapse, uses synapseclient) |
| Bundle size | ~few hundred MB per dataset (tables + configs + small PDFs; `inference_mudata.h5mu` is **not** included — already mirrored under `crispr_pipeline/`) |
| Source | HPC: `/cellar/users/aklie/projects/tf_perturb_seq/datasets/<dataset_id>/results/energy_distance/<run_label>/` |
| Synapse target | `2026_UTSW/datasets/<dataset_id>/energy_distance/` |
| Runner | `scripts/run_energy_distance_pipeline.sh` (wrapped per-dataset by `5_run_energy_distance.sh`) |

**Per-dataset status (2026-05-09):**

| Dataset | Run label | Status | Synapse |
|---|---|---|---|
| Hon WTC11 Cardiomyocyte | TBD | Not run — awaiting full `crispr_pipeline/` bundle from Hon team; source MuData [`syn74522725`](https://www.synapse.org/Synapse:syn74522725) | — |
| Huangfu HUES8 Definitive Endoderm | `muddy_penguin` | ✅ Complete (9h 34m) — all validation layers PASS; ⚠ p-values mis-calibrated (see per-dataset README) | [`syn74883327`](https://www.synapse.org/Synapse:syn74883327) |
| Huangfu HUES8 Embryonic Stem Cell | `sceptre_v1` | ✅ Complete (10h 0m) — all validation layers PASS; ⚠ same calibration concern as DE | [`syn74883475`](https://www.synapse.org/Synapse:syn74883475) |
| Gersbach WTC11 Hepatocyte | TBD | Not run — awaiting canonical run; source MuData [`syn74728027`](https://www.synapse.org/Synapse:syn74728027) | — |
| Engreitz WTC11 Endothelial | — | Blocked — no inference MuData (not on portal yet) | — |

**⚠ Known issue across both Huangfu runs**: p-values are anti-conservative — all 100 negative-control targets have `pval_mean=0` despite negative-control distance distribution overlapping the targeting distribution. The pipeline's permutation null is too tight relative to our ~1-2-orders-of-magnitude-larger distances vs the HTv2 verified reference (likely cause: we kept all genes for PCA, vs HTv2 which used HVG subset). Use raw `distance_mean` as an effect-size proxy until calibration is fixed; do not threshold on the p-value alone.

**Mirror to Synapse (after each run completes):**

```bash
# from the HPC, in the project root, for each completed run:
.venv/bin/python docs/jamborees/2026_UTSW/scripts/mirror_edistance_outputs.py \
  --dataset Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq \
  --source-dir /cellar/users/aklie/projects/tf_perturb_seq/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/results/energy_distance/muddy_penguin
```

The script uploads only the deliverable artifacts (tables, configs, `image/`) — not intermediates (`preprocessed.h5ad`, `gRNA_dict.pickle`, `pca_dataframe.pickle`), the cloned `energy_dist_pipeline/` source, or the downloaded `inference_mudata.h5mu`. After running on the HPC, rsync the updated `synapse_paths.tsv` back to the laptop.

## Tracking

Work for this jamboree is tracked on the [TFP3 GitHub Project](https://github.com/users/adamklie/projects/4) (milestones as drafts, tasks as issues).
