---
name: dataset-scaffolder
description: Scaffold a new TFP3 dataset directory under datasets/&lt;Lab&gt;_&lt;CellLine&gt;-&lt;cond&gt;_TF-Perturb-seq/ with the canonical layout (setup/{scripts,configs,samplesheets}, README.md, dataset_config.yaml) and the Stage 1 driver scripts pre-filled with the IGVF accession and dataset name. Triggers on keywords like new dataset, scaffold dataset, onboard dataset, create dataset directory, dataset template, dataset_config.
user_invocable: true
---

# Dataset Scaffolder

You are an interactive assistant for onboarding a new dataset into TFP3. This skill creates the canonical `datasets/<Lab>_<CellLine>-<cond>_TF-Perturb-seq/` directory tree with the right subdirectories, a `README.md` adapted from the template, a `dataset_config.yaml` skeleton, and the Stage 1 driver scripts (`1_generate_per_sample_metadata.sh`, `2_upload_to_gcp.sh`, `3_patch_gcp_files.sh`) pre-filled with the accession.

For Stage 2 (`4_run_CRISPR_pipeline.sh`) scaffolding, hand off to `crispr-pipeline-runner`'s scaffolder. For Stages 3/4/5 setup, those skills' run-level scaffolders apply per-run, not per-dataset.

## Constants

```
REPO_ROOT:           /Users/adamklie/Desktop/tfp3/tf_perturb_seq  (local Mac)
TEMPLATE:            docs/data/dataset_template/
SCAFFOLDER:          .claude/skills/dataset-scaffolder/scripts/scaffold_dataset.py
SISTER (Stage 1):    .claude/skills/igvf-portal-staging/scripts/scaffold_setup_scripts.py
SISTER (Stage 2):    .claude/skills/crispr-pipeline-runner/scripts/scaffold_run_script.py
```

## Naming convention

```
<Lab>_<CellLine>-<condition>_TF-Perturb-seq[_<chemistry>]
```

Examples:

| Dataset folder | Lab | Cell line | Condition | Chemistry |
|---|---|---|---|---|
| `Hon_WTC11-benchmark_TF-Perturb-seq` | Hon | WTC11 | benchmark | (single tech) |
| `Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq` | Hon | WTC11 | cardiomyocyte-differentiation | (single tech) |
| `Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq` | Huangfu | HUES8 | definitive-endoderm-differentiation | (single tech) |
| `Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3` | Gersbach | WTC11 | benchmark | GEM-Xv3 |
| `Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2` | Gersbach | WTC11 | benchmark | HTv2 |

Chemistry suffix is required when a lab contributes multiple chemistries for the same condition (Gersbach contributes GEM-Xv3 and HTv2 to the benchmark).

The scaffolder validates against a regex but won't refuse a non-canonical name — it just warns. Per the user's global preference, **ask before reorganizing existing datasets** that don't match.

## Step 0: Gather inputs

Ask (or infer):

1. **Lab** — `Hon`, `Huangfu`, `Engreitz`, `Gersbach`, or other.
2. **Cell line** — `WTC11`, `HUES8`, etc.
3. **Condition** — `benchmark`, `cardiomyocyte-differentiation`, `definitive-endoderm-differentiation`, `embryonic-stemcell-differentiation`, `hepatocyte-differentiation`, `endothelial-cells`, etc.
4. **Chemistry suffix?** — Only if multiple chemistries from same lab/condition.
5. **IGVF accession** — `IGVFDSXXXXXXXX` for the analysis set. Required by `1_generate_per_sample_metadata.sh`.

## Step 1: Run the scaffolder

```bash
python3 .claude/skills/dataset-scaffolder/scripts/scaffold_dataset.py \
  --dataset-name <Lab>_<CellLine>-<cond>_TF-Perturb-seq[_<chem>] \
  --accession IGVFDS<XXXXXXXX> \
  --lab <Lab> \
  --cell-line <CellLine> \
  --condition <cond>
```

What it creates:

```
datasets/<DATASET_NAME>/
├── README.md                  # Adapted from docs/data/dataset_template/README.md, header filled in
├── setup/
│   ├── scripts/
│   │   ├── 1_generate_per_sample_metadata.sh   # ACCESSION filled in
│   │   ├── 2_upload_to_gcp.sh                  # DATASET_NAME filled in
│   │   └── 3_patch_gcp_files.sh                # DATASET_NAME filled in
│   ├── configs/               # empty (per-run Nextflow configs go here)
│   └── samplesheets/          # empty (Stage 1 output lands here)
└── dataset_config.yaml         # Adapted from template; defaults for chemistry/QC params
```

The scaffolder is idempotent — re-running with the same `--dataset-name` skips files that exist. Pass `--force` to overwrite.

The Stage 1 setup scripts are written by delegating to `.claude/skills/igvf-portal-staging/scripts/scaffold_setup_scripts.py`. See that skill's docs for details.

## Step 2: Edit `dataset_config.yaml`

The scaffolder writes a config with TFP3-default values. Edit it to match the dataset's chemistry/protocol:

```yaml
dataset_name: <Lab>_<CellLine>-<cond>_TF-Perturb-seq
igvf_accession: IGVFDS<XXX>

enable_data_hashing: true    # HTO dataset?
is_10x3v3: false              # 10x 3' v3?
spacer_tag: "TTAGCTCTTAAAC"   # CROP-seq default; "" for direct-capture; check guide library
reverse_complement_guides: true
dual_guide: false

qc_min_genes_per_cell: 500
qc_min_cells_per_gene: 0.05
qc_pct_mito: 15               # TFP3 uses 15%, not the template's 20%

guide_assignment_method: "sceptre"   # or "cleanser"
multiplicity_of_infection: "high"
```

For the full param meaning, see `.claude/skills/crispr-pipeline-runner/references/01-config-spec.md`. The `dataset_config.yaml` is a **convenience** — it's mirrored into the Nextflow `.config` by hand or via `crispr-pipeline-runner`'s scaffolder.

## Step 3: Fill in `README.md`

The template README has placeholders (`{DATASET_NAME}`, `{LAB_NAME}`, etc.). The scaffolder fills in identity headers but leaves the per-dataset notes section empty. Add:

- One-paragraph description (lab, cell line, condition, special context).
- Links to relevant GitHub issues / TODO entries.
- Synapse ID(s) if the data is already mirrored.
- Per-run table once the first CRISPR pipeline run completes (see `docs/data/DATA.md` §"Per-dataset layout").

## Step 4: Connect to downstream stages

After scaffolding:

1. **Stage 1 (immediate):** Set IGVF credentials, then run `1_generate_per_sample_metadata.sh` through `3_patch_gcp_files.sh` via the `igvf-portal-staging` skill.
2. **Stage 2 (after Stage 1 completes):** Use `crispr-pipeline-runner`'s scaffolder to generate `4_run_CRISPR_pipeline.sh` + the per-run `.config`:

   ```bash
   python3 .claude/skills/crispr-pipeline-runner/scripts/scaffold_run_script.py \
     --dataset-name <DATASET_NAME> \
     --run-label <RUN_LABEL> \
     --data-date <YYYY_MM_DD>
   ```

3. **Stage 3/4/5 (later):** Per-run scripts go into `<RUN_LABEL>/{qc,energy_distance,cnmf}/Script/` after each CRISPR pipeline run lands.
4. **Documentation:** Add the new dataset row to `docs/data/DATA.md` inventory table.
5. **Jamboree (if production):** Update `docs/jamborees/2026_UTSW/synapse_paths.tsv` with the new `dataset_id` row.

## Important notes

- **Don't reorganize existing datasets.** Per user global guideline, only scaffold *new* datasets unless explicitly asked to migrate.
- **Naming is sticky.** The dataset folder name flows into Synapse paths, GCS prefixes, and many references — bikeshed the name once, then commit.
- **IGVF accession must already exist.** The portal-side analysis set must be properly configured (see `igvf-portal-staging`'s `references/01-portal-query.md` §"Prereqs (portal side)") before Stage 1 will succeed.
- **The `bin/` subdirectory** is for per-dataset script variants (per the frozen-scripts rule in memory `[[feedback_dataset_local_scripts]]`). Not auto-created — make it when needed.
- **Numbered scripts (`1_`, `2_`, `3_`, …) go under `setup/scripts/`** for Stage 1–2 drivers; per-run analysis scripts go under `<RUN_LABEL>/<analysis_tier>/scripts/` (energy distance) or `<RUN_LABEL>/<analysis_tier>/Script/` (cNMF, capital-S convention).
- The full per-dataset layout doc: [docs/data/DATA.md](../../../docs/data/DATA.md) §"Per-dataset layout".
