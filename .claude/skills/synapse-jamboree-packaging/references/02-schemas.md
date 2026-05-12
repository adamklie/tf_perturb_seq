# Schemas + canonical layouts

Every bundle uploaded to Synapse must conform to a JSON schema in `docs/jamborees/2026_UTSW/schemas/`. Schemas describe required directories, file names, file types, and (for TSVs) column lists.

## Schema files

```
docs/jamborees/2026_UTSW/schemas/
├── README.md                            (overview of the schema collection)
├── crispr_pipeline.json                 per-dataset
├── cnmf.json                            per-dataset
├── energy_distance.json                 per-dataset
├── experimental_metadata.json           cross-dataset (5 rows)
├── experimental_metadata_simplified.json    same, simplified columns
├── guide_metadata.json                  cross-dataset
├── tf_metadata.json                     cross-dataset (HOCOMOCO / Lambert)
└── tf_metadata_simplified.json          same, simplified columns
```

Each schema has a header block:

```json
{
  "name": "<bundle name>",
  "title": "<human title>",
  "scope": "per-dataset" | "cross-dataset",
  "description": "<what's in this bundle>",
  "bundle_path_local": "datasets/<dataset_id>/...",
  "bundle_path_synapse": "2026_UTSW/datasets/<dataset_id>/...",
  "size_per_dataset_gb_approx": <number>,
  "directories": [...]
}
```

## CRISPR pipeline schema highlights

`crispr_pipeline.json` — ~63 GB per dataset.

**Top-level directories:**
- `pipeline_dashboard/` (~40 GB)
- `pipeline_info/` (~10 KB)
- `pipeline_outputs/` (~20 GB)

**Required files (in `pipeline_dashboard/`):**
- `dashboard.html`
- `inference_mudata.h5mu`
- `additional_qc/{gene,guide,intended_target,trans}_metrics.tsv`
- `evaluation_output/` (bedgraph/bedpe for cis/trans perturbo + sceptre)
- `figures/` (knee plot, scatter, violin, evaluation barplots, precision-recall, volcano)
- `guide_seqSpec_plots/` (per-sample)
- `svg/` (vector versions of dashboard figures)

**Required files (in `pipeline_info/`):**
- `params_<timestamp>.json`
- `nf_core_pipeline_software_versions.yml`

**Required files (in `pipeline_outputs/`):**
- `inference_mudata.h5mu`
- `perturbo_cis_per_element_output.tsv.gz`
- `perturbo_cis_per_guide_output.tsv.gz`
- `perturbo_trans_per_element_output.tsv.gz`
- `perturbo_trans_per_guide_output.tsv.gz`
- `sceptre/` (per-sceptre outputs)

## cNMF schema highlights

`cnmf.json` — ~7–8 GB per dataset (after `--upload-only-k<N>-h5mu` filtering).

**Important:** the **Synapse layout is flat** — no `<cnmf_run_name>/` nesting. Local is `<DS>/<RUN>/cnmf/<cnmf_run_name>/Result/...`; Synapse is `<DS>/cnmf/...`.

**Top-level Synapse layout:**
```
datasets/<DATASET_ID>/cnmf/
├── README.md                                  (per-run summary)
├── inference_mudata_cNMF_<K>_2_0.h5mu         (selected K only)
├── Evaluation_<K>_2_0/                        (all per-K evaluation TSVs for selected K)
│   ├── <K>_perturbation_association_results_all.txt
│   ├── <K>_geneset_enrichment.txt
│   ├── <K>_GO_term_enrichment.txt
│   ├── <K>_trait_enrichment.txt
│   ├── <K>_Explained_Variance.txt
│   └── <K>_fake_perturbation_association_results.txt
├── Plot/
│   ├── k_selection/
│   │   ├── K-selection_panel_2.0.png
│   │   └── K-selection_panel_2.0.svg
│   └── Perturb_gene_<K>_2_0/
│       └── <TF>.pdf (one per perturbed target)
└── Interpretation/Summary_table/<K>_2_0/cNMF_<K>_2_0.xlsx
```

## Energy distance schema highlights

`energy_distance.json` — ~10–50 MB per dataset (small).

**Top-level Synapse layout:**
```
datasets/<DATASET_ID>/energy_distance/
├── pval_edist_full.csv                        (headline)
├── target_by_target_matrix.csv                (if step 3 run)
├── edist_embedding_info.csv                   (if step 3 run)
├── targeting_outlier_table.csv
├── non_targeting_outlier_table.csv
├── discordance_gRNA.csv
├── config1_2.json
├── config3.json
└── *.png                                       (step 2.1 figures)
```

**`pval_edist_full.csv` columns:** `intended_target_name`, `edist_mean`, `edist_std`, `pvalue`, plus per-background breakdowns.

## Cross-dataset reference schemas

- **`tf_metadata.json`** — HOCOMOCO + Lambert harmonized TF master list. Columns: `gene_id`, `gene_name`, `tf_family`, `dna_binding_domain`, plus aliases. Synapse: `syn74834227`.
- **`experimental_metadata.json`** — Per-dataset identity + production status. Columns: `dataset_id`, `lab`, `cell_line`, `differentiation`, `igvf_accession`, `gcs_output_path`, `synapse_id`, status flags. Synapse: `syn74834309`.
- **`guide_metadata.json`** — Per-guide harmonized table (the actual library). Synapse: `syn74834519`.
- **`tf_metadata_simplified.json` / `experimental_metadata_simplified.json`** — Same data with fewer columns; the human-readable summaries that live in the repo under `reference/`.

## How to use a schema

### When uploading a new bundle

1. Read the schema. List every required directory and file.
2. Compare against your local source. Note missing/extra entries.
3. If missing entries: investigate. Either the source is incomplete (fix upstream) or the schema is wrong (file an issue + PR).
4. If extra entries: usually fine; mirror scripts include them. The schema is a minimum, not a maximum (unless explicitly marked).

### When auditing existing Synapse contents

```python
import json, os, synapseclient
syn = synapseclient.Synapse(silent=True); syn.login(authToken=os.environ['SYNAPSE_AUTH_TOKEN'])

schema = json.load(open('docs/jamborees/2026_UTSW/schemas/crispr_pipeline.json'))
syn_id = 'synXXXXX'  # the per-dataset crispr_pipeline folder

# Walk Synapse contents and check against schema['directories']
```

A dedicated `validate_synapse_against_schema.py` is on the TODO list but doesn't exist yet.

### When adding a new bundle type

1. Write the schema first under `schemas/<new_bundle>.json`.
2. Write the mirror script under `scripts/mirror_<new_bundle>.py`.
3. Add the column to `synapse_paths.tsv`.
4. File a PR with all three changes together — schema + script + registry.

## Authoritative source

When a script and a schema disagree, **the schema wins**. File an issue against the script.

When the local source and the schema disagree, the schema wins (in spirit) — but check whether the gap is "this dataset isn't fully cooked yet" (fix the source) vs "the schema requires something we don't produce" (revise the schema).
