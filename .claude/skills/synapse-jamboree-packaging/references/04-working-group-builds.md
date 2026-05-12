# Working-group build scripts

After bundles are on Synapse, a set of `build_wg<N>_<topic>.py` scripts assembles cross-dataset summary TSVs for working-group analyses. These live in `docs/jamborees/2026_UTSW/scripts/` and write to `working_groups/wg<N>_<topic>/`.

## When to run

Each script tolerates partial coverage (datasets without mirrors get empty cells). Re-run a WG script when:

- A new dataset gets a canonical bundle on Synapse.
- An upstream `reference/*.tsv` (TF metadata, experimental metadata, etc.) changes.
- Inputs to a specific script change (e.g., recalibrated e-distance p-values).

## WG1 — Data QC + cross-dataset summaries

`wg1` is the data-quality working group. Outputs feed almost every downstream WG.

### `build_wg1_qc_summary.py`

**Inputs:**
- `reference/cross_dataset_pipeline_summary.tsv` (per-dataset CRISPR pipeline metrics; populated for 3/5 datasets as of 2026-05-12)
- `reference/experimental_metadata_simplified.tsv` (identity + pipeline status; 5/5 datasets)

**Output:** `working_groups/wg1_data_qc/qc_summary.tsv`

**Logic:** left-join metrics onto identity; select a curated set of QC columns; emit one row per production dataset. Rows for datasets without metrics get empty cells with `data_state = pending` / `blocked`.

### `build_wg1_edistance_summary.py`

**Inputs:** per-dataset `pval_edist_full.csv` mirrors on Synapse (under each dataset's `energy_distance/` folder).

**Output:** `working_groups/wg1_data_qc/edistance_summary.tsv` — one row per (dataset × intended_target), columns include `edist_mean`, `pvalue`, BH-adjusted `qvalue`.

**Caveat:** Huangfu DE and Huangfu ESC are mirrored to Synapse with un-calibrated p-values. The build script BH-adjusts on the fly, but if the upstream e-distance step is re-run with calibrated p-values, the build script should re-run too to override the on-the-fly adjustment.

### `build_wg1_significant_tfs.py`

**Inputs:** WG1 e-distance summary + cNMF perturbation association tables (per dataset on Synapse).

**Output:** `working_groups/wg1_data_qc/significant_tfs.tsv` — per (dataset × TF), is it significant by e-distance, by cNMF perturbation association, both?

### `build_wg1_tf_cross_lineage.py`

**Output:** `working_groups/wg1_data_qc/tf_cross_lineage.tsv` — per TF, which lineages it significantly affects. Useful for identifying "lineage-promiscuous" TFs.

### `build_wg1_trans_target_counts.py`

**Output:** `working_groups/wg1_data_qc/trans_target_counts.tsv` — per (dataset × TF), count of significant trans targets at a fixed FDR cutoff. Inputs: per-dataset calibrated trans results from the `deg-calibration` skill output.

## WG3 — Disease + TF convergence

WG3 connects TFs to disease via OpenTargets / HPO.

### `build_wg3_disease_tf_activity.py`

**Inputs:** WG1 cross-lineage TFs + OpenTargets / HPO mappings (`reference/hpo_gene_disease.tsv` from `fetch_hpo_gene_disease.py`).

**Output:** `working_groups/wg3_disease/disease_tf_activity.tsv` — per (disease × TF), evidence score from cross-lineage perturbation activity.

### `build_wg3_tf_convergence_scorecard.py`

**Output:** `working_groups/wg3_disease/tf_convergence_scorecard.tsv` — per TF, summary of cross-lineage convergence (does this TF do the same thing in multiple lineages, or different things?).

## WG4 — Network structure

WG4 builds TF→gene regulatory edges.

### `build_wg4_tf_gene_edges.py`

**Inputs:** per-dataset calibrated DEG results from `deg-calibration` (significant TF→gene pairs).

**Output:** `working_groups/wg4_networks/tf_gene_edges.tsv` — one row per significant (TF, gene, dataset) triple, with effect size + direction.

### `build_wg4_network_structure_by_lineage.py`

**Output:** `working_groups/wg4_networks/network_structure_by_lineage.tsv` — per-lineage network summary stats (hub TFs, mean degree, modularity).

## WG5 — TF family-level

WG5 aggregates by TF family.

### `build_wg5_tf_family_scorecard.py`

**Inputs:** WG1 TF activity tables + `reference/tf_metadata.tsv` (HOCOMOCO+Lambert family labels).

**Output:** `working_groups/wg5_tf_families/tf_family_scorecard.tsv` — per TF family, summary of activity across datasets / lineages.

## Dependency graph

```
reference/*.tsv (TF metadata, exp metadata, harmonized guide library)
    │
    ▼
build_wg1_qc_summary          ← reference/cross_dataset_pipeline_summary.tsv
build_wg1_edistance_summary   ← per-dataset Synapse e-dist mirrors
build_wg1_significant_tfs     ← e-dist + cNMF
build_wg1_tf_cross_lineage    ← significant_tfs (depends on)
build_wg1_trans_target_counts ← per-dataset calibrated DEG (depends on deg-calibration skill outputs)
    │
    ▼
build_wg3_disease_tf_activity ← wg1_tf_cross_lineage + HPO
build_wg3_tf_convergence      ← wg1_tf_cross_lineage
build_wg4_tf_gene_edges       ← per-dataset calibrated DEG
build_wg4_network_structure   ← wg4_tf_gene_edges
build_wg5_tf_family_scorecard ← wg1_significant_tfs + tf_metadata
```

Run WG1 builds first, then WG3/4/5 (they consume WG1 outputs).

## Auth + env

All scripts assume:

- `SYNAPSE_AUTH_TOKEN` in env (run via `zsh -ic '...'` to pull from `~/.zshrc`).
- `.venv` activated (pandas + synapseclient).
- Working directory: repo root (the scripts use `Path(__file__).parent.parent` to resolve `JAMB_ROOT`).

## What to commit

- All WG outputs under `working_groups/wg<N>_<topic>/` are **tracked** (small TSVs). Commit alongside the build-script run.
- The build scripts themselves are also tracked.
- The Synapse-mirrored bulk outputs (per-dataset bundles) are **not** in this repo.

## Building everything from scratch

```bash
cd <REPO_ROOT>
zsh -ic '
source .venv/bin/activate
for wg in build_wg1_qc_summary build_wg1_edistance_summary build_wg1_significant_tfs \
          build_wg1_tf_cross_lineage build_wg1_trans_target_counts; do
  python docs/jamborees/2026_UTSW/scripts/${wg}.py || { echo "FAIL: $wg"; break; }
done

for wg in build_wg3_disease_tf_activity build_wg3_tf_convergence_scorecard \
          build_wg4_tf_gene_edges build_wg4_network_structure_by_lineage \
          build_wg5_tf_family_scorecard; do
  python docs/jamborees/2026_UTSW/scripts/${wg}.py || { echo "FAIL: $wg"; break; }
done
'
```

Failures usually cascade — if WG1 fails, WG3/4/5 will fail too. Read the first error first.

## When working-group analyses fail

| Symptom | Cause | Fix |
|---|---|---|
| `FileNotFoundError: cross_dataset_pipeline_summary.tsv` | Not all CRISPR mirrors are complete | Mirror remaining datasets, or accept partial coverage (WG1 QC summary already handles this) |
| `KeyError: 'qvalue'` in WG1 e-distance | p-value column has different name | Inspect `pval_edist_full.csv` header; older runs may use `pvalue_raw` |
| `synapseclient.core.exceptions.SynapseAuthenticationError` | Token unset | `zsh -ic '...'` |
| Empty output | All datasets in `data_state = pending` | Wait for mirrors |
