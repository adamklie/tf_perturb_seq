# `synapse_paths.tsv` — registry of what's on Synapse

`docs/jamborees/2026_UTSW/synapse_paths.tsv` is the single source of truth for "what's on Synapse for each (dataset, output) pair." Tracked in git, edited on every mirror.

## Schema

```
dataset_id	dataset_name	tf_metadata	igvf_gtf	experimental_metadata	guide_metadata	crispr_pipeline	cnmf	energy_distance	qc
```

| Column | Type | Notes |
|---|---|---|
| `dataset_id` | string | Repo `datasets/` folder name, or `_reference_` for the cross-dataset row |
| `dataset_name` | string | Human-readable label (used in headers / READMEs) |
| `tf_metadata` / `igvf_gtf` / `experimental_metadata` / `guide_metadata` | `syn<id>` or empty | Cross-dataset; populated only in the `_reference_` row |
| `crispr_pipeline` / `cnmf` / `energy_distance` / `qc` | `syn<id>` or empty | Per-dataset; one cell per (dataset, output) |

Empty cell = not yet on Synapse. `-` = explicitly not applicable (e.g., Engreitz Endothelial has no upstream data yet).

## Current state (as of 2026-05-12)

```
dataset_id                                                              tf_metadata  igvf_gtf      experimental  guide_metadata  crispr_pipeline  cnmf         energy_distance  qc
_reference_                                                             syn74834227  syn74834518   syn74834309   syn74834519
Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq                                                                            syn74919102                   syn74897350      syn74917453
Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq                                                                  syn74834952       syn74893844  syn74883327      syn74918479
Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq                                                                   syn74835010       syn74893846  syn74883475      syn74918600
Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq                                                                          syn70518849
Engreitz_WTC11-endothelial-cells_TF-Perturb-seq                                                                                   -
```

(Format with `column -t -s $'\t' synapse_paths.tsv` for readability.)

## When to update

After **every** successful mirror upload. Mirror scripts auto-update the registry, but manual edits are needed when:

- Replacing a non-canonical layout (e.g., Gersbach Hepatocyte `syn70518849` → new canonical syn ID)
- Adding a new bundle type (extend the column list)
- Recording a partial upload (set the cell to the partial syn ID; note in a commit message that it's partial)

## Reading the registry programmatically

```python
import pandas as pd
reg = pd.read_csv('docs/jamborees/2026_UTSW/synapse_paths.tsv', sep='\t', dtype=str).fillna('')

# All datasets with a CRISPR pipeline mirror
ready = reg[(reg.crispr_pipeline.str.startswith('syn')) & (reg.dataset_id != '_reference_')]

# All bundles for one dataset
hon = reg[reg.dataset_id == 'Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq'].iloc[0]
print({c: hon[c] for c in ['crispr_pipeline', 'cnmf', 'energy_distance', 'qc']})
```

## Naming convention

- Production dataset rows use the full repo folder name as `dataset_id` (e.g., `Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq`).
- Synapse IDs are recorded as `synXXXXXX` (no `syn:` prefix, no URL).
- A cell can be empty (not on Synapse), `-` (not applicable), or `synXXXXX` (the parent folder for that bundle).
- For partial / non-canonical mirrors, still record the `synXXXXX` of the existing folder; note the partial status in `docs/jamborees/2026_UTSW/README.md` "Status at a glance".

## Cross-references

- Sister file `docs/jamborees/2026_UTSW/2026_05_09_state.tsv` is a one-time snapshot used to track packaging progress over time. Not auto-updated.
- The `_reference_` row points at cross-dataset Synapse folders (TF metadata, GTF, experimental metadata, guide library).
- The `dataset_name` column is also used as the row label in working-group summary scripts (`build_wg1_*.py`) — keep it short and consistent.

## Commit hygiene

When committing a registry update:

- Use commit message `chore(jamboree): record <bundle> mirror for <dataset> (synXXX)`. Recent commits use this style: `docs(jamboree/Huangfu_ESC): QC mirror (syn74918600) + Synapse-complete state`.
- Include both the registry edit AND the README "Status at a glance" update in the same commit so they don't drift.
- If you're updating the README's status table, also update the dataset-specific note (e.g., "blocked on Weizhou", "Sara to deliver"). Stale notes mislead the team.
