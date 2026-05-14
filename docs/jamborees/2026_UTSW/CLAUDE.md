# 2026 UTSW Jamboree — Claude working notes

Agent-facing conventions for this folder. Participant-facing onboarding is in [`README.md`](README.md); do not duplicate it here.

Inherits everything from [`../../../CLAUDE.md`](../../../CLAUDE.md) (style guide, daily workflow, `module load git`, `uv` not `pip`). This file only adds what's *jamboree-specific*.

## Scope

- **In scope**: the 5 production datasets only — Hon CM, Huangfu DE, Huangfu ESC, Gersbach Hep, Engreitz Endo.
- Working dir is `docs/jamborees/2026_UTSW/`. Run scripts from this folder root.
- Synapse parent is `syn64423137/2026_UTSW/`.

## Environment

- **Python venv**: `/Users/adamklie/Desktop/tfp3/tf_perturb_seq/.venv` (activate with `source /Users/adamklie/Desktop/tfp3/tf_perturb_seq/.venv/bin/activate`, or just `uv run python …` from the repo root).
- **Synapse auth**: `SYNAPSE_AUTH_TOKEN` is exported in `~/.zshrc`. Use `os.environ["SYNAPSE_AUTH_TOKEN"]` in Python; no need to prompt for it.

## Layout (what goes where)

```
README.md              participant landing page
data/                  per-dataset READMEs + schemas/ + scripts/ (mirror/upload)
reference/             cross-dataset refs (TF/exp metadata, guide library, gene annotations) + scripts/
guide/                 interpretation docs (CRISPR / cNMF / ED / formats / FAQ / glossary)
working_groups/        wg1..wg6 — each with README.md + examples/ (build + examples + roll-up TSVs)
manuscript/            placeholder
```

Per-dataset folders under `data/<dataset>/` hold only `README.md`s that link to Synapse. Bulky outputs live on Synapse, **not** in this folder.

## Conventions

- **Synapse-as-we-go**: upload artifacts to Synapse directly from where they live (HPC / GCS / portal). Do **not** stage bulky outputs in this directory.
- **What lives locally**: docs, JSON schemas, small `_simplified` tables, WG roll-up TSVs sized for slide decks.
- **Naming**: simplified / human-readable variants get a `_simplified` suffix (e.g., `tf_metadata.tsv` + `tf_metadata_simplified.tsv`). Both live side-by-side; no audience-based subdirs.
- **Schemas**: every output table and metadata file has a JSON schema in `data/schemas/`. Update the schema when the table changes.
- **Status vocabulary** (per-dataset and WG READMEs): `ready` / `caveat` / `blocked` / `discussion`. Stick to these.

## Working-group roll-ups

- WG content under `working_groups/wg<N>_*/examples/` is **illustrative, not finished analysis**. The `examples/` subdir name is intentional — TSV roll-ups, `build_*.py`, and `examples.py` are all examples of what a participant might do, not deliverables.
- Each `working_groups/wg<N>_*/examples/build_*.py` re-runs an example cross-dataset roll-up off the Synapse-mirrored sources.
- Build scripts **auto-discover datasets by filesystem scan** under `data/<dataset>/<analysis>/`. When a new dataset's per-dataset companion (e.g., `crispr_pipeline/wg4_tf_gene_edges_FDR05.tsv`) lands, re-running the cross-dataset script picks it up — no code edits needed.
- Each script is **self-contained** — `DATASETS` dict + any loader helpers it needs are inlined at the top, no shared library import. Run with `uv run python working_groups/wg<N>_*/examples/examples.py` from the jamboree root.
- Roll-up TSVs land inside `examples/` alongside the scripts that produced them.

## Things to ask before doing

- Don't fill `[FILL IN]` / TBD sections without Adam confirming the content.
- Don't reorganize folders, rename datasets, or restructure WG layouts without asking.
- Don't add new datasets to roll-ups or schemas — only the 5 production ones.
- Mirror scripts in `data/scripts/` write to Synapse; confirm before running anything that uploads.

## Pointers

### In this folder

| Need | Go to |
|---|---|
| Participant landing page, env setup, dataset status | [`README.md`](README.md) |
| Per-dataset cards | [`data/README.md`](data/README.md) |
| Output schemas (CRISPR / cNMF / ED / metadata / guides) | [`data/schemas/`](data/schemas/) + [`data/schemas/README.md`](data/schemas/README.md) |
| Mirror / upload scripts (HPC / GCS → Synapse) | [`data/scripts/`](data/scripts/) + [`data/scripts/README.md`](data/scripts/README.md) |
| What each output column means (interpretation) | [`guide/CRISPR.md`](guide/CRISPR.md), [`guide/CNMF.md`](guide/CNMF.md), [`guide/ENERGY_DISTANCE.md`](guide/ENERGY_DISTANCE.md) |
| File-format primer, FAQ, glossary | [`guide/DATA_FORMATS.md`](guide/DATA_FORMATS.md), [`guide/FAQ.md`](guide/FAQ.md), [`guide/GLOSSARY.md`](guide/GLOSSARY.md) |
| What's in each dataset's Synapse folder | [`guide/COMPLETE_DATASET_CONTENTS.md`](guide/COMPLETE_DATASET_CONTENTS.md) |
| Cross-dataset reference tables (TF / experimental metadata, gene annotations, disease assoc., element + TF universe) | [`reference/`](reference/) (+ `_simplified` variants alongside) + [`reference/README.md`](reference/README.md) |
| Slide-deck-sized cross-dataset summaries | [`reference/cross_dataset_pipeline_summary.tsv`](reference/cross_dataset_pipeline_summary.tsv), [`reference/cross_dataset_edistance_summary.tsv`](reference/cross_dataset_edistance_summary.tsv) |
| Reference-table generators (re-run to refresh) | [`reference/scripts/`](reference/scripts/) |
| Working-group catalog + extension recipe | [`working_groups/README.md`](working_groups/README.md) |

### Parent repo (up one / two levels)

| Need | Go to |
|---|---|
| Project-wide Claude rules + daily workflow | [`../../../CLAUDE.md`](../../../CLAUDE.md) |
| Markdown style for all docs | [`../../STYLE_GUIDE.md`](../../STYLE_GUIDE.md) |
| External links (Synapse, IGVF portal, Slack, Google Docs) | [`../../LINKS.md`](../../LINKS.md) |
| Per-dataset repo layout + driver scripts | [`../../data/DATA.md`](../../data/DATA.md) |
| DACC file-format spec + open submission issues | [`../../data/DACC.md`](../../data/DACC.md) |
| 5-stage pipeline overview | [`../../analysis/ANALYSIS.md`](../../analysis/ANALYSIS.md) |
| CRISPR pipeline run guide + outputs | [`../../analysis/crispr_pipeline/`](../../analysis/crispr_pipeline/) |
| Energy-distance pipeline + outputs | [`../../analysis/energy_dist/`](../../analysis/energy_dist/) |
| cNMF + PerturbNMF run guides + outputs | [`../../analysis/cnmf/`](../../analysis/cnmf/) |
| Manuscript working folders (benchmark, lineage atlas) | [`../../manuscripts/`](../../manuscripts/) |
| IGVF portal snapshots (whole-repo, not jamboree-local) | [`../../data/portal_snapshots/`](../../data/portal_snapshots/) |
| Python package (QC, inference, portal/GCP utils) | [`../../../src/tf_perturb_seq/`](../../../src/tf_perturb_seq/) |

### Memories (auto-loaded; verify before acting on stale facts)

- Synapse-as-we-go strategy + naming for this jamboree: [[2026_utsw_jamboree]]
- HTv2 testbed scope rule: [[htv2_benchmark_role]]
- Synapse project root + auth + idempotent upload helper: [[synapse_project]]
- HPC source paths: [[hpc_project_root]], [[hpc_pipeline_mirror]]

