# tf_perturb_seq

IGVF consortium project: CRISPRi Perturb-seq of ~2000 TFs across multiple human cell lineages.

## Workflow

**GitHub Project**: [TFP3](https://github.com/users/adamklie/projects/4) — milestones as draft items, issues as tasks

**Roles**: Claude leads analysis execution. Adam reviews findings, troubleshoots, synthesizes, and provides direction.

**Daily session**:
1. Read `docs/TODAY.md`. Decide with Adam what to work on.
2. Plan the approach, or pick up from yesterday.
3. Post the plan on the relevant GitHub issue thread.
4. Execute. Minimize friction — don't ask permission for routine operations (file reads/writes, scripts, git, gsutil, uv run).
5. Update the issue thread with results, findings, blockers.
6. Iterate 2-5 until Adam signs off.

Grow documentation organically as you work. Always reference the corresponding GitHub issue.


## Quick orientation

Top-level docs (in [`docs/`](docs/)) to get up to speed:

| Doc | What it covers |
|-----|---------------|
| [TODAY.md](docs/TODAY.md) | Daily task list — what's being worked on now |
| [ROADMAP.md](docs/ROADMAP.md) | Milestones for the UTSW jamboree (May 13-16 2026) and paper Figure 1 |
| [TEAM.md](docs/TEAM.md) | Collaborators, sub-aims, dataset assignments |
| [REFERENCES.md](docs/REFERENCES.md) | External links (Synapse, IGVF portal, Google Docs, Slack, tools) |

Topical subdirs:

| Path | What it covers |
|-----|---------------|
| [docs/data/DATA.md](docs/data/DATA.md) | Per-dataset layout (`setup/` + `<run>/<tier>/`), driver scripts, full dataset inventory + status |
| [docs/data/DACC.md](docs/data/DACC.md) | DACC file-format spec audit + outstanding portal-submission issues |
| [docs/data/dataset_template/](docs/data/dataset_template/) | Scaffolding template for new datasets |
| [docs/analysis/ANALYSIS.md](docs/analysis/ANALYSIS.md) | 5-stage pipeline overview (portal → CRISPR → QC → energy distance → gene programs) |
| [docs/analysis/crispr_pipeline/](docs/analysis/crispr_pipeline/) | CRISPR FG pipeline run guide + output schema |
| [docs/analysis/energy_dist/](docs/analysis/energy_dist/) | Energy-distance pipeline + outputs |
| [docs/analysis/cnmf/](docs/analysis/cnmf/) | cNMF + PerturbNMF run guides + output schema |
| [docs/manuscripts/CRISPRi_tech_benchmark/](docs/manuscripts/CRISPRi_tech_benchmark/) | Benchmark manuscript working folder (Figure 1) |
| [docs/manuscripts/CRISPRi_lineage_atlas/](docs/manuscripts/CRISPRi_lineage_atlas/) | Lineage-atlas manuscript outline |
| [docs/jamborees/2026_UTSW/](docs/jamborees/2026_UTSW/) | UTSW jamboree packaging (Synapse `syn64423137/2026_UTSW/`) — topics, working groups, schemas |
| [docs/meetings/](docs/meetings/) | Dated meeting agendas + notes |

## Environment notes

- **git**: Run `module load git` before using git commands. The bare `git` binary is not on PATH.
- **Python 3.10**, managed with `uv` (see `pyproject.toml`) — `uv add`, not `pip install` (see memory [[feedback_use_uv_not_pip]])
- **Package**: `src/tf_perturb_seq/` (QC modules, inference, utilities)
- **Compute**: SLURM-based HPC for local jobs, GCP for CRISPR pipeline
- **Key data format**: MuData (`.h5mu`) — multimodal single-cell data
- **Analysis ecosystem**: scanpy, anndata, mudata, pandas, numpy

## Key directories

- [`datasets/`](datasets/) — One directory per experiment. Layout: `setup/{scripts,configs,samplesheets}/` + one or more `<run_name>/<tier>/` (where tier = `crispr_pipeline`, `calibration`, `qc`, `cnmf`, `energy_distance`). See [docs/data/DATA.md](docs/data/DATA.md).
- [`scripts/`](scripts/) — Shared (frozen) pipeline runners — `qc_array.sh`, `run_calibration.sh`, `run_energy_distance_pipeline.sh`, etc. Per-dataset patches live in `datasets/<ds>/bin/`.
- [`src/tf_perturb_seq/`](src/tf_perturb_seq/) — Python package (QC, inference, portal/GCP utilities)
- [`ref/`](ref/) — Reference files: `guide_libraries/`, `genome/`, `gene_sets/`, `motifs/`, `opentargets/`, `scE2G_links/`, `lab_files/`
- [`external/`](external/) — Git submodules (energy_dist_pipeline, PerturbNMF/cNMF tools)
- [`scratch/`](scratch/) — Exploratory work (gitignored), incl. `2026_03_25_legacy_archives/` for archived planning docs
