# IGVF TF Perturb-seq Pillar Project (TFP3)

Code and analysis for a multi-lab IGVF consortium effort: CRISPRi Perturb-seq of ~2000 transcription factors across multiple human cell lineages.

## Quick start

```bash
# Clone with submodules
git clone --recurse-submodules https://github.com/adamklie/tf_perturb_seq.git
cd tf_perturb_seq

# Or initialize submodules in an existing clone
git submodule update --init --recursive

# Install dependencies (Python 3.10, managed via uv)
uv sync
```

## Repository layout

```
tf_perturb_seq/
├── datasets/                       # One directory per experiment
│   └── <Lab>_<CellLine>-<cond>_TF-Perturb-seq/
│       ├── README.md
│       ├── setup/                  # Shared inputs (scripts, configs, samplesheets)
│       │   ├── scripts/            # 1_–5_*.sh pipeline drivers
│       │   ├── configs/            # Base Nextflow .config(s)
│       │   └── samplesheets/       # sample_metadata.csv → _gcp_<date>.csv → _patched.csv
│       └── <run_name>/             # One per parameter set / GCS source
│           ├── crispr_pipeline/    # Nextflow outputs (pipeline_info/ tracked; bulk gitignored)
│           ├── calibration/        # FDR-controlled TSVs (gitignored)
│           ├── qc/                 # mapping_gene / mapping_guide / intended_target (gitignored)
│           ├── cnmf/                # cNMF / PerturbNMF (Script/ tracked; Data/, Result/ gitignored)
│           └── energy_distance/    # configs tracked; image/, logs/, results gitignored
├── src/tf_perturb_seq/             # Python package: qc/, inference/, cnmf/, energy_dist/, crispr_pipeline/, portal/, gcp/
├── scripts/                        # Shared (frozen) pipeline runners
│   ├── qc_array.sh                 # SLURM array QC runner
│   ├── run_qc_pipeline.sh
│   ├── run_calibration.sh
│   ├── run_energy_distance_pipeline.sh
│   ├── run_pathways.sh
│   ├── audit_dataset.py
│   └── sync_gcp_run.sh
├── ref/                            # Reference files: guide_libraries/, genome/, gene_sets/, motifs/, ...
├── external/                       # Git submodules
│   ├── energy_dist_pipeline/
│   ├── PerturbNMF/
│   └── cNMF-TF-perturbseq-CMs-Honlab/
├── docs/                           # See "Documentation" below
└── scratch/                        # Exploratory work (gitignored)
```

See [docs/data/DATA.md](docs/data/DATA.md) for the per-dataset layout conventions (tracked vs gitignored, run-provenance, etc.).

## Documentation

| Doc | What it covers |
|-----|----------------|
| [docs/TODAY.md](docs/TODAY.md) | Daily task list |
| [docs/ROADMAP.md](docs/ROADMAP.md) | Milestones (UTSW jamboree May 13–16 2026 + paper Figure 1) |
| [docs/TEAM.md](docs/TEAM.md) | Collaborators, sub-aims, dataset assignments |
| [docs/REFERENCES.md](docs/REFERENCES.md) | External links (Synapse, IGVF portal, Google Docs, Slack) |
| [docs/data/DATA.md](docs/data/DATA.md) | Dataset inventory, layout, driver scripts |
| [docs/data/DACC.md](docs/data/DACC.md) | DACC file-format audit + portal-submission gaps |
| [docs/analysis/ANALYSIS.md](docs/analysis/ANALYSIS.md) | 5-stage pipeline overview |
| [docs/analysis/crispr_pipeline/](docs/analysis/crispr_pipeline/) | CRISPR FG pipeline run guide + output schema |
| [docs/analysis/energy_dist/](docs/analysis/energy_dist/) | Energy-distance pipeline + outputs |
| [docs/analysis/cnmf/](docs/analysis/cnmf/) | cNMF + PerturbNMF run guides + output schema |
| [docs/manuscripts/](docs/manuscripts/) | Tech-benchmark and lineage-atlas manuscript working folders |
| [docs/jamborees/2026_UTSW/](docs/jamborees/2026_UTSW/) | Jamboree packaging (Synapse `syn64423137/2026_UTSW/`) |

## Pipeline overview

The workflow has 5 stages. Top-level overview: [docs/analysis/ANALYSIS.md](docs/analysis/ANALYSIS.md); per-stage docs under [docs/analysis/](docs/analysis/).

```
Raw fastqs → [1] IGVF Portal → [2] CRISPR Pipeline → [3] QC → [4] Energy Distance → [5] Gene Programs (cNMF)
                (upload)         (GCP/Nextflow)       (local)    (local/SLURM)         (local/SLURM)
```

Stage-to-script mapping (per dataset):

| Stage | What happens | How it's run |
|---|---|---|
| **1. IGVF Portal** | Data producers upload raw fastqs to [data.igvf.org](https://data.igvf.org/) (measurement sets + auxiliary sets + construct library set → analysis set). On our side, we query the portal for metadata and stage files on GCS. | `setup/scripts/1_generate_per_sample_metadata.sh` → `2_upload_to_gcp.sh` → `3_patch_gcp_files.sh` |
| **2. CRISPR Pipeline** | Nextflow run of the [IGVF CRISPR Pipeline](https://github.com/IGVF/CRISPR_Pipeline) on GCP Batch. Produces `inference_mudata.h5mu`. | `setup/scripts/4_run_CRISPR_pipeline.sh` → outputs mirrored locally to `<run>/crispr_pipeline/` |
| **3. QC** | Local QC on the inference MuData (mapping_gene, mapping_guide, intended_target). | Shared SLURM driver: `scripts/qc_array.sh` → outputs to `<run>/qc/` |
| **4. Energy Distance** | Permutation tests + phenotype clustering via the [energy_dist_pipeline](https://github.com/Chikara-Takeuchi/energy_dist_pipeline) submodule. | `setup/scripts/5_run_energy_distance.sh` (or shared `scripts/run_energy_distance_pipeline.sh`) → `<run>/energy_distance/` |
| **5. Gene Programs (cNMF)** | Multi-stage PerturbNMF / torch-cNMF: GPU K-sweep → evaluation (perturbation assoc + GO/geneset/trait + EV) → U-test FDR calibration → K-selection panel (group decision) → per-target PDFs → Excel summary. ~3–5 h GPU + ~12–16 h CPU per dataset. | No single driver — per-stage SLURM wrappers under `<run>/cnmf/Script/`. Needs h5mu prep ([[project_perturbnmf_h5mu_prep]]) and per-stage compute budget ([[project_perturbnmf_compute_budget]]). See [docs/analysis/cnmf/PerturbNMF.md](docs/analysis/cnmf/PerturbNMF.md). |

Quick end-to-end (CRISPR pipeline path only):

```bash
DS=datasets/<Lab>_<CellLine>-<cond>_TF-Perturb-seq

# Stage 1: portal query + GCS staging
bash $DS/setup/scripts/1_generate_per_sample_metadata.sh   # → setup/samplesheets/sample_metadata.csv
bash $DS/setup/scripts/2_upload_to_gcp.sh                  # → ..._gcp_<date>.csv
bash $DS/setup/scripts/3_patch_gcp_files.sh                # → ..._gcp_<date>_patched.csv

# Stage 2: launch Nextflow on GCP
bash $DS/setup/scripts/4_run_CRISPR_pipeline.sh

# Stage 4: energy distance (production datasets)
bash $DS/setup/scripts/5_run_energy_distance.sh
```

Stage 3 (QC) is run from shared SLURM scripts (`scripts/qc_array.sh`), not from `setup/scripts/`. Stage 5 (cNMF) is a multi-stage workflow with its own per-stage SLURM scripts and group K-selection step — see [docs/analysis/cnmf/PerturbNMF.md](docs/analysis/cnmf/PerturbNMF.md).

## Datasets

### Technology benchmark (WTC11)

Cross-lab comparison using shared WTC11 iPSC line and ~50-gene guide library.

| Dataset | Lab | Technology | Portal |
|---|---|---|---|
| [Hon_WTC11-benchmark](datasets/Hon_WTC11-benchmark_TF-Perturb-seq/) | Hon | 10x 5' HT v2 w/ HTO | [IGVFDS4761PYUO](https://data.igvf.org/analysis-sets/IGVFDS4761PYUO) |
| [Huangfu_WTC11-benchmark](datasets/Huangfu_WTC11-benchmark_TF-Perturb-seq/) | Huangfu | 10x 3' v3 | [IGVFDS8556LYRW](https://data.igvf.org/analysis-sets/IGVFDS8556LYRW) |
| [Engreitz_WTC11-benchmark](datasets/Engreitz_WTC11-benchmark_TF-Perturb-seq/) | Engreitz | CC Perturb-seq | [IGVFDS5057HJKP](https://data.igvf.org/analysis-sets/IGVFDS5057HJKP) |
| [Gersbach_WTC11-benchmark_GEM-Xv3](datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/) | Gersbach | 10x GEM-X 3' | [IGVFDS6673ZFFG](https://data.igvf.org/analysis-sets/IGVFDS6673ZFFG) |
| [Gersbach_WTC11-benchmark_HTv2](datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/) | Gersbach | 10x HT v2 | [IGVFDS6237URFJ](https://data.igvf.org/analysis-sets/IGVFDS6237URFJ) |

All five share an identical 7-run parameter sweep (`cleanser_500`, `cleanser_800`, `cleanser_extremes_{200,2000}`, `cleanser_knee2`, `scrublet_off_cleanser_800`, `scrublet_on_sceptre_800`, all `mito_15pc`). Benchmark manuscript: [docs/manuscripts/CRISPRi_tech_benchmark/](docs/manuscripts/CRISPRi_tech_benchmark/).

### Production datasets

Full TF library (~2000 targets) on differentiated lineages.

| Dataset | Lab | Lineage | Local runs |
|---|---|---|---|
| [Hon_WTC11-cardiomyocyte-differentiation](datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/) | Hon | Cardiomyocyte | `seqspec_v3/`, `weizhou_syn74520421/`, `2026_04_19_no_spacer/` |
| [Huangfu_HUES8-definitive-endoderm](datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/) | Huangfu | Definitive endoderm | `muddy_penguin/` |
| [Huangfu_HUES8-embryonic-stemcell](datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/) | Huangfu | Embryonic stem cell | `sceptre_v1/` |
| [Gersbach_WTC11-hepatocyte-differentiation](datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/) | Gersbach | Hepatocyte | `sara_synapse_syn74842722/` |

### Public datasets

| Dataset | Description |
|---|---|
| [ENCODE_WTC11_ChIP-seq](datasets/ENCODE_WTC11_ChIP-seq/) | 90 ENCODE TF ChIP-seq experiments in WTC11; binding-site validation |

Full status and external tracking sheets: [docs/data/DATA.md](docs/data/DATA.md) and [docs/REFERENCES.md](docs/REFERENCES.md).

## Adding a new dataset

1. Scaffold from [`docs/data/dataset_template/`](docs/data/dataset_template/) into `datasets/<Lab>_<CellLine>-<cond>_TF-Perturb-seq/` with `setup/{scripts,configs,samplesheets}/` and a `README.md`.
2. Adapt the driver scripts (fill in IGVF accession, dataset name).
3. Run scripts 1 → 4 (and 5 for production datasets). Each writes its output back to `setup/samplesheets/`.
4. After Nextflow finishes on GCP, mirror outputs to `<run_name>/crispr_pipeline/` — commit `pipeline_info/` only.
5. Add the dataset row to [docs/data/DATA.md](docs/data/DATA.md) and update the dataset README with a "Pipeline runs" table.
