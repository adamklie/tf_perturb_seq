# Data

All datasets live under [`datasets/`](../../datasets/) at the repo root. This doc describes the canonical per-dataset layout, the pipeline driver scripts, what's tracked vs gitignored, and the current inventory of datasets.

For DACC submission file specs and outstanding upload issues, see [DACC.md](DACC.md).

## Dataset naming convention

```
<Lab>_<CellLine>-<condition>_TF-Perturb-seq
```

The condition slot covers either the technology benchmark (`WTC11-benchmark`) or a differentiation lineage (e.g., `HUES8-definitive-endoderm-differentiation`). For Gersbach benchmarks the chemistry is appended (`_GEM-Xv3`, `_HTv2`) since the lab contributes two.

## Per-dataset layout

The canonical layout has two levels: a `setup/` directory (shared inputs that generate metadata, configs, samplesheets) and one or more `<run_name>/` directories (one per pipeline parameter set or GCS source).

```
datasets/<Lab>_<CellLine>-<cond>_TF-Perturb-seq/
├── README.md                  # Dataset-specific notes + run table
├── setup/                     # Shared input generation (lives OUTSIDE run dirs)
│   ├── scripts/               # Numbered driver scripts (see below)
│   ├── configs/               # Base Nextflow .config(s)
│   └── samplesheets/          # sample_metadata.csv + _gcp_<date>.csv + _patched.csv
├── <run_name>/                # One per Nextflow run / parameter set / GCS source
│   ├── crispr_pipeline/       # CRISPR FG Nextflow outputs (mirror of GCS outs/)
│   │   ├── pipeline_info/     # params_*.json + versions yml — TRACKED (small)
│   │   ├── pipeline_outputs/  # bulk outputs — gitignored
│   │   ├── pipeline_dashboard/# HTML/PNG QC dashboard — gitignored
│   │   └── anndata/           # h5mu/h5ad — gitignored
│   ├── calibration/           # FDR-controlled TSVs (CRISPR FG outputs) — gitignored
│   ├── qc/                    # mapping_gene / mapping_guide / intended_target / initial_qc — gitignored
│   ├── cnmf/                  # cNMF / PerturbNMF (capital-S/D/R convention)
│   │   ├── Script/            # SLURM wrappers — TRACKED
│   │   ├── Data/              # bulk inputs — gitignored
│   │   └── Result/            # bulk results — gitignored
│   └── energy_distance/       # configs TRACKED; image/, logs/, *.csv/.h5mu/.pickle gitignored
└── <other_run_name>/          # repeat for each run
```

**Key conventions**

- `setup/` is shared across all runs — the same `sample_metadata.csv` flows into every Nextflow invocation. Never duplicate inputs inside a `<run>/`.
- `crispr_pipeline/`, `calibration/`, `qc/`, `cnmf/`, `energy_distance/` are **analysis tiers** that live in parallel inside a run dir. Not every run has every tier.
- Run-provenance lives in `<run>/crispr_pipeline/pipeline_info/params_*.json`. Folder names (especially Lucas's GCS sweep names) are often misleading — always read the params JSON to know what differs between two runs. See memory [[feedback_verify_run_params_from_pipeline_info]].
- Bulk outputs are regenerable from GCS, so they're gitignored. Small provenance + scripts are tracked. See [`.gitignore`](../../.gitignore) for the exact patterns.

## Pipeline driver scripts

Numbered shell scripts under `setup/scripts/` drive the pipeline end-to-end. The canonical sequence is:

| Script | Purpose | Output |
|---|---|---|
| `1_generate_per_sample_metadata.sh` | Query IGVF portal for the analysis set → write CSV | `setup/samplesheets/sample_metadata.csv` |
| `2_upload_to_gcp.sh` | Transfer fastqs (and CSV files) from IGVF S3 → GCS | `setup/samplesheets/sample_metadata_gcp_<date>.csv` |
| `3_patch_gcp_files.sh` | Decompress `.tsv.gz` (barcode_onlist, guide_design); strip i7/i5 reads from seqspecs where needed | `setup/samplesheets/sample_metadata_gcp_<date>_patched.csv` |
| `4_run_CRISPR_pipeline.sh` | Launch the CRISPR FG Nextflow pipeline on GCP Batch | Outputs at `gs://igvf-pertub-seq-pipeline-data/<dataset>/<date>/outs/<run_name>/` |
| `5_run_energy_distance.sh` | Launch the energy-distance pipeline (production datasets only) | `<run>/energy_distance/` |

cNMF runs are kicked off separately via SLURM wrappers placed under `<run>/cnmf/Script/` — there's no standardized `6_run_cnmf.sh`. See [analysis/cnmf/cNMF.md](../analysis/cnmf/cNMF.md) and [analysis/cnmf/PerturbNMF.md](../analysis/cnmf/PerturbNMF.md).

Local QC and calibration are driven from [`scripts/`](../../scripts/) at the repo root (e.g., `qc_array.sh`, `run_calibration.sh`, `run_energy_distance_pipeline.sh`). Per memory [[feedback_dataset_local_scripts]], scripts in `scripts/` are frozen — when a dataset needs a patched variant, copy it into `datasets/<ds>/bin/` and edit there.

## IGVF portal structure

The CRISPR FG pipeline requires an **analysis set** on the portal to drive `1_generate_per_sample_metadata.sh`. An analysis set groups measurement sets (scRNA), auxiliary sets (gRNA, HTO), and a construct library set.

A complete portal setup (see [Huangfu_WTC11-benchmark](https://data.igvf.org/analysis-sets/IGVFDS8556LYRW) as reference):

- **Analysis Set** — groups everything for the dataset
  - **Measurement Sets** — one per 10x lane/pool (scRNA), each with `strand_specificity`, `onlist_files`, `onlist_method`, and R1/R2 sequence files + seqspecs
  - **Auxiliary Sets** — one per measurement set (gRNA sequencing), linked via `measurement_sets`
  - **Construct Library Set** — exactly one, with one `integrated_content_files` entry having `content_type: "guide RNA sequences"`

Portal search for all TF Perturb-seq datasets: [data.igvf.org &nearr;](https://data.igvf.org/search/?type=MeasurementSet&preferred_assay_titles=Perturb-seq&collections=TF+Perturb-seq+Project)

## Technology benchmark datasets

Cross-lab comparison using the same WTC11 iPSC line and shared ~50-gene guide library (see memory [[project_benchmark_guide_library]] for the 30 NT + 54 OR + 8 PC + 324 TF breakdown). Per-dataset details in each README.

**Canonical CRISPR-paper runs (current):** each dataset now has **two** canonical runs that differ *only* in `GUIDE_ASSIGNMENT_method` — a `basic_threshold_sceptre` run and a `basic_threshold_cleanser` run (all other params identical within a dataset). These supersede the earlier 7-run sweep and are the ones synced locally under `datasets/<bench>/basic_threshold_{sceptre,cleanser}/crispr_pipeline/`. Full provenance + parameters for all 10 runs are in [docs/manuscripts/CRISPRi_tech_benchmark/docs/pipeline_runs.tsv](../manuscripts/CRISPRi_tech_benchmark/docs/pipeline_runs.tsv). Two gotchas (always confirm via `pipeline_info/params_*.json`, per memory [[feedback_verify_run_params_from_pipeline_info]]): the cleanser Engreitz run lives in a misleadingly-named GCS folder `Engreitz_200umi_total_ccPerturb_20mito` (its params are min_genes=800/mito=15pc, same as sceptre), and the **cleanser Gersbach_HTv2 run is incomplete** (stalled after guide assignment, no `inference_mudata.h5mu`).

The earlier **7-run parameter sweep** driven by Lucas (`cleanser_500`, `cleanser_800`, `cleanser_extremes_200`, `cleanser_extremes_2000`, `cleanser_knee2`, `scrublet_off_cleanser_800`, `scrublet_on_sceptre_800`, all `mito_15pc`) is retained in each dataset dir but is superseded for paper analyses.

| Dataset | Lab | Technology | Portal Accession |
|---|---|---|---|
| [Hon_WTC11-benchmark](../../datasets/Hon_WTC11-benchmark_TF-Perturb-seq/) | Hon | 10x 5' HT v2 w/ HTO | [IGVFDS4761PYUO](https://data.igvf.org/analysis-sets/IGVFDS4761PYUO) |
| [Huangfu_WTC11-benchmark](../../datasets/Huangfu_WTC11-benchmark_TF-Perturb-seq/) | Huangfu | 10x 3' v3 | [IGVFDS8556LYRW](https://data.igvf.org/analysis-sets/IGVFDS8556LYRW) |
| [Engreitz_WTC11-benchmark](../../datasets/Engreitz_WTC11-benchmark_TF-Perturb-seq/) | Engreitz | CC Perturb-seq | [IGVFDS5057HJKP](https://data.igvf.org/analysis-sets/IGVFDS5057HJKP) |
| [Gersbach_WTC11-benchmark_GEM-Xv3](../../datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/) | Gersbach | 10x GEM-X 3' | [IGVFDS6673ZFFG](https://data.igvf.org/analysis-sets/IGVFDS6673ZFFG) |
| [Gersbach_WTC11-benchmark_HTv2](../../datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/) | Gersbach | 10x HT v2 | [IGVFDS6237URFJ](https://data.igvf.org/analysis-sets/IGVFDS6237URFJ) |

For cross-tech analysis (Figure 1 of the paper), milestone tracking lives in [ROADMAP.md](../ROADMAP.md) (#15).

## Production datasets

Full TF library (~2000 targets) applied to differentiated lineages.

| Dataset | Lab | Lineage | Local runs on disk | Status |
|---|---|---|---|---|
| [Hon_WTC11-cardiomyocyte-differentiation](../../datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/) | Hon | Cardiomyocyte | `seqspec_v3/` (prod, GCS 2026_04_15), `weizhou_syn74520421/` (QC-only), `2026_04_19_no_spacer/` (ED-only) | CRISPR + cNMF Stage 1 done; tracked in issue [#20](https://github.com/adamklie/tf_perturb_seq/issues/20) |
| [Huangfu_HUES8-definitive-endoderm](../../datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/) | Huangfu | Definitive endoderm (DE) | `muddy_penguin/` (GCS 2026_04_09) | CRISPR + cNMF + ED + QC complete; portal analysis set still TBD (see below) |
| [Huangfu_HUES8-embryonic-stemcell](../../datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/) | Huangfu | Embryonic stem cell (ESC) | `sceptre_v1/` (GCS 2026_04_13) | CRISPR + cNMF + ED + QC complete; portal analysis set still TBD |
| [Gersbach_WTC11-hepatocyte-differentiation](../../datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/) | Gersbach | Hepatocyte | `sara_synapse_syn74842722/` (Sara's h5mu, QC-only locally) | CRISPR ran on Sara's side; local QC re-run done. See issue [#28](https://github.com/adamklie/tf_perturb_seq/issues/28) |

The Huangfu HUES8 production runs were processed without going through `1_generate_per_sample_metadata.sh` (the portal analysis set is still missing) — they were driven from raw fastqs supplied directly by the lab.

### Huangfu HUES8 portal status (audited 2026-03-25)

Both HUES8 datasets have measurement sets, auxiliary sets, and a construct library set on the portal; the **only missing piece for the standard workflow** is the analysis set that groups them. Local production runs (`muddy_penguin/`, `sceptre_v1/`) are processed and downstream analyses (cNMF, ED) are complete, but DACC submission still depends on resolving the portal gaps below.

|  | Measurement Sets | Auxiliary Sets (gRNA) | Construct Library Set | Analysis Set | Seqspecs |
|---|:---:|:---:|:---:|:---:|:---:|
| HUES8 Stemcell | 8 | 8 | [IGVFDS3299AXST](https://data.igvf.org/construct-library-sets/IGVFDS3299AXST) (guide: [IGVFFI8270UPKB](https://data.igvf.org/tabular-files/IGVFFI8270UPKB)) | **None** | **Missing** on files |
| HUES8 Endoderm | 8 | 8 | same as above | **None** | **Missing** on files |
| WTC11 Benchmark (for reference) | 4 | 4 | separate CLS | [IGVFDS8556LYRW](https://data.igvf.org/analysis-sets/IGVFDS8556LYRW) | Present |

**Measurement-set metadata** (`strand_specificity`, `onlist_files`, `onlist_method`) looks good across both datasets.

**HUES8 Stemcell measurement → auxiliary pairs:**
| Measurement Set | Auxiliary Set (gRNA) |
|---|---|
| IGVFDS0746MYRH | IGVFDS8263ODBC |
| IGVFDS0956UUQN | IGVFDS2137GIIX |
| IGVFDS2908DIHX | IGVFDS2178ZGIR |
| IGVFDS2940LYGK | IGVFDS3927XOZY |
| IGVFDS3348KRMQ | IGVFDS0173JYTZ |
| IGVFDS5934EZKT | IGVFDS9550TKIZ |
| IGVFDS6247FKDR | IGVFDS6962YCLR |
| IGVFDS8623ONYE | IGVFDS4608EUUL |

**HUES8 Endoderm measurement → auxiliary pairs:**
| Measurement Set | Auxiliary Set (gRNA) |
|---|---|
| IGVFDS0788TUIL | IGVFDS7673CCYF |
| IGVFDS1313VGHL | IGVFDS5105RAVI |
| IGVFDS1403KZYV | IGVFDS0367BQIW |
| IGVFDS1437XFJK | IGVFDS6182MVVF |
| IGVFDS2129VHBD | IGVFDS2661AFQU |
| IGVFDS2520ZEYH | IGVFDS7765TYSM |
| IGVFDS3434HGEI | IGVFDS5564SHUB |
| IGVFDS4079VXPX | IGVFDS1454WEUZ |

**Action items** (coordinate with Huangfu lab / Denis Torre / DACC):
1. **Create analysis sets** — one for stemcell (8 MS + 8 aux), one for endoderm (8 MS + 8 aux). Needed for portal-driven submission, not for local processing.
2. **Add seqspecs** to R1 sequence files (currently none linked). Either upload to the portal or provide fallback YAML paths at runtime (the benchmark used the fallback route).
3. Tag the auxiliary sets with `collections: TF Perturb-seq Project` for discoverability.

Legacy internal processing outputs from before the standardized pipeline have been archived to `scratch/2026_03_25_legacy_archives/`.

## Public datasets

| Dataset | Description |
|---|---|
| [ENCODE_WTC11_ChIP-seq](../../datasets/ENCODE_WTC11_ChIP-seq/) | 90 ENCODE TF ChIP-seq experiments in WTC11; used for binding-site validation of perturbation targets. Flat layout (no `setup/` or `<run>/` — just manifests + `downloads/`). |

## Shared data locations

- **Synapse**: [syn63675917](https://www.synapse.org/Synapse:syn63675917) (project root; per-dataset mirrors linked from each README)
- **GCS**: `igvf-pertub-seq-pipeline` project, `gs://igvf-pertub-seq-pipeline-data/` bucket
- **IGVF Portal**: [data.igvf.org](https://data.igvf.org/)
- **Guide metadata on portal**: [IGVFFI5765HMZH](https://data.igvf.org/tabular-files/IGVFFI5765HMZH/)
- **Lucas's parameter-sweep root** (benchmark datasets): `gs://igvf-pertub-seq-pipeline-data/scratch/bioinfolucas/<lucas_dir>/<dataset_subdir>/`

## Reference files (`ref/`)

Located in [`ref/`](../../ref/) at the repo root:

- `guide_libraries/` — Guide-to-target mapping files used by pipelines (incl. `finalized_annotation_files/`)
- `genome/` — GTF genome annotation files (GENCODE V43)
- `gene_sets/` — Control guide lists, gene sets, target lists
- `motifs/` — TF motif PWMs
- `opentargets/`, `scE2G_links/`, `lab_files/` — Auxiliary references

## Adding a new dataset

1. Create `datasets/<Lab>_<CellLine>-<condition>_TF-Perturb-seq/` with `setup/{scripts,configs,samplesheets}/` and a `README.md`.
2. Copy and adapt driver scripts from [`dataset_template/`](dataset_template/) into `setup/scripts/` — fill in the IGVF accession, dataset name, and Nextflow paths.
3. Run scripts 1 → 4 (and 5 for production datasets). Each script writes its output back to `setup/samplesheets/`.
4. After Nextflow finishes on GCP, mirror outputs to `<run_name>/crispr_pipeline/` — commit `pipeline_info/` only; the bulk subdirs are gitignored.
5. Add the dataset row to the table above, update the README with a "Pipeline runs" table.
