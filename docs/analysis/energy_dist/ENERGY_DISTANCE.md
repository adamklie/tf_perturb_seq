# Running the Energy Distance Pipeline

This guide walks through running the TF Perturb-seq energy-distance pipeline on an inference MuData produced by the IGVF CRISPR pipeline.

## Overview

The pipeline takes the per-dataset `inference_mudata.h5mu` (output of the CRISPR pipeline) and:

1. Filters outlier gRNAs (DISCO test + hypergeometric ranks; K-means on non-targeting guides).
2. Computes per-target **energy distance** vs. the non-targeting baseline, with a permutation test for significance (default: 20 random non-targeting backgrounds × 1000 permutations each).
3. (Step 3, separate sbatch) Builds a pairwise energy-distance matrix between all significant targets and a 2D t-SNE embedding for visualization.

Outputs land in a single `OUTPUT_FOLDER` (see [ENERGY_DISTANCE_OUTPUTS.md](ENERGY_DISTANCE_OUTPUTS.md)).

## Two ways to run

There are two implementations:

| Implementation | Repo | When to use |
|---|---|---|
| **Containerized wrapper (current)** | [`Chikara-Takeuchi/energy_dist_TFperturb`](https://github.com/Chikara-Takeuchi/energy_dist_TFperturb) | **Preferred.** Pulls a MuData from Synapse by syn ID and runs the pipeline inside an apptainer container. |
| **Direct (legacy)** | `tf_perturb_seq/external/energy_dist_pipeline/` | Reference implementation; useful for understanding internals or running on a local h5ad. |

Both produce the same output filenames; the containerized wrapper just adds a Synapse fetch + preprocess step on top.

## Prerequisites

### 1. Software requirements

| Software | Version | Notes |
|---|---|---|
| apptainer / singularity | recent | for the container; HPC modules typically have this |
| Python | 3.10+ | for the wrapper's `preprocess_mudata.py` (or use the container) |
| `git` | any | to clone both pipeline + wrapper |

### 2. Synapse personal access token

Get a PAT from [synapse.org](https://www.synapse.org) (Settings → Personal Access Tokens). The wrapper takes the token as a CLI arg, so you can either pass `$SYNAPSE_AUTH_TOKEN` directly (already in `~/.bashrc` on the HPC) or paste a fresh one.

### 3. GPU access (recommended)

Step 1 (gRNA filtering) and Step 2 (energy distance) use PyTorch. The container is GPU-aware (`apptainer exec --nv`). On the UCSD HPC, request a GPU partition (`#SBATCH -p GPUv100s` or similar).

## Quick start (this repo's setup)

There's a unified runner at `scripts/run_energy_distance_pipeline.sh` plus per-dataset SLURM submission scripts at `datasets/<dataset_id>/5_run_energy_distance.sh`. Submit one per dataset:

```bash
# From the carter-gpu submitter (where sbatch is available):
cd /cellar/users/aklie/projects/tf_perturb_seq

sbatch datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/5_run_energy_distance.sh
sbatch datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/5_run_energy_distance.sh
sbatch datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/5_run_energy_distance.sh
```

What each per-dataset script does:

1. Requests an A30 GPU + 200 GB RAM on `carter-gpu` for up to 2 days.
2. Sources MuData from GCS (Huangfu DE/ES) or Synapse (Hon CM) — the runner downloads to `OUTPUT_DIR/inference_mudata.h5mu` (skipped on resubmissions).
3. Calls the shared runner `scripts/run_energy_distance_pipeline.sh`, which:
   - Pulls the apptainer container from `docker.io/takechikara/energy_distance_env:latest` to `/cellar/users/aklie/opt/containers/edist_pipeline.sif` (one-time).
   - Clones the pipeline source from `Chikara-Takeuchi/energy_dist_pipeline`.
   - Generates `config1_2.json` and `config3.json` per run.
   - Runs `preprocess_mudata_local.py` (an in-script variant of the wrapper's preprocess that takes a local MuData), then steps 1, 2, and 2.1.
   - **Step 3 is skipped by default** — run it separately after picking cutoffs in `config3.json`.

Outputs land at:
- Hon CM: `datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/results/energy_distance/2026_04_19_no_spacer/`
- Huangfu DE: `datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/results/energy_distance/muddy_penguin/`
- Huangfu ES: `datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/results/energy_distance/sceptre_v1/`

Gersbach hepatocyte was already run by Sara (Geraghty); see [`syn74381167`](https://www.synapse.org/Synapse:syn74381167) for the layout we should match. Engreitz endothelial is gated on inference MuData being produced.

## Quick start (using the upstream wrapper directly)

If you prefer to call the upstream Chikara-Takeuchi wrapper without our scripts:

```bash
git clone https://github.com/Chikara-Takeuchi/energy_dist_TFperturb.git
cd energy_dist_TFperturb
# Edit TARGET_FILE_ID in run_all_energy_distance_pipeline.sh to point at a Synapse syn ID
sbatch run_all_energy_distance_pipeline.sh "$SYNAPSE_AUTH_TOKEN"
```

What this does end-to-end:

1. Pulls `docker.io/takechikara/energy_distance_env:latest` to `./edist_pipeline.sif`.
2. `git clone https://github.com/Chikara-Takeuchi/energy_dist_pipeline.git` — the underlying scripts.
3. Runs `preprocess_mudata.py` inside the container:
   - Fetches the MuData from Synapse.
   - Normalizes, log1p, scales, runs PCA (n_comps=50).
   - Builds an annotation table with `intended_target_name|chr:start-end` promoter labels.
4. Runs `bin/1_filtereing_gRNA.py` (filter outlier gRNAs).
5. Runs `bin/2_e_distance_nontargeting.py` (energy distance + permutation test).
6. Runs `bin/2_1_Plot_figure.py` (diagnostic plots).
7. Step 3 (`3_e_distance_among_regions.py`) is launched as a **separate** sbatch after picking p-value / energy-distance cutoffs.

Outputs land in `./data/` (`OUTPUT_FOLDER` from `config1_2.json`). See [ENERGY_DISTANCE_OUTPUTS.md](ENERGY_DISTANCE_OUTPUTS.md).

## Configuration knobs you'll usually edit

The container's `config1_2.json` controls Steps 1 & 2; `config3.json` controls Step 3.

| Knob | Default | Effect |
|---|---|---|
| `gRNA_filtering.threshold_gRNA_num` | 6 | Above this many gRNAs per target, use `combi_count` subset comparisons rather than all-vs-all. |
| `gRNA_filtering.combi_count` | 4 | Subset size for energy-distance comparisons among gRNAs targeting the same region. |
| `permutation_test.num_of_bg` | 20 | Number of random non-targeting backgrounds to permute against. Higher → more robust p-values. |
| `permutation_test.permute_per_bg` | 1000 | Permutations per background. Sets the p-value resolution (e.g. 1/1000). |
| `permutation_test.non_target_pick` | 2000 | Cells per non-targeting background sample. |
| `permutation_test.use_matched_bg` | false | If clonal effects are substantial, set true to match co-transfected gRNA composition. |
| `aggregate.downsampling_maximum` | 10000 | Max cells per target when building the pairwise similarity matrix in Step 3. |

Cutoffs for Step 3 live in `config_clustering.json` — adjust after looking at `figures/e-dist_cutoff_value.pdf` from Step 2.

## Where to find the input MuData per dataset

`docs/jamborees/2026_UTSW/synapse_paths.tsv` `crispr_pipeline` column → for each production dataset, the Synapse folder containing `pipeline_outputs/inference_mudata.h5mu`. Pass that file's syn ID (or the parent folder's, if your fetch logic walks into it) as `TARGET_FILE_ID`.

## Troubleshooting notes

- **Container can't see GPUs**: confirm `apptainer exec --nv` and that the SLURM allocation actually requested GPUs.
- **`Error: empty sequence list`** during preprocess: usually a corrupted MuData; re-pull from Synapse.
- **Permutation test runs out of memory on GPU**: lower `permutation_test.batch_num_basic` or `target_cell_num_max`.
