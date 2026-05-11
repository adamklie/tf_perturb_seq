# cNMF Output Directory Structure

Reference output paths:

- **Reference run (verified)**: Hon WTC11 benchmark — `aklie@nrnb-login.ucsd.edu:/cellar/users/aklie/projects/tf_perturb_seq/datasets/Hon_WTC11-benchmark_TF-Perturb-seq/PerturbNMF/`. Run name: `030726_20iter_5KHVG_torch_halsvar_batch_e7`. Total run size: ~72 GB.
- **Pipeline source**: `external/cNMF_benchmarking/cNMF_benchmarking_pipeline/` (Inference / Evaluation / Plotting subfolders).
- **Production-dataset runs**: not yet generated as of 2026-05-09 (selected-k decisions pending — see status table below).

> **Verification status — VERIFIED for the file inventory at the (k, density-threshold) level against the Hon benchmark run.** Names + per-(k,dt) sweep structure of `gene_spectra_score`, `gene_spectra_tpm`, `spectra.consensus`, `starcat_spectra`, `usages.consensus`, `clustering.<k>.<dt>.png`, `k_selection.png`, `k_selection_stats.df.npz`, `overdispersed_genes.txt`, and the 12-13 evaluation TXT filenames inside each `Eval/<k>_<dt>/` subdir are confirmed. **SAMPLED** (column-level not yet enumerated): `Annotation/<k>_<dt>.xlsx`, `Interpretation/Summary_table/<k>_<dt>/`, the integrated MuData (`adata/cNMF_<k>_<dt>.h5mu`) internal structure, and the `Plot/{Program_*,Perturb_gene_*,k_selection_*}/` PDF inventories.

The output is organized to track the four pipeline stages from `docs/analysis/cNMF.md`. Each stage writes files into `Result/<run_name>/` (and into specific subdirectories within it). The sections below walk through the stages in order.

## Top-level layout

```
PerturbNMF/
├── Data/                              # Inputs to cNMF
│   ├── inference_mudata_cleaned.h5mu     # Source MuData (gene + guide + hashing) post CRISPR pipeline
│   ├── inference_mudata_cleaned.h5ad     # Same data converted for cNMF (X = unnormalized counts)
│   ├── inference_mudata_structure.txt    # Modality/key dump of the .h5mu (auto-generated)
│   └── ensembl_to_symbol.csv             # Gene-id ↔ symbol map used in plotting
├── Result/<run_name>/                 # All cNMF outputs (see stage breakdown below)
└── Script/                            # Per-run shell + .py + .ipynb (Convert, k-selection, evaluation, plotting)
```

## Stage 1: Inference (torch-cNMF across the k sweep)

Runner: `external/cNMF_benchmarking/cNMF_benchmarking_pipeline/Inference/torch-cNMF/Slurm_Version/` (per-dataset entrypoint `datasets/<dataset_id>/6_run_cnmf.sh`). Container: `docker://igvf/torch-cnmf:v01`.

What runs: cNMF is run across a sweep of k values (benchmark: 30 values from 5..200; production: 32 values from 5..500) at multiple density thresholds (0.01, 0.05, 2.0). Each (k, dt) combination produces program loadings, cell usages, and consensus-clustering diagnostics.

**Per-(k, dt) flat outputs in `Result/<run_name>/`:**

| File | Contents | Useful for |
|------|----------|------------|
| `<run_name>.gene_spectra_score.k_<X>.dt_<Y>.txt` | Z-scored gene loadings (programs × genes). Rows = programs (1..k), cols = Ensembl IDs. | **Primary loadings**. Cross-lineage program-similarity (WG2). Top-loaded genes per program. |
| `<run_name>.gene_spectra_tpm.k_<X>.dt_<Y>.txt` | TPM-normalized gene loadings (programs × genes). Same shape as `gene_spectra_score`. | Ranking genes by absolute (rather than z-scored) contribution within a program. |
| `<run_name>.spectra.k_<X>.dt_<Y>.consensus.txt` | Consensus W matrix (programs × genes) — direct cNMF output, no normalization. | Reproducibility / advanced re-analysis that wants raw W. |
| `<run_name>.starcat_spectra.k_<X>.dt_<Y>.txt` | STARCAT-normalized spectra (alternative normalization of W). | Optional alternate loading scale. |
| `<run_name>.usages.k_<X>.dt_<Y>.consensus.txt` | Cell × program usages (H matrix). Rows = cells, cols = programs. **Bulky** at high k (157 MB at k=200). | **Primary cell-level activations**. WG1 (which programs are affected per perturbation). WG2 (program activation across cell states). |
| `<run_name>.clustering.k_<X>.dt_<Y>.png` | Per-(k,dt) clustergram + local-density histogram (cNMF's built-in QC). | Visual inspection of program separation and consensus stability. |

**Per-(k, dt) outputs in `adata/`, `loading/`, `prog_data/`:**

| File | Contents | Useful for |
|------|----------|------------|
| `adata/cNMF_<k>_<dt>.h5mu` | Integrated MuData: gene mod + guide mod + cNMF mod (cells × programs, with `varm['loadings']` = programs × genes). ~1.7 GB each. | **Single canonical handoff for downstream**. WG1/WG2 read this directly to get loadings + usages + guide assignment together. |
| `loading/cNMF_loadings_<k>_<dt>.txt` | Same loadings as `gene_spectra_score` in a slightly different format. | Duplicate of the spectra files; not mirrored. |
| `prog_data/NMF_<k>_<dt>.h5ad` | Program-level AnnData (program × gene). | Duplicate of info already in `adata/cNMF_<k>_<dt>.h5mu`; not mirrored. |

## Stage 2: Evaluation (per-(k, dt) statistical evaluation)

Runner: `external/cNMF_benchmarking/cNMF_benchmarking_pipeline/Evaluation/cNMF_evaluation_pipeline.ipynb` / `Script/cNMF_evaluation_pipeline.sh`.

What runs: For each (k, dt) combination, the evaluation step computes explained variance, GO/gene-set/trait enrichment per program, perturbation-association tests per program × per batch, and categorical-covariate associations.

**Per-(k, dt) outputs land in `Eval/<k>_<dt>/`** (one subdirectory per (k, dt); 68 in the Hon benchmark = 34 k-values × 2 density thresholds):

| File (inside `Eval/<k>_<dt>/`) | Contents | Useful for |
|------|----------|------------|
| `<k>_Explained_Variance.txt` | Per-program explained variance. | Program-level fit quality. |
| `<k>_Explained_Variance_Summary.txt` | Cumulative explained variance summary at this k. | k-selection (one of the criteria). |
| `<k>_GO_term_enrichment.txt` | GO term enrichment per program. | WG2: program biology annotation. |
| `<k>_geneset_enrichment.txt` | MSigDB / curated gene-set enrichment per program. | WG2: program biology annotation. |
| `<k>_trait_enrichment.txt` | GWAS trait enrichment per program. | WG3 (Disease / GWAS). |
| `<k>_categorical_association_results.txt` | Program-vs-categorical-covariate association tests (e.g., batch, sample). | QC: programs that are batch-confounded. |
| `<k>_categorical_association_posthoc.txt` | Post-hoc pairwise tests for the categorical associations. | QC follow-up. |
| `<k>_perturbation_association_results_<batch>.txt` | Per-program perturbation-association test results — **one file per batch** (e.g., `IGVFDS6244NAXC`, `IGVFDS8721BKRO`, ..., plus an aggregated `WTC` file in the Hon benchmark). | **Primary regulators-per-program input** (WG2). |
| `<k>_fake_perturbation_association_calibration.txt` | Optional. Negative-control perturbation calibration. Present at a subset of k values in the Hon benchmark (k = 10, 35, 40, 50, 80, 90). | Calibration / null check. |

> Why every (k, dt) bundle is kept: these tables drive the k-selection figures (Stage 3). Mirroring the full sweep here is what makes the k decision auditable — a reader can re-derive enrichment-by-k or perturbation-recovery-by-k plots at any point.

## Stage 3: k-selection

Runner: `external/cNMF_benchmarking/cNMF_benchmarking_pipeline/Evaluation/cNMF_k_selection.ipynb` / `Script/cNMF_k_selection*.sh`.

What runs: A reduced set of 10-15 candidate k values is chosen from the stability-error plot, and the k-selection notebook produces summary plots over those k values. The final k is chosen by group consensus ("clinical review board" style, per `cNMF.md`).

**Run-level outputs:**

| File | Contents | Useful for |
|------|----------|------------|
| `<run_name>.k_selection.png` | Stability-vs-error curve across the k sweep. | The canonical visual used to pick the top 10-15 candidate k values. |
| `<run_name>.k_selection_stats.df.npz` | Raw stats array behind `k_selection.png` (one row per k). | Re-plot or apply alternate k-selection criteria without re-running cNMF. |
| `<run_name>.overdispersed_genes.txt` | HVG list cNMF was run on (top-N by overdispersion). | Required to interpret loadings and to re-run with the same gene set. |
| `Plot/k_selection_<run_id>/` | Folder of k-selection PDFs from `cNMF_k_selection.ipynb`: stability-error, GO/genesets/trait enrichment by k, perturbation sensitivity by k, explained variance by k, program dot plot by conditions. | **The k-decision evidence pack.** |

**Group decision is recorded in `README.txt`** (one per run), per `cNMF.md`:

> A README.txt file with (1) the density threshold you used, (2) the k values you selected (below), and (3) any other notes.

`README.txt` is the human-readable record of *why* this k was picked. It's what lets a future reader decide if the pick still holds without re-running anything.

## Stage 4: Program plotting + Excel summarization (at the selected k)

Runner: program-QC plotting + perturbed-gene plotting + Excel-summarization notebooks (`Script/cNMF_program_analysis_<sel>_<dt>.sh`, `Script/run_perturbed_gene_analysis_<sel>_<dt>.sh`, `Script/cNMF_compile_excel_table.ipynb`).

What runs: At the selected k, generate program-level QC figures, per-perturbed-gene figures, and an integrated Annotation Excel workbook + Summary Table.

**Selected-k outputs:**

| Path | Contents | Useful for |
|------|----------|------------|
| `Annotation/<sel>_<dt>.xlsx` | Annotated workbook integrating program loadings + top genes + enrichment summaries. | Human-readable per-program summary. (Sheet/column inventory not yet documented — see TODOs.) |
| `Interpretation/Summary_table/<sel>_<dt>/` | Compiled summary tables (TSV / XLSX / TXT) integrating mdata + evaluation results at the selected k. | One-stop summary for the selected k. |
| `Plot/Program_<sel>_<dt>/` | Program-QC PDFs: program UMAP, program violin, loading correlations, top GO term plot, top loading genes, regulated-program volcano + dot + waterfall + bar plots. | WG2: program biology figures. |
| `Plot/Perturb_gene_<sel>_<dt>/` | Perturbed-gene PDFs: gene UMAP, guide UMAP, gene dotplot, gene loading correlations, top loading programs, regulated-program volcano + dot + waterfall + bar plots. | WG1: perturbation effect figures. |

The Hon benchmark example used `<sel>_<dt> = 50_2_0` so the folders are literally `Plot/Program_50_2_0/`, `Plot/Perturb_gene_50_2_0/`, etc. Production runs should rename the selected-k folder to match the chosen k.

## Subdirectory inventory (full reference)

| Subdirectory | What's in it | Stage | Mirrored? |
|---|---|---|---|
| `adata/` | `cNMF_<k>_<dt>.h5mu` per (k, dt) — integrated MuData | 1 | Selected k only |
| `loading/` | `cNMF_loadings_<k>_<dt>.txt` per (k, dt) — duplicate of spectra | 1 | No (duplicate) |
| `prog_data/` | `NMF_<k>_<dt>.h5ad` per (k, dt) — duplicate of MuData program info | 1 | No (duplicate) |
| `Eval/` | `<k>_<dt>/` per (k, dt) — 12-13 evaluation TXTs | 2 | **All (k, dt)** |
| `Evaluation/` | Alternative/legacy naming, mostly empty | 2 | No (legacy) |
| `Plot/k_selection_<run_id>/` | k-selection figure folder | 3 | Yes |
| `Plot/Program_<sel>_<dt>/` | Selected-k program-QC PDFs | 4 | Selected k only |
| `Plot/Perturb_gene_<sel>_<dt>/` | Selected-k perturbed-gene PDFs | 4 | Selected k only |
| `Annotation/` | `<k>_<dt>.xlsx` per (k, dt) — annotated workbooks | 4 | Selected k only |
| `Interpretation/Summary_table/<k>_<dt>/` | Per-(k, dt) summary tables | 4 | Selected k only |
| `cnmf_tmp/` | NPZ cache, regeneratable | 1 (intermediate) | No |
| `Inference/` | Empty / duplicate of `cnmf_tmp` | 1 (legacy) | No |
| `logs/` | SLURM `.err` / `.out` / `resource_monitor` | All | Yes (small, useful for audit) |
| `config_*.yml` | SLURM job configs (~5 per run) | All | Yes (reproducibility) |

## MuData shape after cNMF (Stage-1 output schema)

The integrated MuData written by cNMF (one per (k, dt)):

```
mdata
├── rna/                                    # Mirror of the input gene mod
│   ├── X (cells × genes, raw counts)
│   ├── obs / cell_ID, sample, ...
│   ├── var / gene_names (Ensembl IDs)
│   ├── uns / guide_names, guide_targets
│   └── obsm / guide_assignment, X_PCA, X_umap
└── cNMF/
    ├── X (cells × programs)                # Cell-level usages (H matrix, k columns)
    ├── obs (direct copy from rna.obs)
    ├── obsm (direct copy from rna.obsm; guide_assignment, X_PCA, X_umap)
    ├── uns (direct copy from rna.uns; guide_names, guide_targets)
    ├── var_names (direct copy from rna.var_names)
    └── varm / loadings (programs × genes)  # Gene-level loadings (W matrix)
```

The `<run_name>.spectra.k_<X>.dt_<Y>.consensus.txt` flat file is the same data as `cNMF.varm['loadings']`; the `<run_name>.usages.k_<X>.dt_<Y>.consensus.txt` flat file is the same as `cNMF.X`. The `gene_spectra_score` and `gene_spectra_tpm` files are post-hoc normalizations of `cNMF.varm['loadings']`.

### Required input for cNMF (Stage-1 input schema)

The `Convert_file_adata.py` script converts `inference_mudata.h5mu` → AnnData with:

```
adata
├── X (cells × genes, unnormalized counts; cNMF normalizes internally — TPM)
├── obs / cell_ID, sample (sample is required even if there's only one condition)
├── var / gene_names (Ensembl ID or symbol — must match var_names exactly)
├── uns / guide_names, guide_targets
└── obsm / guide_assignment, X_PCA, X_umap
```

`X_PCA` and `X_umap` are required when running with `shuffle_cells=True` so cells stay aligned with their guide assignments and embeddings.

## Bundle curation rule (for the 2026 UTSW jamboree)

The full per-(k, dt) sweep is ~72 GB per dataset (5 datasets × 72 GB exceeds what we can reasonably mirror to Synapse). The mirrored bundle has two purposes:

1. **Selected-k full data for downstream analysis** (WG1 + WG2): the integrated MuData, all loading variants (score / tpm / consensus / starcat), cell usages, full Eval/<sel>/ TXT bundle, selected-k Plot + Annotation + Interpretation folders.
2. **Sweep-as-provenance** so the k decision is auditable and revisitable without re-running cNMF: `k_selection.png` + raw stats, all-k clustering pngs, all-k `gene_spectra_score` loadings (small files, lets a reader spot-check programs at alternate k), all-k `Eval/<k>_<dt>/` TXT bundles, the `Plot/k_selection_<run_id>/` figure folder, and a `README.txt` with the selection rationale.

Estimated mirrored size: **~5-7 GB per dataset**. Full schema lives in `docs/jamborees/2026_UTSW/schemas/cnmf.json`.

## Per-dataset run status (2026-05-09)

| Dataset | Selected k | Status | Synapse |
|---|---|---|---|
| Hon WTC11 Cardiomyocyte | TBD | ⏳ cNMF not run yet (gated on full `crispr_pipeline/` bundle from Hon team) | — |
| Huangfu HUES8 Definitive Endoderm | TBD | ⏳ cNMF not run yet | — |
| Huangfu HUES8 Embryonic Stem Cell | TBD | ⏳ cNMF not run yet | — |
| Gersbach WTC11 Hepatocyte | TBD | ⏳ cNMF not run yet (awaiting canonical run from Gersbach team) | — |
| Engreitz WTC11 Endothelial | — | ☐ Blocked — no inference MuData (not on portal yet) | — |
| Hon WTC11 benchmark (reference, not a production dataset) | 50 | ☑ Local at `/cellar/.../PerturbNMF/Result/030726_20iter_5KHVG_torch_halsvar_batch_e7/` (used to verify schema) | — |

## TODOs

- [ ] Hold the k-selection group review for each production dataset and record the selected k + rationale in each run's `README.txt`.
- [ ] Document the `Annotation/<sel>_<dt>.xlsx` sheet/column structure (sample one workbook from the Hon benchmark).
- [ ] Document the `Interpretation/Summary_table/<sel>_<dt>/` file inventory.
- [ ] Build `docs/jamborees/2026_UTSW/scripts/mirror_cnmf_outputs.py` (HPC → Synapse) once the first production run lands.
- [ ] Build a cross-dataset program-similarity TSV (cosine similarity of `gene_spectra_score` loadings across all datasets at each dataset's selected k).
- [ ] Build a per-dataset top-20-genes-per-program TSV + regulators-per-program TSV (from `Eval/<sel>_<dt>/<sel>_perturbation_association_results_*.txt`).
