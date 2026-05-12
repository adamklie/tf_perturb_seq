# cNMF / PerturbNMF — output directory reference

What's in a per-dataset PerturbNMF run directory, and which downstream consumer cares about each file. Pair with [`PerturbNMF.md`](PerturbNMF.md), which is the how-to-run companion.

**Reference runs** (both completed end-to-end through Stage 3e at K=200, 2026-05-09):

- DE: `datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/muddy_penguin/cnmf/042926_huangfu_de_torchcnmf_KskillA/`
- ESC: `datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/sceptre_v1/cnmf/042926_huangfu_esc_torchcnmf_KskillA/`

Every run uses the **same K sweep** (8 values: 30, 50, 60, 80, 100, 200, 250, 300) and a **single density threshold** `dt = 2.0` (the canonical downstream value — multiple-dt sweeps from the legacy `cNMF_benchmarking` runs are no longer produced). Selected-K plots and tables use the suffix `<K>_2_0`.

## Top-level layout

```
datasets/<dataset>/<run>/cnmf/<run_name>/
├── README.md                              # per-run summary (selected K, status table, K rationale)
├── Data/                                  # inputs to cNMF
│   ├── <dataset>_<run>_perturbnmf.h5ad    # inference input AnnData (X = raw counts; built by Convert_file_adata.py)
│   └── guide_annotation.tsv               # guide_id renamed to guide_names; required by U-test calibration
├── Script/                                # per-stage SLURM submission scripts + helper .py
│   ├── <run_name>_inference.sh            # Stage 1
│   ├── prepare_h5mu_for_eval.{py,sh}      # adds obs['sample']='all' + remaps NT targets
│   ├── inject_umap_into_h5mu.{py,sh}      # pre-computes UMAP into h5mu (remediation only — see PerturbNMF.md)
│   ├── cNMF_evaluation_pipeline.sh        # Stage 2a
│   ├── cNMF_evaluation_trait_only.sh      # Stage 2a follow-up after OpenTargets is available
│   ├── U-test_perturbation_calibration.sh # Stage 2b (needs fix/utest-oom-leak branch)
│   ├── cNMF_k_selection.sh                # Stage 3a
│   ├── cNMF_program_analysis_<K>_2_0.sh   # Stage 3b — currently OOMs (upstream #7)
│   ├── cNMF_perturbed_gene_analysis_<K>_2_0.sh   # Stage 3c
│   ├── cNMF_compile_excel_summary.{py,sh} # Stage 3e
│   └── upload_to_synapse.py
└── Result/<run_name>/                     # all outputs (see stage-by-stage breakdown below)
```

## Stage 1: Inference (torch-cNMF across the K sweep)

Runner: [`external/PerturbNMF/src/Stage1_Inference/`](../../../external/PerturbNMF/src/Stage1_Inference/) via per-run `Script/<run_name>_inference.sh`. Container: `docker://igvf/torch-cnmf:v01`.

What runs: torch-cNMF is run across the 8-K sweep at `dt = 2.0`. Each (K, dt) yields program loadings, cell usages, and consensus-clustering diagnostics. Run-level diagnostics (stability/error across the full K sweep) are written once.

**Per-K outputs in `Result/<run_name>/Inference/`** (the file basename is the run-level `Inference.` prefix):

| File | Contents | Useful for |
|------|----------|------------|
| `Inference.gene_spectra_score.k_<K>.dt_2_0.txt` | Z-scored gene loadings (programs × genes). | **Primary loadings.** Cross-lineage program similarity (WG2). Top-loaded genes per program. |
| `Inference.gene_spectra_tpm.k_<K>.dt_2_0.txt` | TPM-normalized gene loadings (programs × genes). | Ranking genes by absolute (not z-scored) contribution. |
| `Inference.spectra.k_<K>.dt_2_0.consensus.txt` | Consensus W matrix (programs × genes), raw cNMF output. | Reproducibility / advanced re-analysis. |
| `Inference.starcat_spectra.k_<K>.dt_2_0.txt` | STARCAT-normalized spectra (alternate normalization of W). | Optional alternate loading scale. |
| `Inference.usages.k_<K>.dt_2_0.consensus.txt` | Cell × program usages (H matrix). Rows = cells, cols = programs. **Bulky** at high K (~150 MB at K=200). | **Primary cell-level activations.** WG1 (programs affected per perturbation). WG2 (activations across cell states). |
| `Inference.clustering.k_<K>.dt_2_0.png` | Per-(K, dt) clustergram + local-density histogram. | Visual inspection of program separation / consensus stability. |

**Run-level outputs in `Result/<run_name>/Inference/`:**

| File | Contents | Useful for |
|------|----------|------------|
| `Inference.k_selection.png` | Stability-vs-error curve across the full K sweep. | The canonical visual used to pick a defensible selected K. |
| `Inference.k_selection_stats.df.npz` | Raw stats array behind `k_selection.png` (one row per K). | Re-plot or re-apply alternate selection criteria. |
| `Inference.overdispersed_genes.txt` | HVG list cNMF was run on. | Required to interpret loadings or re-run on the same gene set. |

**Subdirectories in `Result/<run_name>/Inference/`:**

| Subdir | Contents | Mirrored to Synapse? |
|--------|----------|----------------------|
| `adata/` | `cNMF_<K>_2_0.h5mu` (×8 K values) — **integrated MuData**: rna mod + cNMF mod (cells × programs, with `varm['loadings']` = programs × genes). ~1.7 GB each. | Selected K only (~5–6 GB) |
| `cnmf_tmp/` | NPZ cache from factorize/refit steps. Regeneratable. | No |
| `loading/`, `prog_data/` | Duplicates of `gene_spectra_score` / per-K program data. | No (duplicates) |
| `diagnosis_plots/` | Per-K clustering pngs + scatter plots. | All K (small) |
| `Annotation/` | Auto-generated per-K xlsx workbooks. | No (superseded by Stage 3e Summary_table) |

The `<run>/adata` symlink at the `Result/<run_name>/` level points to `Inference/adata/` — workaround for upstream PerturbNMF [issue #6](https://github.com/EngreitzLab/PerturbNMF/issues/6) (path mismatch between `<run>/adata` and `Inference/adata`).

## Stage 2a: Evaluation (per-K statistical tests + enrichments)

Runner: [`external/PerturbNMF/src/Stage2_Evaluation/`](../../../external/PerturbNMF/src/Stage2_Evaluation/) via `Script/cNMF_evaluation_pipeline.sh`. Outputs land in `Result/<run_name>/Evaluation/<K>_2_0/`, one subdir per K.

| File (per `<K>_2_0/` subdir) | Contents | Useful for |
|------|----------|------------|
| `<K>_perturbation_association_results_all.txt` | Per-program perturbation-association test results — target × program × log2FC × q-value. **Single file** because the project convention sets `obs['sample']='all'` (pooled), so there's no per-batch stratification. | **Primary regulators-per-program input** (WG2). |
| `<K>_geneset_enrichment.txt` | MSigDB / curated gene-set enrichment per program. | WG2 program biology annotation. |
| `<K>_GO_term_enrichment.txt` | GO term enrichment per program. | WG2 program biology annotation. |
| `<K>_trait_enrichment.txt` | GWAS trait enrichment per program. Requires OpenTargets resource at `external/PerturbNMF/src/Stage2_Evaluation/Resources/OpenTargets_L2G_Filtered.csv.gz` — added in a follow-up `cNMF_evaluation_trait_only.sh` pass. | WG3 disease / GWAS. |
| `<K>_Explained_Variance.txt` | Per-program explained variance at this K. | Program-level fit quality. |
| `<K>_Explained_Variance_Summary.txt` | Cumulative explained variance summary at this K. | One of the K-selection criteria. |

**Evaluations we deliberately skip** (project-specific):

- **Categorical association** — skipped because `obs['sample']` is single-valued (`'all'`).
- **Motif enrichment** — requires `hg38.fa` + cell-type-specific enhancer-gene links (e.g., scE2G); the latter is unavailable for HUES8 endoderm/ESC, so we skip.

## Stage 2b: U-test fake-target calibration

Runner: [`external/PerturbNMF/src/Stage2_Evaluation/`](../../../external/PerturbNMF/src/Stage2_Evaluation/) via `Script/U-test_perturbation_calibration.sh`. **Requires** the local `fix/utest-oom-leak` branch (see [PerturbNMF.md](PerturbNMF.md#required-fixes-local-branches-not-yet-merged-upstream)).

| File (per `<K>_2_0/` subdir) | Contents | Useful for |
|------|----------|------------|
| `<K>_fake_perturbation_association_results.txt` | 50 fake-targeting iterations using random subsets of NT controls; per-program null distribution. | Calibrating the per-K FDR threshold — is the perturbation-hit set noise-distinguishable? |

## Stage 3a: K-selection panel

Runner: `Script/cNMF_k_selection.sh`. Outputs land in `Result/<run_name>/Plot/k_selection/`.

| File | Contents | Useful for |
|------|----------|------------|
| `K-selection_panel_2.0.{png,svg}` | The canonical 6-panel summary used to pick selected K. | **Primary K-decision artifact.** |
| `Stability_Error_{stability,error}.{png,svg}` | The two underlying stability + error curves. | Inputs to the panel. |
| `Explained_Variance_2.0.{png,svg}` | Explained variance vs K curve. | One panel of the K-selection summary. |
| `Enrichment_2.0_{genesets,go_terms,traits}.{png,svg}` | Per-metric enrichment-count vs K curve. | Enrichment panels. |
| `Perturbation_2.0_{all_samples,per_sample}.{png,svg}` | Number of significant perturbations recovered vs K. | Perturbation-recovery panel. |
| `Program_dotplot_<K>_2.0.png` | Program × condition dotplot, one per K. | Visual check of program-condition specificity per K. |

K-selection is a group decision ("clinical review board" style historically). The rationale is captured in the run's `README.md`.

## Stage 3b: Per-program PDFs — **currently deferred**

Runner: `Script/cNMF_program_analysis_<K>_2_0.sh`. Target output dir: `Result/<run_name>/Plot/Program_<K>_2_0/`.

**Status: skipped on all production runs.** OOMs in the pre-loop correlation precompute at full data and is glacially slow even on a 10% subsample (~30 min/program × K = days). Tracked upstream as [issue #7](https://github.com/EngreitzLab/PerturbNMF/issues/7). The per-program view is covered by the Stage 3e Summary sheet; the per-TF view is covered by Stage 3c.

The `Plot/Program_<K>_2_0/` folder may exist but is empty. A `Program_<K>_2_0_thinned/` variant has been used as a 10%-subsample workaround on the DE run.

## Stage 3c: Per-target (perturbed-gene) PDFs

Runner: `Script/cNMF_perturbed_gene_analysis_<K>_2_0.sh`. Outputs land in `Result/<run_name>/Plot/Perturb_gene_<K>_2_0/`.

| File | Contents | Useful for |
|------|----------|------------|
| `<TF>.pdf` (one per perturbed target) | Volcano + gene UMAP + guide UMAP + per-program log2FC + correlation waterfall for the knockdown. | WG1 — the canonical "what happens when this TF is knocked down" figure. |
| `merged_perturbed_genes_<dataset>_K<K>.pdf` | All per-TF PDFs glued together (~580 MB for DE / ~888 MB for ESC at K=200). Built with `pdfunite` because upstream's PyPDF2 merge hangs on thousands of PDFs ([issue #8](https://github.com/EngreitzLab/PerturbNMF/issues/8)). | One-file skim across the whole library. |

Occasional per-target PDFs fail to render in the parallel matplotlib pass (e.g., 7 empty PDFs on DE) — they're re-run single-threaded via `Script/cNMF_perturbed_gene_analysis_<K>_2_0_redo<N>.sh`.

## Stage 3e: Excel summary

Runner: `Script/cNMF_compile_excel_summary.sh`. Output lands in `Result/<run_name>/Interpretation/Summary_table/<K>_2_0/`.

| File | Sheets | Useful for |
|------|--------|------------|
| `cNMF_<K>_2_0.xlsx` | **Summary** — one row per program: top loaded genes, top regulators, top GO / geneset / trait terms, explained variance. **Targets Summary** — one row per TF: which programs it most strongly regulates, expression baseline, cell count. **Perturbation Association** — full target × program × log2FC × q-value table (split across sheets when >1M rows). | **The primary human-readable per-program / per-target view of the whole run.** Review starts here. |

## Subdirectory inventory (full reference)

| Subdirectory under `Result/<run_name>/` | Contents | Stage |
|---|---|---|
| `Inference/adata/` | `cNMF_<K>_2_0.h5mu` per K — integrated MuData | 1 |
| `Inference/` (flat files) | Per-K spectra / usages / clustering + run-level k_selection / overdispersed_genes | 1 |
| `Inference/{loading,prog_data}/` | Duplicate of spectra / per-K program data | 1 (intermediate) |
| `Inference/cnmf_tmp/` | NPZ factorize/refit cache | 1 (intermediate) |
| `Inference/diagnosis_plots/` | Per-K consensus diagnostics | 1 |
| `Inference/Annotation/` | Auto-generated per-K annotated xlsx (superseded by Stage 3e) | 1 |
| `adata -> Inference/adata` | Symlink — workaround for upstream issue [#6](https://github.com/EngreitzLab/PerturbNMF/issues/6) | 1 |
| `Evaluation/<K>_2_0/` | Per-K eval TXTs (perturbation, geneset, GO, trait, EV, fake-test) | 2a + 2b |
| `Plot/k_selection/` | K-selection panel + per-metric K curves + per-K program dotplots | 3a |
| `Plot/Program_<K>_2_0/` | Selected-K per-program PDFs — currently deferred | 3b |
| `Plot/Perturb_gene_<K>_2_0/` | Selected-K per-TF PDFs + merged combined PDF | 3c |
| `Interpretation/Summary_table/<K>_2_0/` | Excel summary workbook | 3e |
| `logs/`, `config_*.yml` | SLURM job logs + per-stage Hydra YAML configs | All (small, kept for audit) |

## MuData shape after Stage 1 (the integrated `cNMF_<K>_2_0.h5mu`)

```
mdata
├── rna/                                    # mirror of the inference-input gene modality
│   ├── X                                   # cells × genes, raw counts
│   ├── obs                                 # cell_ID, sample (= 'all' after h5mu prep), ...
│   ├── var                                 # gene_names (Ensembl IDs)
│   ├── uns                                 # guide_names, guide_targets
│   └── obsm                                # guide_assignment, X_PCA, X_umap
└── cNMF/
    ├── X                                   # cells × programs — cell usages (H matrix, K columns)
    ├── obs / obsm / uns                    # copies from rna (incl. guide_assignment + UMAP)
    ├── var_names                           # copy from rna.var_names
    └── varm['loadings']                    # programs × genes — gene loadings (W matrix)
```

The flat `Inference.spectra.k_<K>.dt_2_0.consensus.txt` and `Inference.usages.k_<K>.dt_2_0.consensus.txt` files are the same data as `cNMF.varm['loadings']` and `cNMF.X` respectively. The `gene_spectra_score` and `gene_spectra_tpm` files are post-hoc normalizations of the loadings.

### Stage-1 input schema (what `Convert_file_adata.py` produces)

The script converts the project's `inference_mudata.h5mu` → an AnnData with:

```
adata
├── X                                       # cells × genes, raw counts (cNMF TPM-normalizes internally)
├── obs                                     # cell_ID, sample (required even if single-valued)
├── var                                     # gene_names — must match var_names exactly
├── uns                                     # guide_names, guide_targets
└── obsm                                    # guide_assignment, X_PCA, X_umap
```

`X_PCA` and `X_umap` are required when running with `shuffle_cells=True` so cells stay aligned with their guide assignments and embeddings — pass `--compute_umap` to populate them at this step (see [PerturbNMF.md](PerturbNMF.md#project-specific-h5mu-prep-convention)).

## Per-dataset run status (2026-05-12)

| Dataset | Selected K | Status | Synapse |
|---|---|---|---|
| Huangfu HUES8 Definitive Endoderm | 200 (250 may be preferred) | ✅ Stage 1, 2a, 2b, 3a, 3c, 3e | [`syn74893844`](https://www.synapse.org/Synapse:syn74893844) |
| Huangfu HUES8 Embryonic Stem Cell | 200 | ✅ Stage 1, 2a, 2b, 3a, 3c, 3e | [`syn74893846`](https://www.synapse.org/Synapse:syn74893846) |
| Hon WTC11 Cardiomyocyte | TBD | ⏳ Setup only — Stage 1 not yet kicked off | — |
| Gersbach WTC11 Benchmark HTv2 (testbed) | TBD | Stage 1 run; Stage 2/3 pending. Verifies the pipeline before launching remaining production runs. | — |
| Gersbach WTC11 Hepatocyte | — | ☐ Blocked — awaiting canonical run from Gersbach team | — |
| Engreitz WTC11 Endothelial | — | ☐ Blocked — no inference MuData (not on portal yet) | — |

Stage 3b is skipped on all runs (upstream issue [#7](https://github.com/EngreitzLab/PerturbNMF/issues/7)). The Hon WTC11 and Huangfu WTC11 *benchmark* datasets each have an older `030726_20iter_5KHVG_torch_halsvar_batch_e7` run on disk; those used the predecessor `cNMF_benchmarking` tool with a 3-density-threshold sweep (`dt = 0.01, 0.05, 2.0`) — superseded for production work, retained for historical comparison.
