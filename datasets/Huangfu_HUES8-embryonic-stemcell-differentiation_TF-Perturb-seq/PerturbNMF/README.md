# PerturbNMF — Huangfu HUES8 Embryonic Stem Cell

End-to-end PerturbNMF run on the Huangfu HUES8 stem cell differentiation dataset (~189,595 cells × ~36,000 genes, 14,151 guides → 2,063 unique targets including 600 NT controls). Run started 2026-04-29 (inference) and completed 2026-05-09 (Stage 3).

## Run identifiers

| Field | Value |
|---|---|
| Run name | `042926_huangfu_esc_torchcnmf_KskillA` |
| Input h5ad | `Data/ESC_sceptre_v1_perturbnmf.h5ad` |
| K values | 30, 50, 60, 80, 100, 200, 250, 300 |
| sel_thresh | 2.0 |
| Inference algo | torch-cNMF, halsvar, batch mode, tol=1e-4, 10 iter, 2000 HVG |
| Categorical key (perturbation tests) | `sample` (single value `'all'` — pooled across batches; not stratified) |
| Default K for plotting | **K = 200** (revisable; see K-selection panel) |

## How to review this run

Top-down skim, ~15 minutes:

1. Open `Plot/k_selection/K-selection_panel_2.0.png` — confirm K=200 is reasonable.
2. Open the Summary sheet of `Interpretation/Summary_table/200_2_0/cNMF_200_2_0.xlsx` — one row per program with top genes, top regulators, top enrichment terms, explained variance. This is your "what programs are there + which TFs control them" overview.
3. For any TF of interest, open `Plot/Perturb_gene_200_2_0/<TF>.pdf` — volcano + UMAP + per-program log2FC + correlation waterfall. The merged PDF (`merged_perturbed_genes_ESC_K200.pdf`, 888 MB) is all of these glued together.
4. For raw analysis, `Evaluation/200_2_0/200_perturbation_association_results_all.txt` is the full target × program × log2FC × q-value table.

See [docs/analysis/PerturbNMF.md](../../../docs/analysis/PerturbNMF.md#what-to-look-at--review-checklist) for the broader rationale.

## What's been run

| Stage | Status | Outputs |
|---|---|---|
| 1. Inference (torch-cNMF) | ✅ | `Inference/adata/cNMF_<K>_2_0.h5mu` (× 8 K values) |
| h5mu prep (sample col + NT remap + UMAP inject) | ✅ | applied to K=200 h5mu |
| 2a. Evaluation — perturbation, geneset, GO, explained variance | ✅ | `Evaluation/<K>_2_0/<K>_*.txt` |
| 2a. Evaluation — trait enrichment | ✅ | `Evaluation/<K>_2_0/<K>_trait_enrichment.txt` (added in follow-up after OpenTargets file became available) |
| 2a. Evaluation — categorical association | ⏭️ Skipped by design (single-value `sample`) |
| 2a. Evaluation — motif enrichment | ⏭️ Skipped — no scE2G enhancer-gene links for HUES8 stem cell |
| 2b. U-test fake calibration | ✅ | `Evaluation/<K>_2_0/<K>_fake_perturbation_association_results.txt` (× 8 K, 50 iter each) |
| 3a. K-selection panel | ✅ | `Plot/k_selection/` — `K-selection_panel_2.0.png/svg`, plus per-metric panels |
| 3b. Per-program PDFs | ❌ Cancelled | OOMs in pre-loop correlation precompute even at 700 GB; tracked upstream as [issue #7](https://github.com/EngreitzLab/PerturbNMF/issues/7) |
| 3c. Per-target PDFs | ✅ | `Plot/Perturb_gene_200_2_0/<TF>.pdf` (1948 PDFs; 0 empty in this dataset) |
| 3c. Merged combined PDF | ✅ | `Plot/Perturb_gene_200_2_0/merged_perturbed_genes_ESC_K200.pdf` (~888 MB) — built with `pdfunite` since upstream's PyPDF2 merge hangs |
| 3e. Excel summary | ✅ | `Interpretation/Summary_table/200_2_0/cNMF_200_2_0.xlsx` |
| 3d. LLM annotation | ⏭️ Skipped per Adam |

## Directory layout

```
PerturbNMF/
├── README.md                                # this file
├── Data/
│   ├── ESC_sceptre_v1_perturbnmf.h5ad     # inference input (raw counts)
│   └── guide_annotation.tsv                 # built from ref/finalized_annotation_files/harmonized_guide_file_poolabcd.tsv with `guide_id` renamed to `guide_names`; required by U-test fake-test path
├── Result/042926_huangfu_esc_torchcnmf_KskillA/
│   ├── Inference/                           # Stage 1 output
│   │   ├── adata/cNMF_<K>_2_0.h5mu          # ×8 K values
│   │   ├── cnmf_tmp/                        # intermediate factorize/refit artifacts
│   │   ├── Inference.k_selection_stats.df.npz
│   │   ├── diagnosis_plots/                 # per-K clustering pngs + scatter
│   │   └── Annotation/                      # auto-generated per-K xlsx
│   ├── adata -> Inference/adata             # symlink (workaround for upstream issue #6)
│   ├── Evaluation/<K>_2_0/                  # Stage 2a per-K outputs
│   │   ├── <K>_perturbation_association_results_all.txt
│   │   ├── <K>_geneset_enrichment.txt
│   │   ├── <K>_GO_term_enrichment.txt
│   │   ├── <K>_trait_enrichment.txt
│   │   ├── <K>_Explained_Variance.txt
│   │   ├── <K>_Explained_Variance_Summary.txt
│   │   └── <K>_fake_perturbation_association_results.txt   # Stage 2b (after OOM fix branch applied)
│   ├── Plot/
│   │   ├── k_selection/                     # Stage 3a panels
│   │   ├── Program_200_2_0/                 # Stage 3b — empty (cancelled, see issue #7)
│   │   └── Perturb_gene_200_2_0/            # Stage 3c per-TF PDFs + merged combined PDF
│   └── Interpretation/Summary_table/200_2_0/
│       └── cNMF_200_2_0.xlsx                # Stage 3e
└── Script/                                  # SLURM submission scripts (this dataset)
    ├── 042926_huangfu_esc_torchcnmf_KskillA_inference.sh
    ├── prepare_h5mu_for_eval.sh             # adds obs['sample']='all' + remaps NT targets
    ├── inject_umap_into_h5mu.sh             # gene-based UMAP for K=200 h5mu
    ├── cNMF_evaluation_pipeline.sh
    ├── cNMF_evaluation_trait_only.sh        # follow-up trait-only pass
    ├── U-test_perturbation_calibration.sh
    ├── cNMF_k_selection.sh
    ├── cNMF_program_analysis_200_2_0.sh     # Stage 3b — currently OOMs
    ├── cNMF_perturbed_gene_analysis_200_2_0.sh
    ├── cNMF_perturbed_gene_analysis_200_2_0_redo7.sh   # single-thread redo for genes that failed in parallel pass
    └── cNMF_compile_excel_summary.sh
```

## How this run was put together

The PerturbNMF source on `external/PerturbNMF/` is checked out on local branch `fix/utest-oom-leak`, which is 2 commits ahead of `origin/main`:

1. `9665129` — fixes `args.reference_targets` AttributeError in U-test fake-test path ([PR #4](https://github.com/EngreitzLab/PerturbNMF/pull/4))
2. `520cb15` — adds `del + gc.collect()` per fake-test iter / per K to keep the deep-copy memory bounded ([PR #5](https://github.com/EngreitzLab/PerturbNMF/pull/5))

External resources downloaded from [EngreitzLab/gene_network_evaluation](https://github.com/EngreitzLab/gene_network_evaluation/tree/main/smk/resources):
- `OpenTargets_L2G_Filtered.csv.gz` → `external/PerturbNMF/src/Stage2_Evaluation/Resources/`
- `hocomoco_meme.meme` → same (motif enrichment not yet run)

Upstream issues filed:
- [#2](https://github.com/EngreitzLab/PerturbNMF/issues/2), [#3](https://github.com/EngreitzLab/PerturbNMF/issues/3) — already closed by PRs #4/#5
- [#6](https://github.com/EngreitzLab/PerturbNMF/issues/6) — `<run>/adata` vs `Inference/adata` path mismatch (workaround: symlink)
- [#7](https://github.com/EngreitzLab/PerturbNMF/issues/7) — Stage 3b precompute OOM
- [#8](https://github.com/EngreitzLab/PerturbNMF/issues/8) — `merge_pdfs_in_folder` hang on thousands of PDFs (workaround: `pdfunite`)
- [#9](https://github.com/EngreitzLab/PerturbNMF/issues/9) — Stage 3c HVG breaks on raw counts (workaround: pre-inject UMAP)

## h5mu prep convention (project-wide)

For any TFP3 dataset post-Stage-1, before Stage 2/3 run cleanly:
1. `prepare_h5mu_for_eval.py` — add `obs['sample']='all'` + remap NT `guide_targets` from `'nan'` → `'non-targeting'`.
2. `inject_umap_into_h5mu.py` — gene-based UMAP on `mdata['rna']` (normalize_total → log1p → HVG → scale → PCA → neighbors → UMAP) written into both `mdata['rna'].obsm` and `mdata['cNMF'].obsm`.
3. `Data/guide_annotation.tsv` — built from `ref/finalized_annotation_files/harmonized_guide_file_poolabcd.tsv` with `guide_id` renamed to `guide_names`.

Going forward, UMAP should be computed once on the inference INPUT h5ad so the embedding propagates through cNMF naturally.

## K choice and refinement plan

K = 200 was chosen as a defensible mid-range value across DE and ESC for first-pass plotting. The K-selection panel (`Plot/k_selection/K-selection_panel_2.0.png`) suggests:
- Stability is best at K = 30 (silhouette 0.42) and recovers at K ≥ 200; trough at K = 80–100
- Explained variance peaks at K = 250, then drops at K = 300 (overfitting signal)
- Enrichment metrics (GO, geneset, trait) increase monotonically with K

**For DE specifically: K = 250 may be slightly preferred** based on explained-variance peak. Worth re-running 3c at K = 250 for a publication pass. Stage 3b would need either an upstream fix or a workaround.

## Synapse / IGVF portal upload status

- **Synapse**: ✅ uploaded to [`syn74893846`](https://www.synapse.org/Synapse:syn74893846) (`2026_UTSW/datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/cnmf/`) on 2026-05-10. Layout follows [`docs/jamborees/2026_UTSW/schemas/cnmf.json`](../../../docs/jamborees/2026_UTSW/schemas/cnmf.json) (no run_name nesting; tighter curation per `mirror_cnmf_outputs.py`). Bundle includes `adata/`, `Evaluation/`, `Plot/k_selection/`, `Plot/Perturb_gene_200_2_0/merged_*.pdf`, `Interpretation/Summary_table/200_2_0/cNMF_200_2_0.xlsx`, run-level `Inference.*` files, `Config/` with labeled per-stage SLURM YAML, and `README.txt` with K rationale. **Excluded** (intentional): per-target individual PDFs, `gene_spectra_tpm`, `starcat_spectra`, `Annotation/<sel>_2.0.xlsx`, per-K MuDatas other than K=200.
- **IGVF DACC**: not yet uploaded. Catalog spec requires (per `docs/data/DACC.md`) gene universe (HVG list with Ensembl IDs), gene programs, gene program regulators (perturbation tables) — files need column-name remapping before submission.
