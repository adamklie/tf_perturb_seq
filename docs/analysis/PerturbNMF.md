# PerturbNMF — Project conventions and run notes

This doc captures the project-specific conventions and known quirks for running [PerturbNMF](https://github.com/EngreitzLab/PerturbNMF) on TFP3 datasets. It complements the upstream README and skill docs.

PerturbNMF supersedes the older `cNMF_benchmarking` tool ([docs/analysis/cNMF.md](cNMF.md)). The repo is checked into this project at [external/PerturbNMF/](../../external/PerturbNMF/).

---

## Pipeline at a glance

```
Stage 1 — Inference (GPU torch-cNMF)
    Input: counts h5ad → 8× cNMF_<K>_2_0.h5mu (one per K)
    Time: 3–5 h on A30; 96–128 GB RAM

Stage 2a — Evaluation (CPU; per K, 5 of 9 metrics for our run)
    Inputs: h5mu + Stage 2a resources (OpenTargets, etc.)
    Outputs: <K>_perturbation_association_results_<sample>.txt
             <K>_geneset_enrichment.txt
             <K>_GO_term_enrichment.txt
             <K>_trait_enrichment.txt          (needs OpenTargets file)
             <K>_Explained_Variance.txt
    Time: ~9 h for 8 K values, 128 GB RAM, 20 CPUs
    Trait-only follow-up rerun: ~3 min once perturbation outputs exist

Stage 2b — U-test calibration (CPU; 50 fake-targeting iterations × 8 K)
    Inputs: h5mu + guide_annotation.tsv + Stage 2a perturbation results
    Output: <K>_fake_perturbation_association_results.txt
    Time: ~50–80 min, 256 GB RAM, 4 CPUs
    Requires local fix branch — see "Required fixes" below

Stage 3a — K-selection plot (CPU)
    Inputs: all per-K eval CSVs + h5mu (for program_dotplot panels)
    Output: K-selection_panel_2.0.png/svg + per-metric panels
    Time: 1.5–4 min, 96 GB RAM

Stage 3b — Program analysis (per-program PDFs)  [BLOCKED for ~270k cells × K=200]
Stage 3c — Perturbed-gene analysis (per-TF PDFs)
    Inputs: h5mu + perturbation results
    Output: <TF>.pdf per perturbed target + merged combined PDF
    Time: 1–2 h on 256 GB RAM, 20 CPUs

Stage 3e — Excel summary (compile per-K results)
    Inputs: same as 3a
    Output: cNMF_<K>_<sel>.xlsx
    Time: ~10 min, 128 GB RAM
```

---

## Project-specific h5mu prep convention

After Stage 1 inference, **before** Stage 2/3 runs, three prep steps make the
upstream code work cleanly on our datasets:

1. **Add `obs['sample'] = 'all'` + remap NT `guide_targets`** — `prepare_h5mu_for_eval.py`. The single-value sample column means perturbation tests are pooled across batches (not stratified). The NT remap fixes the `'nan'` → `'non-targeting'` convention so the U-test fake-test path can find reference guides.

2. **Pre-compute UMAP from the gene matrix** — `inject_umap_into_h5mu.py`. Standard scanpy recipe (`normalize_total → log1p → HVG → scale → PCA → neighbors → UMAP`) on a copy of `mdata['rna']`, written into both `mdata['rna'].obsm` and `mdata['cNMF'].obsm`. Sidesteps Stage 3b's slow rna-PCA and Stage 3c's broken HVG call (upstream issue [#9](https://github.com/EngreitzLab/PerturbNMF/issues/9)). UMAP depends only on the rna matrix (identical across K), so for a given dataset you only need to inject the *selected K* h5mu — defer this step until Stage 3a K-selection has picked a k.

3. **Generate `Data/guide_annotation.tsv`** — copy from `ref/finalized_annotation_files/harmonized_guide_file_poolabcd.tsv` with `guide_id` renamed to `guide_names`. Required by the U-test fake-test code path.

**Going forward (read this before the next dataset's Stage 1 run)** — compute UMAP on the inference INPUT (the h5ad / inference_mudata before Stage 1) so the embedding propagates through cNMF naturally and step 2 becomes unnecessary. `Convert_file_adata.py` (the per-dataset script that builds the cNMF AnnData input) now takes a `--compute_umap` flag that runs the same `normalize_total → log1p → HVG → scale → PCA → neighbors → UMAP` recipe on a temp copy and stashes `obsm['X_pca']` + `obsm['X_umap']` in the output AnnData. Stage 1 carries them through into every per-K h5mu's rna modality, so `inject_umap_into_h5mu.py` becomes a remediation-only tool for h5mu files produced before this convention. **Pass `--compute_umap` whenever you invoke Convert_file_adata.py for a new dataset.** Reference implementation: [`datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/PerturbNMF/Script/Convert_file_adata.py`](../../datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/PerturbNMF/Script/Convert_file_adata.py).

---

## Required fixes (local branches not yet merged upstream)

We ran on `external/PerturbNMF/` checkout on local branch `fix/utest-oom-leak`, two commits ahead of `origin/main`:

| Commit | What | Upstream PR |
|---|---|---|
| `9665129` | Fix `args.reference_targets` AttributeError in U-test fake-test | [#4](https://github.com/EngreitzLab/PerturbNMF/pull/4) |
| `520cb15` | Add `del + gc.collect()` per fake-test iter and per K (memory leak fix) | [#5](https://github.com/EngreitzLab/PerturbNMF/pull/5) |

Both are also pushed to [adamklie/PerturbNMF](https://github.com/adamklie/PerturbNMF).

Open issues filed (no fix yet):
- [#6](https://github.com/EngreitzLab/PerturbNMF/issues/6) — `<run>/adata` vs `Inference/adata` path mismatch — workaround: per-dataset `<run>/adata → Inference/adata` symlink
- [#7](https://github.com/EngreitzLab/PerturbNMF/issues/7) — Stage 3b infeasible at production scale (precompute OOM at full data; per-program plotting glacially slow even on a 10%-subsampled h5mu — ~30 min/program → ~5 days for K=200). Stage 3b deliberately skipped. Per-program view is covered by the Stage 3e Excel Summary sheet; per-TF view is covered by Stage 3c merged Perturb_gene PDF.
- [#8](https://github.com/EngreitzLab/PerturbNMF/issues/8) — `merge_pdfs_in_folder` hangs on thousands of PDFs — workaround: `pdfunite`
- [#9](https://github.com/EngreitzLab/PerturbNMF/issues/9) — Stage 3c HVG breaks on raw counts — workaround: pre-inject UMAP

---

## External resources

Downloaded from [EngreitzLab/gene_network_evaluation/smk/resources](https://github.com/EngreitzLab/gene_network_evaluation/tree/main/smk/resources) and placed at `external/PerturbNMF/src/Stage2_Evaluation/Resources/`:

- `OpenTargets_L2G_Filtered.csv.gz` — required for `--Perform_trait`
- `hocomoco_meme.meme` — for motif enrichment (not yet wired in upstream)

Motif enrichment also requires `hg38.fa` and an enhancer-gene linking file (e.g., scE2G); the latter is cell-type specific and unavailable for HUES8 endoderm/ESC, so we currently skip motif enrichment.

---

## What to look at — review checklist

PerturbNMF answers four questions per dataset:

1. **What gene programs exist?** — cNMF discovers K coregulated gene programs. Each is a vector of gene loadings + per-cell scores.
2. **Which TFs regulate which programs?** — for each TF perturbation, test whether knockdown shifts program scores vs NT controls. Output: target × program × log2FC × q-value.
3. **What does each program mean biologically?** — GO/geneset/trait enrichment on top-loaded genes per program.
4. **Is the FDR real?** — U-test fake-target calibration generates a null distribution from random NT-guide subsets to validate that perturbation hits aren't noise.

When opening the outputs, work top-down:

| Priority | Artifact | What it tells you |
|---|---|---|
| 1 | `Plot/k_selection/K-selection_panel_2.0.png` | Was K appropriate? Stability, explained variance, enrichment counts vs K |
| 2 | `Interpretation/Summary_table/<K>_2_0/cNMF_<K>_2_0.xlsx` — **Summary** sheet | One row per program: top loaded genes, top regulators, top GO/geneset/trait terms, explained variance |
| 3 | Same xlsx — **Targets Summary** sheet | One row per TF: which programs it most strongly regulates, expression baseline, # cells |
| 4 | Same xlsx — **Perturbation Association** sheets | Full target × program × log2FC × q-value table (split across sheets if >1M rows) |
| 5 | `Plot/Perturb_gene_200_2_0/<TF>.pdf` for TFs of interest | Volcano + UMAP + per-program effects for that knockdown |
| 6 | `Evaluation/<K>_2_0/<K>_perturbation_association_results_all.txt` | Raw stats for custom analysis |

The merged combined PDF (`merged_perturbed_genes_<dataset>_K<K>.pdf`) glues all the per-TF PDFs into one — handy for skimming the whole library at once.

For **cross-dataset comparison** (e.g., DE vs ESC):
- Same K → same number of programs → compare program-level signatures (top-loaded genes, perturbation effects)
- Use the Summary sheet's top-loaded genes to find program correspondences across datasets
- Or compute pairwise correlation between DE and ESC program-score matrices on shared cells (none here since different cells, but on shared TFs the perturbation log2FCs are comparable)

---

## Runs in flight / completed

See per-dataset PerturbNMF READMEs:

- DE: [datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/PerturbNMF/README.md](../../datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/PerturbNMF/README.md)
- ESC: [datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/PerturbNMF/README.md](../../datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/PerturbNMF/README.md)

Both datasets have completed Stage 1, Stage 2a (5/9 metrics), Stage 2b U-test, Stage 3a K-selection, Stage 3c per-target PDFs at K=200. Stage 3b is deferred per upstream issue #7. Stage 3e Excel summary in progress.

---

## Compute budget reference

Empirical SLURM resource budgets for ~190–270k-cell datasets at K = {30,50,60,80,100,200,250,300}, sel_thresh = 2.0 — see project memory `project_perturbnmf_compute_budget.md` (or `cellar/.../memory/project_perturbnmf_compute_budget.md`).

| Stage | CPUs | Memory | Time |
|---|---|---|---|
| Stage 1 (GPU) | 4 + A30 | 96–128 GB | 3–5 h actual / 10–12 h alloc |
| h5mu prep + UMAP inject | 4 | 256 GB | 8–15 min |
| Stage 2a all metrics (no trait/motif) | 20 | 128 GB | ~9 h |
| Stage 2a trait-only follow-up | 10 | 96 GB | ~3 min |
| Stage 2b U-test (with OOM fix) | 4 | 256 GB | 50–80 min |
| Stage 3a K-selection plot | 2 | 96 GB | 1.5–4 min |
| Stage 3b program analysis | 20 | **>700 GB OOMs** | infeasible without upstream fix |
| Stage 3c perturbed-gene | 20 | 256 GB | 1–2 h (parallel matplotlib) |
| Stage 3e Excel summary | 4 | 128 GB | ~10 min |

---

## Synapse upload

Per-dataset upload script: `<dataset>/PerturbNMF/Script/upload_to_synapse.py`. Synapse project id `syn64423137`, layout `<project>/PerturbNMF/<dataset>/<run_name>/`.

```bash
SYNAPSE_AUTH_TOKEN=...  # from ~/.bashrc or rotated PAT
python <ds>/PerturbNMF/Script/upload_to_synapse.py --dataset DE --upload-only-k200-h5mu
# add --dry-run first to preview
```

Per-K h5mu's are large (~5–6 GB at K=200; total ~30+ GB if all 8 K). The default skip-all-h5mu-except-K=200 mode keeps total upload to ~7–8 GB per dataset.

---

## Related docs

- [docs/analysis/cNMF.md](cNMF.md) — older `cNMF_benchmarking` instructions (predecessor of PerturbNMF)
- [docs/data/DACC.md](../data/DACC.md) — IGVF DACC catalog file format requirements
- [docs/REFERENCES.md](../REFERENCES.md) — external links (Synapse, Google Sheets, Slack)
