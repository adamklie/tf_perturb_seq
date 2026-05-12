# Compute budget (empirical, TFP3 production scale)

Empirical SLURM resource budgets for ~190–270 k cell datasets at K = {30, 50, 60, 80, 100, 200, 250, 300} and sel_thresh = 2.0. Measured on UCSD nrnb in May 2026.

For canonical reference, see project memory `project_perturbnmf_compute_budget.md`.

## Per-stage table

| Stage | CPUs | Memory | GPU | Wall time (actual) | Wall time (alloc) |
|---|---|---|---|---|---|
| **Stage 1 — Inference** (torch-cNMF) | 4 | 96–128 GB | 1× A30 | 3–5 h | 10–12 h |
| **h5mu prep + UMAP inject** | 4 | 256 GB | — | 8–15 min | 1 h |
| **Stage 2a — Evaluation (all 5 of 9 metrics)** | 20 | 128 GB | — | ~9 h | 12 h |
| **Stage 2a — Trait-only follow-up** | 10 | 96 GB | — | ~3 min | 1 h |
| **Stage 2b — U-test calibration** (with fixes) | 4 | 256 GB | — | 50–80 min | 4 h |
| **Stage 3a — K-selection plot** | 2 | 96 GB | — | 1.5–4 min | 30 min |
| **Stage 3b — Program analysis** | 20 | **>700 GB OOM** | — | infeasible | — |
| **Stage 3c — Perturbed-gene PDFs** | 20 | 256 GB | — | 1–2 h | 4 h |
| **Stage 3e — Excel summary** | 4 | 128 GB | — | ~10 min | 1 h |

**Total wall time end-to-end** (excluding Stage 3b): ~15–20 h compute spread across ~2 days clock time (Stage 1 GPU queue + Stage 2/3 serialization).

## Notes on each stage

### Stage 1 (GPU, 3–5 h)

- Bottleneck: GPU memory for the torch-cNMF iterations. A30 with 24 GB fits ~270 k cells × 5K HVG comfortably; bump to A100 only if you're going past 500 k cells.
- 8 K values × default 20 iterations × ~14 s/iter ≈ 45 min — but only 1 K runs at a time on the GPU. Time scales linearly with K-count.
- RAM: holds the AnnData in CPU memory throughout; 96 GB tight for production-scale, 128 GB comfortable.
- Alloc 10–12 h to absorb queue+startup overhead; actual ~4 h.

### Stage 2a (CPU, ~9 h)

- The 5 metrics actually computed: `perturbation_association`, `geneset_enrichment`, `GO_term_enrichment`, `trait_enrichment`, `Explained_Variance`. Motif enrichment is skipped (no scE2G for HUES8).
- Trait enrichment requires `OpenTargets_L2G_Filtered.csv.gz` in `external/PerturbNMF/src/Stage2_Evaluation/Resources/`.
- Per-K work is independent; the SLURM script processes Ks sequentially. 20 CPUs help inside the perturbation association test (parallelized over guides).
- 128 GB tight at K=300 with full TF library; raise to 256 GB if K extends past 300.

### Stage 2a trait-only follow-up (CPU, ~3 min)

- Used to re-compute trait enrichment alone after dropping in a new OpenTargets file. Reuses cached perturbation outputs.

### Stage 2b U-test calibration (CPU, 50–80 min)

- **Requires the `fix/utest-oom-leak` branch** — see `02-required-fixes.md`. Without it, OOMs at iteration ~150/400.
- 50 fake-targeting iterations × 8 Ks = 400 inner loops. Each iter samples NT subsets as fake targets and re-runs the perturbation test.
- 256 GB needed even with fixes — production-scale fake-test runs are memory-hungry by design.
- 4 CPUs is fine; bottleneck is single-threaded NumPy / pandas.

### Stage 3a K-selection plot (CPU, 1.5–4 min)

- Reads all per-K eval CSVs + the h5mus (for program_dotplot panels).
- 96 GB tight only if you load all 8 h5mus at once for the dotplot. The skill defaults to lazy loading.
- Output: `K-selection_panel_2.0.png` / `.svg` + per-metric panels.

### Stage 3b — DON'T RUN

- Upstream issue [#7](https://github.com/EngreitzLab/PerturbNMF/issues/7). Pre-compute OOMs at full data (>700 GB); per-program plotting glacial even on a 10%-subsampled h5mu (~30 min/program → 5 days for K=200).
- Coverage gap is filled by Stage 3e Excel (per-program summary sheets) + Stage 3c (per-TF PDFs).

### Stage 3c perturbed-gene PDFs (CPU, 1–2 h)

- One PDF per perturbed target. K=200 ⇒ ~2000 PDFs for full TF library, ~300 for benchmark.
- Parallelism: matplotlib backend is fork-safe; 20 CPUs help if PerturbNMF passes `--n_jobs 20` (check the SLURM submit script).
- 256 GB needed: each worker loads its own copy of the h5mu.
- Final merge step uses `pdfunite` (not the upstream `merge_pdfs_in_folder` which hangs — see upstream issue #8).

### Stage 3e Excel summary (CPU, ~10 min)

- Reads all per-K eval CSVs + h5mus, writes one xlsx with sheets per K + Summary + Targets Summary.
- 128 GB tight; bump if K extends past 300.

## Sizing benchmark datasets

Benchmark datasets (~50-gene library, 2–5 measurement sets) are ~10× smaller across the board:

| Stage | Approx scaling vs production |
|---|---|
| Stage 1 | 30–60 min on A30 |
| Stage 2a | ~1 h |
| Stage 2b | 5–15 min |
| Stage 3a | <1 min |
| Stage 3c | 5–15 min |
| Stage 3e | 1–2 min |

Most stages comfortable in 64 GB. Default SLURM submission scripts in `Script/` size for production; benchmark runs just finish faster.

## When to override the budget

- **K extends past 300:** double the memory request on all stages that load h5mus.
- **>500 k cells:** consider A100 over A30 for Stage 1; bump Stage 2a memory to 256 GB.
- **K count below 4:** halve every wall time (most stages are linear in K-count).
- **More fake-test iterations** (`--n_fake_iter 100`): linear scaling on Stage 2b time.

## Reference runs

Both completed end-to-end at K=200, 2026-05-09:

- **DE:** `datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/muddy_penguin/cnmf/042926_huangfu_de_torchcnmf_KskillA/`
- **ESC:** `datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/sceptre_v1/cnmf/042926_huangfu_esc_torchcnmf_KskillA/`

Inspect their `README.md` for the per-stage status table and any per-stage anomalies.
