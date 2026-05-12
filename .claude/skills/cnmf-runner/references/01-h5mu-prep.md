# h5mu prep (TFP3-specific)

Stage 1 of PerturbNMF (torch-cNMF inference) expects a single AnnData, not a MuData. And Stages 2/3 expect a few project-specific obs/var conventions. These prep steps bridge the gap.

## Step 1: Convert MuData → AnnData

Wraps `src/tf_perturb_seq/cnmf/h5mu_to_perturbnmf_h5ad.py`.

### Column mapping

| Source (`inference_mudata.h5mu`) | Destination (output `.h5ad`) |
|---|---|
| `mod['gene'].X` | `X` (cells × genes counts) |
| `mod['gene'].var['symbol']` (column) | `var.index` (gene SYMBOLS — required by spec) |
| `mod['gene'].var` (whole) | `var` (preserves mt/ribo flags) |
| `mod['gene'].obs` | `obs` (must include a categorical batch key, default `'batch'`) |
| `mod['guide'].layers['guide_assignment']` | `obsm['guide_assignment']` (sparse cells × guides) |
| `mod['guide'].var[<guide_id_col>]` | `uns['guide_names']` |
| `mod['guide'].var[<target_col>]` | `uns['guide_targets']` |

Authoritative spec: `external/PerturbNMF/.claude/skills/perturbNMF-runner/references/data-format-spec.md`.

### Run

```bash
source <REPO_ROOT>/.venv/bin/activate
python <REPO_ROOT>/src/tf_perturb_seq/cnmf/h5mu_to_perturbnmf_h5ad.py \
  --in_h5mu  <REPO_ROOT>/datasets/<DS>/<RUN>/crispr_pipeline/pipeline_outputs/inference_mudata.h5mu \
  --out_h5ad <REPO_ROOT>/datasets/<DS>/<RUN>/cnmf/<cnmf_run_name>/Data/<DS>_<RUN>_perturbnmf.h5ad
```

### Validation

The conversion script calls `external/PerturbNMF/.claude/skills/perturbNMF-runner/scripts/validate_data.py` automatically after writing. If validation fails:

- `'symbol' missing from mod['gene'].var` — older pipeline output; map a different column manually with `--gene-var-symbol-col <name>`.
- `'guide_assignment' missing from mod['guide'].layers` — guide assignment didn't run; re-check Stage 2.
- `gene var index not unique` — duplicate gene symbols. The script's option `--gene-var-symbol-col` can pick a different column (e.g., ensembl IDs) but PerturbNMF expects symbols downstream.

## Step 2: The new convention — `--compute_umap`

**Use this for new runs.** Each dataset's `<cnmf_run_name>/Script/Convert_file_adata.py` takes a `--compute_umap` flag that runs the standard scanpy pipeline on a temp copy and stashes `obsm['X_pca']` + `obsm['X_umap']` in the output AnnData before Stage 1 writes the inference inputs.

```bash
# Per-dataset, in the cNMF run dir's Script/
python Convert_file_adata.py \
  --in_h5ad <RUN>/cnmf/<cnmf_run_name>/Data/raw.h5ad \
  --out_h5ad <RUN>/cnmf/<cnmf_run_name>/Data/converted.h5ad \
  --compute_umap
```

The recipe is `normalize_total → log1p → HVG → scale → PCA → neighbors → UMAP` on a temp copy. PCA + UMAP propagate through every per-K h5mu produced by Stage 1.

Reference implementation: `datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/cleanser_800_mito_15pc/cnmf/050926_HTv2_20iter_5KHVG_torch_halsvar_batch/Script/Convert_file_adata.py`.

**When you can skip the post-Stage-1 prep entirely:** if Stage 1 input was produced with `--compute_umap`, the `inject_umap_into_h5mu` step is no longer needed (UMAP is already in every h5mu's rna modality).

## Step 3: Legacy prep (h5mu's produced before `--compute_umap`)

For DE / ESC / older runs, run these **after Stage 1 completes** and **before Stage 2a**:

### 3a. `prepare_h5mu_for_eval.py`

What it does:

1. Sets `obs['sample'] = 'all'` — single-value sample column makes perturbation tests pool across batches (not stratify).
2. Remaps NT `guide_targets`: `'nan'` → `'non-targeting'` so Stage 2b U-test fake-test path can find reference guides.

```bash
sbatch <run>/cnmf/<cnmf_run_name>/Script/prepare_h5mu_for_eval.sh
```

Runs on the selected-K h5mu (or all K h5mus; the modification is idempotent).

### 3b. `inject_umap_into_h5mu.py`

Pre-computes UMAP from the rna matrix and writes it into `mdata['rna'].obsm` AND `mdata['cNMF'].obsm`. Sidesteps Stage 3b's slow rna-PCA and Stage 3c's broken HVG call (upstream issue #9).

Recipe: standard scanpy (`normalize_total → log1p → HVG → scale → PCA → neighbors → UMAP`) on a copy of `mdata['rna']`.

```bash
sbatch <run>/cnmf/<cnmf_run_name>/Script/inject_umap_into_h5mu.sh
```

**Run only on the selected K's h5mu.** UMAP depends only on the rna matrix (identical across K), so doing it for one K is enough. Defer until Stage 3a has picked the K.

### 3c. `Data/guide_annotation.tsv`

Required by the U-test fake-test code path in Stage 2b. Generate from the canonical harmonized guide file:

```bash
cp ref/guide_libraries/harmonized/harmonized_guide_file_poolabcd.tsv \
   <run>/cnmf/<cnmf_run_name>/Data/guide_annotation.tsv
# then in the copy: rename column `guide_id` to `guide_names`
```

The column-rename can be done with awk/sed; the rest of the file passes through.

## When to use which path

| Scenario | Path |
|---|---|
| New dataset / new run starting from scratch | Step 1 (conversion) → Step 2 (`--compute_umap`) → Stage 1 onwards; skip Step 3 |
| Already-completed Stage 1 inference (DE / ESC / Hon CM seqspec_v3) | Step 1 (if not yet done) → Stage 1 (existing) → Step 3a/3b/3c → Stage 2 onwards |
| Patching an old h5mu mid-stream | Just Step 3a/3b/3c on the affected h5mu |

## Why these exist

Each prep step works around an upstream PerturbNMF bug or omission. Once upstream fixes land:

- Step 3a (single-sample remap) — fixed once perturbation tests support natural batch stratification.
- Step 3b (UMAP inject) — fixed by `--compute_umap` upstream of Stage 1 (us-side); no longer needed for new runs.
- Step 3c (guide_annotation.tsv) — upstream-required input; no fix planned.
