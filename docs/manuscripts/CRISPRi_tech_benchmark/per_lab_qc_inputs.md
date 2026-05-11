# Per-lab QC inputs — unfiltered scRNA h5ads on GCS

**Goal:** distribute the raw, unfiltered scRNA-seq h5ad for each benchmark dataset to its source lab so each lab can run QC locally and propose appropriate filter parameters. Once labs converge on params, the CRISPR pipeline runs end-to-end with those fixed — saving compute cost and turnaround.

The pointer table is in [per_lab_qc_inputs.tsv](per_lab_qc_inputs.tsv).

---

## What's at each pointer (verified 2026-04-29)

For each dataset, two GCS paths are provided. **Default is `scrna_h5ad_unfiltered_gs`** — that's the raw, truly pre-filter scRNA matrix.

| File | Modality | Cells (Hon example) | Pre-filter status |
|---|---|---:|---|
| `anndata/concatenated_adata.h5ad` | gene only (scRNA) | **2,131,101** | ✅ TRULY pre-filter — all kallisto-bustools barcodes, no cell-calling / mito / scrublet applied. Includes empty droplets. |
| `createmudata/mudata.h5mu` | gene + guide | varies (101k @ 800 UMI) | ⚠ Already POST cell-call filter — same cells as final inference_mudata. Use for guide-aware QC, not cell-calling exploration. |
| ❌ `preprocessanndata/filtered_anndata.h5ad` | gene only | 125,971 (Hon) | ❌ POST-FILTER (UMI ≥ 800, mito ≤ 15%). Do NOT use as a "pre-QC" starting point — the filter floor is baked in. |

**Verified for Hon @ cleanser_800_mito_15pc:**
- Raw `anndata/concatenated_adata.h5ad`: 2.13M cells × 62,757 genes (sparse CSR, float64, ~38× the post-filter cell count of 56K)
- Post-filter `filtered_anndata.h5ad`: 126K cells × 28,183 genes; min total_counts = 968, max pct_counts_mt = 15.0%

## What the raw h5ad has and doesn't have

**Has:**
- Cell × gene UMI count matrix (sparse CSR, float64)
- All cell barcodes from kallisto-bustools mapping (including empty droplets)
- 3 batches per dataset (one per measurement set / lane)

**Does NOT have (labs would need to compute):**
- Gene symbols — `var.index` is Ensembl gene IDs (`ENSG...`); no `symbol` column. Labs need to map ENSG → symbol themselves (e.g., via biomart, pyensembl, or the GTF in `ref/`)
- MT / ribo flags — no boolean columns; labs need to identify mito genes themselves (after symbol mapping, by `^MT-` prefix). For Hon there are 26 mito genes in the gencode reference.
- Per-cell QC metrics — no `total_counts`, `n_genes_by_counts`, `pct_counts_mt` etc. Run `sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'])` after symbol/flag mapping.

## Suggested per-lab QC workflow

```python
import anndata as ad
import scanpy as sc

# 1. Download (auth as your project credentials or igvf-pertub-seq-pipeline service account)
# gsutil cp gs://...Benchmark_cleanser_800_mito_15pc/<DATASET>/anndata/concatenated_adata.h5ad ./

# 2. Load
a = ad.read_h5ad('concatenated_adata.h5ad')
# a.shape = (2.13M, 62757) for Hon; var.index is Ensembl IDs

# 3. Map gene symbols (e.g., from a GTF or pyensembl)
#    ensembl_to_symbol = {...}  # build from GTF
#    a.var['symbol'] = a.var.index.map(ensembl_to_symbol)
#    a.var['mt'] = a.var['symbol'].str.startswith('MT-')
#    a.var['ribo'] = a.var['symbol'].str.startswith(('RPL', 'RPS'))

# 4. Compute QC metrics
sc.pp.calculate_qc_metrics(a, qc_vars=['mt', 'ribo'], inplace=True, percent_top=[50])

# 5. Explore filter choices
#   - knee plot on a.obs['total_counts']
#   - flat UMI floors: 200 / 500 / 800 / 2000
#   - mito threshold sweep (10/15/20/25%)
#   - scrublet (sc.external.pp.scrublet)

# 6. Report back recommended filter parameters in shared doc / GitHub issue.
```

## Source pipeline run

All paths are taken from the canonical `Benchmark_cleanser_800_mito_15pc` run on GCS. The upstream `anndata/concatenated_adata.h5ad` stage is functionally identical across the parameter-sweep configs since the cell-calling / mito / scrublet differences happen later in the pipeline.

## Why this matters

Each full pipeline rerun on GCP costs ~$25/dataset × 5 datasets = ~$125 per parameter-sweep config. By having labs run QC on the unfiltered h5ad locally, we can iterate on filter choices for free and only pay for pipeline runs we actually want as canonical. The current parameter sweep on GCS (200, 500, 800, 2000, knee2 + scrublet pair = 7+ runs × 5 datasets) cost on the order of ~$1k. Distributing the unfiltered data eliminates the need for that kind of brute-force sweep.

## Open question — should we pre-process for labs?

The raw h5ad needs gene-symbol mapping + QC metric computation before lab can do QC. Options:

1. **Ship as-is** with the workflow above — labs do the symbol mapping themselves (~30 min of scanpy)
2. **Pre-process and re-host** — Adam runs `calculate_qc_metrics` + symbol mapping once per dataset, hosts a `lab_qc_input.h5ad` (~3-4 GB each) with symbols + mt/ribo flags + per-cell QC metrics + no cell filter applied. Labs get a one-step `ad.read_h5ad` → ready to explore.

Option 2 is more polished (~1 hour to produce all 5 datasets) and probably the right move if we expect 5+ labs to consume this. Option 1 is faster to ship.

## Reference: dataset overview

See [dataset_overview.tsv](dataset_overview.tsv) for biological system, CRISPRi construct, sgRNA delivery, sequencing tech, and contact per dataset.

## GCS access notes (verified 2026-04-29)

`gsutil iam get gs://igvf-pertub-seq-pipeline-data` shows access is gated by membership in the IGVF Perturb-Seq GCP project:
- **storage.admin** (full access) — Chikara (`chikarabs52@gmail.com`), Gary, Isha, Mohsen Razavi, Olga Pushkarev, Weizhou Qian, Zahra Heidary
- **legacyBucketReader / legacyObjectReader** (read) — anyone with `projectViewer:igvf-pertub-seq-pipeline`

**Per-lab access status:**

| Lab | Contact | GCS access? | Action needed |
|---|---|---|---|
| Hon | Chikara Takeuchi (`chikarabs52@gmail.com`) | ✅ already storage.admin | None — can pull directly |
| Huangfu | Nan Zhang | ❌ not in IAM | Add to project as projectViewer, or upload h5ad to Synapse |
| Gersbach (×2) | Boxun Li | ❌ not in IAM | Same |
| Engreitz | Tri Nguyen | ❌ not in IAM | Same |

**Access strategy options:**
1. **Add to GCP project** — Lucas / Gary grant `projectViewer` role to each lab contact's Google account. Cleanest long-term.
2. **Upload to Synapse** — copy the per-dataset h5ads into the consortium Synapse project (`syn63675917`). Most labs already have Synapse access. Larger upload (~30 GB total) but no IAM friction.
3. **Signed URLs** — `gsutil signurl` for time-limited access. OK for a one-off but not for iterative QC.

**Recommendation**: option 1 for the Hon lab (Chikara has access). Option 2 (Synapse) for the others — same pattern that's already used for Hon cardio + Gersbach hep production h5mu uploads.
