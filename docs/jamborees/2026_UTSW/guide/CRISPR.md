# Understanding CRISPR-pipeline outputs

The IGVF CRISPR pipeline ([github](https://github.com/IGVF/CRISPR_Pipeline)) takes raw single-cell data + a guide library and emits a structured set of QC, MuData, and differential-expression files. This doc explains what's in the output set and how to read each piece. Companion: [`ENERGY_DISTANCE.md`](ENERGY_DISTANCE.md), [`CNMF.md`](CNMF.md).

For the deeply technical pipeline walkthrough (run commands, GCP config, intermediate directories), see [`../../../analysis/crispr_pipeline/CRISPR_PIPELINE_OUTPUTS.md`](../../../analysis/crispr_pipeline/CRISPR_PIPELINE_OUTPUTS.md). This doc is the **interpretation** layer.

---

## The three folders you'll see

The CRISPR-pipeline output has three top-level folders:

| Folder | What it is | Size | Open it when… |
|---|---|---|---|
| `pipeline_dashboard/` | HTML dashboard + per-modality QC tables + figures | ~40 GB | …you want to see whether the experiment "worked." |
| `pipeline_outputs/` | The MuData + perturbo cis/trans differential-expression TSVs | ~23 GB | …you want per-perturbation results to analyze. |
| `pipeline_info/` | Pipeline run parameters + software versions | ~10 KB | …you need to reproduce or audit the run. |

We always upload all three together. If a dataset is missing `pipeline_info/`, the output set is incomplete — flag it.

---

## 1. The dashboard (`pipeline_dashboard/dashboard.html`)

Open this in a browser. It's a multi-page report covering:

### Cell counts
A table of how many cells passed each filter. Look for:
- Big drop-offs ("we started with 1M cells, ended with 200K"?). Some attrition is expected; >70% loss is a flag.
- Per-measurement-set balance — every batch should contribute roughly equal cell counts.

### Guide recovery
A bar/heatmap of how often each guide was successfully assigned to a cell. Look for:
- TFs with all guides ≥ 50 cells = robust.
- TFs with 0–5 cells per guide = unreliable; you can find them but don't trust the perturbation result.

### Knockdown QC (intended-target)
The fraction of cells with each guide where the intended target gene is suppressed vs. NTC cells. Look for:
- Median knockdown >50% — typical CRISPRi performance.
- AUROC on intended-target binary detection — typical values 0.65–0.85.

### Mapping QC
Read-alignment rates per modality:
- scRNA alignment: >85% expected for healthy 10x runs.
- Guide alignment: chemistry-dependent — 10x 3' v3 can be as low as 15%; 10x 5' and HT-like chemistries should be 60–95%.

### UMAPs
- Cells colored by measurement set — should look intermixed (no obvious batch effects).
- Cells colored by lineage/cell-type — should cluster as expected for the differentiation system.
- (Optional) Cells colored by guide identity for a chosen TF — strong perturbations will pull cells off the main cluster.

If any of these are off, ask before treating the downstream results as biology.

---

## 2. The inference MuData (`pipeline_outputs/inference_mudata.h5mu`)

The "every cell + every guide + every gene" object. See [`DATA_FORMATS.md`](DATA_FORMATS.md) for how to open. Modalities:

- `mdata.mod["gene"]` — RNA. `obs` has per-cell metadata; `var` has per-gene metadata.
- `mdata.mod["guide"]` — guides. `obs` has the same cells; `X` is cell × guide assignment (binary or counts depending on guide-assignment method).
- `mdata.mod["hashing"]` — optional, HTO multiplexing.

Important `obs` columns:
- `intended_target_name` — the gene symbol of the TF a cell's guide targets. NTC cells get `non_targeting`.
- `intended_target_id` — ENSG of the same.
- `guide_id` — the specific guide barcode the cell carries.
- Lineage/cell-type calls if computed.

You won't normally crack open the MuData unless you're doing a custom analysis — the perturbo TSVs cover most questions.

---

## 3. The perturbo differential-expression TSVs (`pipeline_outputs/`)

These are the **headline-numbers** file for the experiment.

Each TSV has one row per (perturbation, gene) pair. Two flavors:

| File | Rows | What it tells you |
|---|---|---|
| `perturbo_cis_per_element.tsv.gz` | n_perturbations × n_target_genes | **Did the knockdown work?** Compares the intended target gene's expression in guide-carrying cells vs NTC cells. Used for QC. |
| `perturbo_trans_per_element.tsv.gz` | n_perturbations × n_all_genes | **Where did the knockdown's effect propagate?** Genome-wide log2FC + p-value per perturbation. The headline file. |

Key columns:

- `intended_target_name` — TF being knocked down.
- `gene_id` / `gene_symbol` — the gene being tested.
- `log2_fc` — log2 fold change. Negative = the gene went down; positive = up.
- `p_value` — raw p-value from perturbo.
- `fdr_bh` — Benjamini-Hochberg FDR (per-perturbation). The right column to threshold on (usually `fdr_bh < 0.05`).

### "How do I find significant trans targets for my TF?"

```python
import pandas as pd
trans = pd.read_csv("perturbo_trans_per_element.tsv.gz", sep="\t")
my_tf = trans[(trans["intended_target_name"] == "<MY_TF>") & (trans["fdr_bh"] < 0.05)]
print(my_tf.shape)
print(my_tf.sort_values("log2_fc").head(20))   # most-down-regulated
```

### "How many TFs gave any signal at all?"

A useful first cut: count TFs with ≥10 trans targets at FDR < 0.05.

```python
counts = (trans[trans["fdr_bh"] < 0.05]
          .groupby("intended_target_name").size())
print((counts >= 10).sum())     # n TFs with detectable trans signal
```

---

## 4. The pipeline-info folder (`pipeline_info/`)

Two files:

- `params_<run>.json` — every Nextflow parameter used (chemistry, MOI, references, software versions).
- `software_versions.yml` — exact versions of perturbo, sceptre, kallisto, etc.

You'll only consult these when reproducing a run or comparing across pipeline versions.

---

## Issues

### Guide alignment rates are chemistry-dependent
Don't compare guide-alignment percentages across datasets that used different 10x chemistries. The pipeline reports them honestly; a 16% guide-alignment rate on 10x 3' v3 isn't a failure.

### Perturbo p-values are uncalibrated by default
The raw perturbo p-values come from a model assumption that may not hold for every dataset. The project maintains a calibration pipeline (`src/tf_perturb_seq/inference/calibrate.py`) that re-derives empirical p-values from NTC cells. If the dataset has a `calibrated_*` TSV alongside the raw perturbo TSV, **use the calibrated one**.

### Trans-effect density varies by pipeline version + cell line
A newer pipeline version (e.g. `seqspec_v3`) or a different cell line can shift the absolute number of significant trans targets per TF by 5–10×. Be careful before reading absolute counts as biology — focus on relative rankings within the same dataset.

---

## What to look for (the takeaways)

These are the panels and rankings most useful to highlight in a presentation:

1. **The TFs with the broadest trans effects** in the dataset. These are candidate lineage-defining TFs. Pull the top 10–20 by trans-target count.
2. **Whether your favorite TF actually moved its target gene** (cis QC). If not, the trans signal isn't trustworthy.
3. **Whether you can recapitulate a known biology result** — e.g. knocking down a lineage-defining TF in its on-state lineage should down-regulate that lineage's marker genes. If yes, the experiment is grounded; if not, ask.
4. **TFs with strong trans effects but unexpected functions** — candidates for deorphanization. Cross-reference the trans-DE table with TF annotations (Lambert 2018, JASPAR).
5. **Reproducibility across guides per TF** — if 4 guides for the same TF agree, the result is robust. If 1 of 4 disagrees, it's likely an off-target effect.

---
