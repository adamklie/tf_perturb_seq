# DACC-format files — Huangfu HUES8 Definitive Endoderm

Per-dataset deliverables in the IGVF CPN FG file-format spec (see [`../../../guide/COMPLETE_DATASET_CONTENTS.md`](../../../guide/COMPLETE_DATASET_CONTENTS.md) for the full inventory + interpretation guidance).

## What's here

| File | Spec | Source | Notes |
|---|---|---|---|
| [`gene_universe.tsv`](gene_universe.tsv) | DACC `Gene Universe` | cNMF `Inference.overdispersed_genes.txt` from `muddy_penguin` run | 2,000 HVGs, all resolved to ENSG (0 unmapped). |
| [`gene_universe.misses.tsv`](gene_universe.misses.tsv) | (diagnostic) | side-output of the generator | Should be empty; if it's not, those rows blocked portal submission. |

Library-wide DACC files (apply to every dataset using the TF Perturb-seq pool A–D library):

- [`../../../reference/tf_universe.tsv`](../../../reference/tf_universe.tsv) — 1,951 TFs.
- [`../../../reference/element_universe.bed`](../../../reference/element_universe.bed) — 2,260 promoter elements.

## Still to land for this dataset

| Spec | Source | Status |
|---|---|---|
| `Gene Programs` | `Inference.gene_spectra_score.k_<sel>.dt_2_0.txt` | Reformat queued; depends on group-selected k. |
| `Gene Program Regulators` | `<sel>_perturbation_association_results_all.txt` | Reformat queued; same gate. |
| `Global differential expression` | `perturbo_trans_per_element.tsv.gz` | Reformat queued (pySpade or generic spec — decision in [`../../../../data/DACC.md`](../../../../data/DACC.md)). |
| `Local differential expression` | `perturbo_cis_per_element.tsv.gz` | Reformat queued. |

## Regenerating

```bash
PYTHONPATH=src uv run python -m tf_perturb_seq.dacc.build_gene_universe \
    --hvg datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/muddy_penguin/cnmf/Result/Inference/Inference.overdispersed_genes.txt \
    --gtf ref/genome/IGVFFI9573KOZR.gtf.gz \
    --out docs/jamborees/2026_UTSW/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/dacc/gene_universe.tsv \
    --misses-out docs/jamborees/2026_UTSW/datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/dacc/gene_universe.misses.tsv
```
