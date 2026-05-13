# DACC-format files — Hon WTC11 Cardiomyocyte

Per-dataset deliverables in the IGVF CPN FG file-format spec (see [`../../../guide/COMPLETE_DATASET_CONTENTS.md`](../../../guide/COMPLETE_DATASET_CONTENTS.md)).

Library-wide DACC files (apply to every dataset):

- [`../../../reference/tf_universe.tsv`](../../../reference/tf_universe.tsv) — 1,951 TFs.
- [`../../../reference/element_universe.bed`](../../../reference/element_universe.bed) — 2,260 promoter elements.

## Per-dataset files — status

| Spec | Source | Status |
|---|---|---|
| `Gene Universe` | cNMF `Inference.overdispersed_genes.txt` | ☐ blocked — production cNMF not yet run (gated on Hon CM CRISPR bundle from Weizhou). |
| `Gene Programs` | cNMF `Inference.gene_spectra_score.k_<sel>.dt_2_0.txt` | ☐ blocked — same gate. |
| `Gene Program Regulators` | cNMF `<sel>_perturbation_association_results_all.txt` | ☐ blocked — same gate. |
| `Global differential expression` | `perturbo_trans_per_element.tsv.gz` | ⏳ runnable once Hon CRISPR pipeline re-finishes. Reference Hon submission: [`IGVFFI5989UAVX`](https://data.igvf.org/tabular-files/IGVFFI5989UAVX/) (pySpade-flavor). |
| `Local differential expression` | `perturbo_cis_per_element.tsv.gz` | ⏳ runnable once Hon CRISPR pipeline re-finishes. |

Reference IGVF submission for this dataset: [IGVFDS6332VCTO](https://data.igvf.org/analysis-sets/IGVFDS6332VCTO/).

## Regenerating

Once cNMF + CRISPR pipeline land, run the generators in [`src/tf_perturb_seq/dacc/`](../../../../../../src/tf_perturb_seq/dacc/). See the [Huangfu DE dacc README](../../Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/dacc/README.md) for the command template.
