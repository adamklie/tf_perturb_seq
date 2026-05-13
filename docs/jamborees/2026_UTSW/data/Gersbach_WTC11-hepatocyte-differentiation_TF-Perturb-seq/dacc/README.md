# DACC-format files — Gersbach WTC11 Hepatocyte

Per-dataset deliverables in the IGVF CPN FG file-format spec (see [`../../../guide/COMPLETE_DATASET_CONTENTS.md`](../../../guide/COMPLETE_DATASET_CONTENTS.md)).

Library-wide DACC files (apply to every dataset):

- [`../../../reference/tf_universe.tsv`](../../../reference/tf_universe.tsv) — 1,951 TFs.
- [`../../../reference/element_universe.bed`](../../../reference/element_universe.bed) — 2,260 promoter elements.

## Per-dataset files — status

| Spec | Source | Status |
|---|---|---|
| `Gene Universe` | cNMF `Inference.overdispersed_genes.txt` | ☐ blocked — awaiting canonical bundle from Gersbach team. |
| `Gene Programs` | cNMF `Inference.gene_spectra_score.k_<sel>.dt_2_0.txt` | ☐ blocked — same gate. |
| `Gene Program Regulators` | cNMF `<sel>_perturbation_association_results_all.txt` | ☐ blocked — same gate. |
| `Global differential expression` | `perturbo_trans_per_element.tsv.gz` | ☐ blocked — awaiting canonical CRISPR bundle from Sara (Gersbach team). |
| `Local differential expression` | `perturbo_cis_per_element.tsv.gz` | ☐ blocked — same gate. |

See [the project issues](https://github.com/adamklie/tf_perturb_seq/issues) for the full delivery checklist.
