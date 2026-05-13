# DACC-format files — Gersbach HTv2 benchmark (testbed)

Per-dataset deliverables in the IGVF CPN FG file-format spec (see [`../../../guide/COMPLETE_DATASET_CONTENTS.md`](../../../guide/COMPLETE_DATASET_CONTENTS.md)).

This is the benchmark / testbed dataset — not part of the 5 production datasets. We may produce DACC files for it as a sanity check before mirroring production runs, but it isn't enshrined in the canonical schemas.

Library-wide DACC files (apply to every dataset):

- [`../../../reference/tf_universe.tsv`](../../../reference/tf_universe.tsv) — 1,951 TFs (pools A–D; the HTv2 testbed library is a subset — flag if you produce a dataset-specific TF universe).
- [`../../../reference/element_universe.bed`](../../../reference/element_universe.bed) — 2,260 promoter elements.

## Per-dataset files — status

Gene universe / programs / regulators: cNMF has been run end-to-end on this testbed (`cleanser_800_mito_15pc`). The DACC reformats haven't been produced yet because the testbed isn't being submitted to the portal — it's a structure-validation set.
