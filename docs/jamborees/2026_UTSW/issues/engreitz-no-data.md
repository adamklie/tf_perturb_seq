# Issue 4 — Engreitz endothelial data not on portal

**Status**: ☐ blocked — no fastqs / no inference MuData / no nothing for `Engreitz_WTC11-endothelial-cells_TF-Perturb-seq`. Every downstream output is gated on this.

**Owner of fix**: Engreitz team (data upload to IGVF portal).

## TL;DR

The Engreitz endothelial dataset is one of the 5 production datasets for the jamboree, but as of the latest IGVF portal snapshot (2026-05-07), no measurement sets exist for it under the `TF Perturb-seq Project` collection. We can't run any of the three pipelines (CRISPR, cNMF, energy distance) without raw data → inference MuData → analysis.

## Evidence

[`scripts/query_igvf_portal.py`](../scripts/query_igvf_portal.py) snapshots the IGVF portal state into `portal_snapshots/<utc-iso>/`. The most recent snapshot's `manifest.tsv` has entries for Hon CM (and only Hon CM under the `TF Perturb-seq Project` collection) — no Engreitz.

```bash
# Confirm current portal state:
.venv/bin/python docs/jamborees/2026_UTSW/scripts/query_igvf_portal.py
# Then check portal_snapshots/latest/manifest.tsv for Engreitz rows.
```

[`docs/jamborees/2026_UTSW/README.md`](../README.md) — Portal snapshots note records the same: "as of the first snapshot (2026-05-07), only Hon measurement sets currently appear under the `TF Perturb-seq Project` collection on the portal."

## What we want

The Engreitz team to upload the WTC11 endothelial Perturb-seq data to the IGVF portal under the `TF Perturb-seq Project` collection. Specifically:

- Raw fastqs, properly metadata'd into:
  - One or more **measurement sets** (scRNA-seq, split by cellular sub-pools)
  - One **auxiliary set** per measurement set (gRNA-seq)
  - Optionally a second auxiliary set per measurement set (HTO-seq)

The data layout matches what we did for the other production datasets (see [`docs/analysis/PIPELINES.md`](../../../analysis/PIPELINES.md) → Stage 1).

## Acceptance criteria

- [ ] Engreitz endothelial measurement sets visible in [`scripts/query_igvf_portal.py`](../scripts/query_igvf_portal.py) snapshot under `TF Perturb-seq Project`.
- [ ] [`reference/experimental_metadata.tsv`](../reference/experimental_metadata.tsv) row for Engreitz updated with IGVF accessions, `gcs_output_path`, `canonical_run_label`, etc.
- [ ] CRISPR pipeline can be launched (Stage 2 of [`docs/analysis/PIPELINES.md`](../../../analysis/PIPELINES.md)).

After that the rest cascades: CRISPR → MuData → cNMF + energy distance.

## Pointers

| Object | Path |
|---|---|
| Portal snapshot tool | [`scripts/query_igvf_portal.py`](../scripts/query_igvf_portal.py) |
| Portal snapshots dir | [`portal_snapshots/`](../portal_snapshots/) (rolling history) |
| Pipeline stages overview | [`docs/analysis/PIPELINES.md`](../../../analysis/PIPELINES.md) |
| Experimental metadata generator | [`scripts/generate_experimental_metadata.py`](../scripts/generate_experimental_metadata.py) |
| Experimental metadata TSV | [`reference/experimental_metadata.tsv`](../reference/experimental_metadata.tsv) |
