# Jamboree data — for experimentalists

A 1-page map of the 5 production datasets in the 2026 UTSW jamboree, with one-click access to each output. Built for someone who wants to look at results without writing code.

> **Status:** v0 draft (2026-05-11). Will fill out once more outputs land on Synapse.

---

## At-a-glance status

| Dataset | CRISPR pipeline | cNMF | Energy distance |
|---|:---:|:---:|:---:|
| Hon WTC11 Cardiomyocyte | ⚠ partial | ☐ blocked | ✅ |
| Huangfu HUES8 Definitive Endoderm | ✅ | ✅ | ⚠ calibration |
| Huangfu HUES8 Embryonic Stem Cell | ✅ | ✅ | ⚠ calibration |
| Gersbach WTC11 Hepatocyte | ⚠ partial | ☐ blocked | ☐ blocked |
| Engreitz WTC11 Endothelial | ☐ blocked | ☐ blocked | ☐ blocked |

**Legend** — ✅ complete, on Synapse, ready to use • ⚠ on Synapse with a known caveat (read the dataset card before relying on it) • ☐ blocked / not yet available

---

## Dataset cards

### 🫀 Hon WTC11 Cardiomyocyte

- **Lab:** Hon (UTSW) · **Cell line:** WTC11 · **Differentiation:** Cardiomyocyte (12-day)
- **Tech:** 10x 5'-Perturb (Sigma backbone, HT-like chemistry) · **Multiplexing:** HTO
- **Library:** TF guides pools A–D + F · **Measurement sets:** 28
- **Perturbation:** CRISPRi

| Output | Status | Synapse | Notes |
|---|:---:|---|---|
| CRISPR pipeline | ⚠ partial | [`syn74520421`](https://www.synapse.org/Synapse:syn74520421) | Bundle missing `pipeline_info/`; we're re-running upstream pipeline ourselves |
| cNMF | ☐ blocked | — | Gated on full CRISPR bundle |
| Energy distance | ✅ | [`syn74897350`](https://www.synapse.org/Synapse:syn74897350) | Calibration looks healthy |

---

### 🧬 Huangfu HUES8 Definitive Endoderm

- **Lab:** Huangfu (MSKCC) · **Cell line:** HUES8 · **Differentiation:** Definitive endoderm
- **Tech:** 10x 3' v3 · **Multiplexing:** none
- **Library:** TF guides pools A–D · **Measurement sets:** 8
- **Perturbation:** CRISPRi

| Output | Status | Synapse | Notes |
|---|:---:|---|---|
| CRISPR pipeline | ✅ | [`syn74834952`](https://www.synapse.org/Synapse:syn74834952) | `muddy_penguin` run (13 bp spacer_tag) |
| cNMF | ✅ | [`syn74893844`](https://www.synapse.org/Synapse:syn74893844) | |
| Energy distance | ⚠ calibration | [`syn74883327`](https://www.synapse.org/Synapse:syn74883327) | **All 100 NCs have `pval_mean=0` — do not threshold on p-value alone.** Use `distance_mean` as effect-size proxy; the run is otherwise structurally valid. |

---

### 🧫 Huangfu HUES8 Embryonic Stem Cell

- **Lab:** Huangfu (MSKCC) · **Cell line:** HUES8 · **Differentiation:** Embryonic stem cell (undifferentiated)
- **Tech:** 10x 3' v3 · **Multiplexing:** none
- **Library:** TF guides pools A–D · **Measurement sets:** 8
- **Perturbation:** CRISPRi

| Output | Status | Synapse | Notes |
|---|:---:|---|---|
| CRISPR pipeline | ✅ | [`syn74835010`](https://www.synapse.org/Synapse:syn74835010) | `sceptre_v1` run |
| cNMF | ✅ | [`syn74893846`](https://www.synapse.org/Synapse:syn74893846) | |
| Energy distance | ⚠ calibration | [`syn74883475`](https://www.synapse.org/Synapse:syn74883475) | Same calibration caveat as the DE sibling above |

---

### 🫁 Gersbach WTC11 Hepatocyte

- **Lab:** Gersbach (Duke) · **Cell line:** WTC11 · **Differentiation:** Hepatocyte (22-day)
- **Tech:** 10x Perturb-seq (NovaSeq X Plus, 25B; consistent with 10x 3' v3) · **Multiplexing:** none
- **Library:** TF guides pools A–D + F · **Measurement sets:** 47
- **Perturbation:** CRISPRi

| Output | Status | Synapse | Notes |
|---|:---:|---|---|
| CRISPR pipeline | ⚠ partial | [`syn70518849`](https://www.synapse.org/Synapse:syn70518849) | Non-canonical layout; awaiting canonical bundle from Gersbach team |
| cNMF | ☐ blocked | — | Awaiting Gersbach team's canonical run |
| Energy distance | ☐ blocked | — | Awaiting Gersbach team's canonical run |

---

### 🩸 Engreitz WTC11 Endothelial

- **Lab:** Engreitz (Stanford) · **Cell line:** WTC11 · **Differentiation:** Endothelial
- **Tech:** CC Perturb-seq · **Multiplexing:** TBD
- **Library:** TBD · **Measurement sets:** 0 (not yet uploaded)
- **Perturbation:** CRISPRi

| Output | Status | Synapse | Notes |
|---|:---:|---|---|
| CRISPR pipeline | ☐ blocked | — | Raw data not yet on the IGVF portal |
| cNMF | ☐ blocked | — | Gated on CRISPR pipeline |
| Energy distance | ☐ blocked | — | Gated on CRISPR pipeline |

---

## Where to go next

- **Cross-dataset summary tables** (the most useful starting point if you want numbers fast):
  - [`reference/cross_dataset_pipeline_summary.tsv`](reference/cross_dataset_pipeline_summary.tsv) — cell counts, UMI medians, knockdown stats per dataset
  - [`reference/cross_dataset_edistance_summary.tsv`](reference/cross_dataset_edistance_summary.tsv) — target counts, distance medians, calibration-robust significance flags
  - [`reference/tf_metadata_simplified.tsv`](reference/tf_metadata_simplified.tsv) — every TF in the library with HGNC + Lambert + JASPAR annotations
  - [`reference/experimental_metadata_simplified.tsv`](reference/experimental_metadata_simplified.tsv) — per-dataset wet-lab metadata
- **Browsing on Synapse**: click any of the `syn...` links above. Synapse will prompt you to log in; once you do you can browse the bundle in your browser and download files individually.
- **If you want to know what each output means biologically** — [TBD section, will land once interpretation guide is written]
- **If you're stuck or anything's confusing** — drop a comment on [Issue #26](https://github.com/adamklie/tf_perturb_seq/issues/26) (the master TODO) or ping Adam directly.
