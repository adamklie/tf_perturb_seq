# Overview

Prepare data and documentation for the 2026 UTSW jamboree. We are working with the 5 production datasets shown in `state/2026_05_07_state.png` (this directory).

**Where we are (2026-05-12, day before jamboree):** Steps 1, 2, 3 done. **Step 6** (Synapse uploads) covers 3 of 5 production datasets across CRISPR pipeline + energy distance + QC; Hon CM CRISPR awaits Weizhou's bundle, Gersbach Hep awaits Sara's deliverables, Engreitz blocked on portal. **Step 5** human-readable outputs landed for WG1 (qc / e-dist / tf_cross_lineage), WG3 (disease / convergence), WG4 (network structure / edges), WG5 (family scorecard) — see [`working_groups/`](working_groups/). cNMF cross-dataset outputs (WG2) still gated on production runs. **Step 7** (Google Sheet + final docs) still open.

Old daily AGENDA notes archived under `scratch/2026_05_12/UTSW_Jamboree/`; the active plan now lives in this file + the [`issues/`](issues/) reports.

# Process notes

- **Always double-check paths before using them.** Confirm whether a path is local, on GCS, on Synapse, or on the HPC.
- **A lot of the source data lives on the UCSD HPC.** Connection can be slow, so expect to organize things over SSH:
  - Host: `aklie@nrnb-login.ucsd.edu`
  - Project root: `/cellar/users/aklie/projects/tf_perturb_seq`
- When in doubt about which dataset, repo, or platform a task targets, ask before acting.
- **GitHub Project**: [TFP3](https://github.com/users/adamklie/projects/4) — milestones as draft items, issues as tasks
- **Synapse-as-we-go**: We upload each artifact to Synapse **directly from where it lives** (HPC / GCS / IGVF portal) as soon as it's ready, and log the path in `synapse_paths.tsv`. We do **not** stage everything locally first. This folder only holds docs + small simplified outputs.

# Steps

## 1. Set up the jamboree folder in this repo ✅

Working dir: `tf_perturb_seq/docs/jamborees/2026_UTSW/`

- [x] Build a TSV capturing the state shown in `state/2026_05_07_state.png` → `state/2026_05_07_state.tsv`.
- [x] Start a `README.md` in this folder that we extend as decisions get made → live in [`README.md`](README.md), kept up to date.

## 2. Decide on local and cloud organization ✅

- [x] Local layout: `datasets/<dataset>/<analysis>/` (per-analysis subdir per dataset).
- [x] Matching Synapse layout: `syn64423137/2026_UTSW/datasets/<dataset>/<analysis>/`.

## 3. Define the machine-readable outputs ✅

Schemas live in [`schemas/`](schemas/) (one JSON per output table; index in [`schemas/README.md`](schemas/README.md)). Documented in [`README.md`](README.md).

**Reference data**
- [x] IGVF GTF — single canonical file (no schema; standard format)
- [x] TF metadata — [`schemas/tf_metadata.json`](schemas/tf_metadata.json)
- [x] Experimental metadata — [`schemas/experimental_metadata.json`](schemas/experimental_metadata.json)
- [x] Guide metadata — [`schemas/guide_metadata.json`](schemas/guide_metadata.json) (used as-is from the IGVF portal)

**Perturb-seq outputs**
- [x] Inference MuData → covered as part of `crispr_pipeline/pipeline_dashboard/inference_mudata.h5mu` ([`schemas/crispr_pipeline.json`](schemas/crispr_pipeline.json))

**cNMF outputs**
- [x] [`schemas/cnmf.json`](schemas/cnmf.json) — curation rule = selected-k full data + sweep-as-provenance (revisitable without re-running)

**Energy distance outputs**
- [x] [`schemas/energy_distance.json`](schemas/energy_distance.json) — verified against HTv2 reference run ([`syn74381167`](https://www.synapse.org/Synapse:syn74381167))

## 4. Define the human-readable outputs

Outputs aimed at general scientists — simpler artifacts for higher-level figures and exploration.

**Reference data**
- [x] Simplified TF metadata → `reference/tf_metadata_simplified.tsv`
- [x] Simplified experimental metadata → `reference/experimental_metadata_simplified.tsv`
- [x] Simplified guide metadata → not needed (IGVF release used as-is)

**Perturb-seq outputs**
- [x] Cross-dataset summary TSV (cell counts, UMI medians, knockdown stats, perturbo significance counts) → [`reference/cross_dataset_pipeline_summary.tsv`](reference/cross_dataset_pipeline_summary.tsv) (3 datasets so far: Hon CM, Huangfu DE, Huangfu ESC; Gersbach Hep skipped because of non-canonical Synapse layout). Generator: [`src/tf_perturb_seq/crispr_pipeline/cross_dataset_pipeline_summary.py`](../../../src/tf_perturb_seq/crispr_pipeline/cross_dataset_pipeline_summary.py). Re-run when more datasets get canonical bundles.

**cNMF outputs** (deferred — depend on production cNMF runs)
- [ ] Cross-dataset program-similarity heatmap (cosine similarity of `gene_spectra_score` across all 5 datasets at each dataset's selected k)
- [ ] Per-dataset top-20-genes-per-program TSV
- [ ] Per-dataset regulators-per-program TSV (from `Eval/<sel>_<dt>/<sel>_perturbation_association_results_*.txt`)

**Energy distance outputs** (deferred — partial data so far)
- [ ] Cross-dataset TFs-significant-by-FDR summary TSV (rolled up from `pval_edist_full.csv`)
- [ ] Cross-dataset clustered heatmap of per-target energy distances (joins `target_by_target_matrix.csv` on shared targets)

Storage considerations: simplified TSVs (small, ~KB-MB) live in this repo under `reference/` or `datasets/<name>/<analysis>/`. Larger derived artifacts (heatmap PNGs, cross-dataset MuData) go to Synapse alongside the comprehensive outputs.

## 5. Stage simplified outputs locally

Only the small, human-readable artifacts (e.g., `*_simplified.tsv` in `reference/`) live in this repo. Bulky files do not get staged locally — they go straight to Synapse.

- [x] Drop simplified reference tables (TF / experimental) into `reference/`.
- [x] Per-dataset READMEs under `datasets/<name>/` (2026-05-09; top-level + per-analysis subdirs for all 5 production datasets + HTv2 testbed). Will be refreshed in place as each output's status changes — see [Issue 6](issues/per-dataset-readmes.md).
- [ ] Add small simplified summaries per analysis where they make sense (deferred — depends on Step 4 outputs landing).

## 6. Upload to Synapse as artifacts come in (interleaved with 3 & 4)

Done iteratively, not at the end. Tracking lives in [`synapse_paths.tsv`](synapse_paths.tsv).

**Reference data** ✅
- [x] TF metadata → [`syn74834227`](https://www.synapse.org/Synapse:syn74834227)
- [x] Experimental metadata → [`syn74834309`](https://www.synapse.org/Synapse:syn74834309)
- [x] IGVF GTF → [`syn74834518`](https://www.synapse.org/Synapse:syn74834518)
- [x] Guide library → [`syn74834519`](https://www.synapse.org/Synapse:syn74834519)

**CRISPR pipeline** (mirror script: [`scripts/mirror_pipeline_outputs.py`](scripts/mirror_pipeline_outputs.py) for GCS source, [`scripts/mirror_pipeline_outputs_hpc.py`](scripts/mirror_pipeline_outputs_hpc.py) for HPC source)
- [x] Hon WTC11 Cardiomyocyte — Weizhou added `pipeline_info/` to [`syn74520421`](https://www.synapse.org/Synapse:syn74520421); full bundle copied into jamboree `crispr_pipeline/` → [`syn74919102`](https://www.synapse.org/Synapse:syn74919102) (2026-05-12, flat layout matching DE/ESC).
- [x] Huangfu HUES8 Definitive Endoderm → [`syn74834952`](https://www.synapse.org/Synapse:syn74834952)
- [x] Huangfu HUES8 Embryonic Stem Cell → [`syn74835010`](https://www.synapse.org/Synapse:syn74835010)
- [ ] Gersbach WTC11 Hepatocyte — non-canonical at [`syn70518849`](https://www.synapse.org/Synapse:syn70518849); awaiting canonical bundle from **Sara** (Gersbach team)
- [ ] Engreitz WTC11 Endothelial — ☐ blocked: no data on the IGVF portal yet

**cNMF** (mirror script: [`scripts/mirror_cnmf_outputs.py`](scripts/mirror_cnmf_outputs.py); takes `--selected-k` and applies the curation rule from `schemas/cnmf.json`)
- [ ] HTv2 testbed (job 10577039 running) — verify pipeline structure end-to-end before launching production runs
- [ ] Hon WTC11 Cardiomyocyte — gated on full CRISPR bundle from Hon team (Weizhou)
- [ ] Huangfu HUES8 Definitive Endoderm — runnable now; awaits group k-selection
- [ ] Huangfu HUES8 Embryonic Stem Cell — runnable now; awaits group k-selection
- [ ] Gersbach WTC11 Hepatocyte — bug **Sara** to deliver in our format ([`schemas/cnmf.json`](schemas/cnmf.json)); Sara likely has a complete run, we just need it shaped to match our curation rule (selected-k bundle + sweep-as-provenance) and uploaded under `2026_UTSW/datasets/<id>/cnmf/`
- [ ] Engreitz WTC11 Endothelial — blocked: no data
- [ ] Group k-selection meeting per dataset (clinical-review-board style per [`docs/analysis/cNMF.md`](../../analysis/cNMF.md))

**Energy distance** (mirror script: [`scripts/mirror_edistance_outputs.py`](scripts/mirror_edistance_outputs.py))
- [x] Hon WTC11 Cardiomyocyte — both Adam's run → [`syn74897350`](https://www.synapse.org/Synapse:syn74897350) and Sara's comparison run → [`syn74910330`](https://www.synapse.org/Synapse:syn74910330) (`energy_distance_gersbach_comp/`), both on Weizhou's MuData
- [x] Huangfu HUES8 Definitive Endoderm → [`syn74883327`](https://www.synapse.org/Synapse:syn74883327) ⚠ p-value calibration concern; Sara comparison run → [`syn74910358`](https://www.synapse.org/Synapse:syn74910358)
- [x] Huangfu HUES8 Embryonic Stem Cell → [`syn74883475`](https://www.synapse.org/Synapse:syn74883475) ⚠ same calibration concern; Sara comparison run → [`syn74910472`](https://www.synapse.org/Synapse:syn74910472)
- [ ] Re-run Huangfu DE + ESC with HVG-subset PCA to fix anti-conservative p-values, then re-mirror
- [ ] Gersbach WTC11 Hepatocyte — bug **Sara** to deliver in our format ([`schemas/energy_distance.json`](schemas/energy_distance.json)). Sara likely has a complete run; we just need the deliverables shaped to match our schema and uploaded under `2026_UTSW/datasets/<id>/energy_distance/`.

**QC** (mirror script: in-place — `synapseclient` upload from `<run>/qc/`; intended_target + mapping_gene + mapping_guide subdirs)
- [x] Hon WTC11 Cardiomyocyte (our re-run on Weizhou's MuData) → [`syn74917453`](https://www.synapse.org/Synapse:syn74917453) (mirrored 2026-05-12)
- [x] Huangfu HUES8 Definitive Endoderm → [`syn74918479`](https://www.synapse.org/Synapse:syn74918479) (mirrored 2026-05-12)
- [x] Huangfu HUES8 Embryonic Stem Cell → [`syn74918600`](https://www.synapse.org/Synapse:syn74918600) (mirrored 2026-05-12)
- [x] Gersbach WTC11 Hepatocyte (our re-run on Sara's MuData) → [`syn74918946`](https://www.synapse.org/Synapse:syn74918946) (mirrored 2026-05-12)
- [ ] Engreitz WTC11 Endothelial — blocked: no data

**Calibration** (DEG empirical-null calibration — `scripts/run_calibration.sh`; SLURM submission per `.claude/skills/deg-calibration`)
- [x] Huangfu HUES8 Definitive Endoderm → [`syn74920615`](https://www.synapse.org/Synapse:syn74920615) (4 TSVs: all + cis + direct_target + trans; FDR<0.05 = 14,521; 2026-05-12)
- [x] Huangfu HUES8 Embryonic Stem Cell → [`syn74920616`](https://www.synapse.org/Synapse:syn74920616) (4 TSVs, 2026-05-12)
- [x] Hon WTC11 Cardiomyocyte (Weizhou's run) → [`syn74920617`](https://www.synapse.org/Synapse:syn74920617) (4 TSVs, 2026-05-12)
- [ ] Gersbach WTC11 Hepatocyte — calibration running (SLURM 10847496); upload pending
- [ ] Engreitz WTC11 Endothelial — blocked: no data

**Tracking**
- [x] Log Synapse paths in [`synapse_paths.tsv`](synapse_paths.tsv) as artifacts land (one row per dataset / `_reference_`, one column per output type)
- [x] Refreshed snapshot at [`state/2026_05_09_state.tsv`](state/2026_05_09_state.tsv) (2026-05-09; original `state/2026_05_07_state.tsv` kept as historical). Refresh again whenever a major status change lands.

## 7. Final documentation

- [ ] Create a Google Sheet capturing the same state (likely a friendlier view of `synapse_paths.tsv`).
- [ ] Link the jamboree planning doc to the Google Sheet.
- [x] Onboarding doc for jamboree participants → [`GETTING_STARTED.md`](GETTING_STARTED.md) (2026-05-09; refresh as outputs land).
- [ ] Make sure this repo is up to date end-to-end. Build skills (in addition to the GETTING_STARTED doc) so jamboree participants can explore the data and run analyses with agents.

# What's next (top 5)

Each item below has a detailed report under [`issues/`](issues/) — that's the artifact to read before having the conversation, and the doc to hand to the relevant collaborator.

In priority order:

1. **Bug Weizhou** (Hon team) for the rest of the Hon CM CRISPR outputs — partial mirror at [`syn74520421`](https://www.synapse.org/Synapse:syn74520421) is missing `pipeline_info/`. → [Issue 2](issues/hon-cm-crispr-bundle.md)
2. **Bug Sara** (Gersbach team) to deliver Gersbach Hep `crispr_pipeline/`, `cnmf/`, and `energy_distance/` in our schema-defined formats. → [Issue 3](issues/gersbach-hep-deliverables.md)
3. **Verify HTv2 cNMF testbed** (job 10577039) completes cleanly, then launch production cNMF on Huangfu DE + Huangfu ESC. → [Issue 5](issues/htv2-cnmf-testbed.md)
4. **Fix energy-distance p-value calibration** — re-preprocess with HVG-subset PCA in `src/tf_perturb_seq/crispr_pipeline/preprocess_mudata_local.py` and re-run on Huangfu DE + ESC. → [Issue 1](issues/edistance-calibration.md)
5. **Per-dataset READMEs** under `datasets/<production-dataset>/`. → [Issue 6](issues/per-dataset-readmes.md)

Plus the standing blocker: [Issue 4](issues/engreitz-no-data.md) — Engreitz endothelial data not on the IGVF portal yet.

Index of all open issues: [`issues/README.md`](issues/README.md).
