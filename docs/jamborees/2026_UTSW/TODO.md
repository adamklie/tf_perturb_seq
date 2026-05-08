# Overview

Prepare data and documentation for the 2026 UTSW jamboree. We are working with the 5 production datasets shown in `2026_05_07_state.png` (this directory).

# Process notes

- **Always double-check paths before using them.** Confirm whether a path is local, on GCS, on Synapse, or on the HPC.
- **A lot of the source data lives on the UCSD HPC.** Connection can be slow, so expect to organize things over SSH:
  - Host: `aklie@nrnb-login.ucsd.edu`
  - Project root: `/cellar/users/aklie/projects/tf_perturb_seq`
- When in doubt about which dataset, repo, or platform a task targets, ask before acting.
- **GitHub Project**: [TFP3](https://github.com/users/adamklie/projects/4) — milestones as draft items, issues as tasks
- **Synapse-as-we-go**: We upload each artifact to Synapse **directly from where it lives** (HPC / GCS / IGVF portal) as soon as it's ready, and log the path in `synapse_paths.tsv`. We do **not** stage everything locally first. This folder only holds docs + small simplified outputs.

# Steps

## 1. Set up the jamboree folder in this repo

Working dir: `tf_perturb_seq/docs/jamborees/2026_UTSW/`

- [ ] Build a TSV capturing the state shown in `2026_05_07_state.png` (one row per dataset, columns for each artifact/status).
- [ ] Start a `README.md` in this folder that we extend as decisions get made in the steps below.

## 2. Decide on local and cloud organization

- [ ] Pick a local layout. Candidate: one folder per dataset, with each analysis as a subdirectory inside.
- [ ] Decide on the matching layout on Synapse so the mirror is 1:1.

## 3. Define the machine-readable outputs

Outputs that computational collaborators and Claude instances will consume to run analyses and generate figures. Document each in the README.

**Reference data**
- [ ] IGVF GTF file
- [ ] TF metadata
- [ ] Experimental metadata
- [ ] Guide metadata

**Perturb-seq outputs**
- [ ] Inference MuData

**cNMF outputs**
- [ ] _TBD — fill in_

**Energy distance outputs**
- [ ] _TBD — fill in_

## 4. Define the human-readable outputs

Outputs aimed at general scientists — simpler artifacts for higher-level figures and exploration.

**Reference data**
- [ ] Simplified TF metadata
- [ ] Simplified experimental metadata
- [ ] Simplified guide metadata

**Perturb-seq outputs**
- [ ] _TBD — fill in_

**cNMF outputs**
- [ ] _TBD — fill in_

**Energy distance outputs**
- [ ] _TBD — fill in_

Storage considerations: _TBD — note size limits, where each lives (repo vs. Synapse vs. Drive)._

## 5. Stage simplified outputs locally

Only the small, human-readable artifacts (e.g., `*_simplified.tsv` in `reference/`) live in this repo. Bulky files do not get staged locally — they go straight to Synapse.

- [ ] Drop simplified reference tables (TF / experimental / guide metadata) into `reference/`.
- [ ] Add per-dataset READMEs under `datasets/<name>/` describing what's on Synapse and linking to it.
- [ ] Add small simplified summaries per analysis where they make sense.

## 6. Upload to Synapse as artifacts come in (interleaved with 3 & 4)

Done iteratively, not at the end.

- [ ] For each artifact defined in step 3 / 4: transfer directly from source (HPC / GCS / IGVF portal) to Synapse.
- [ ] Log the result in `synapse_paths.tsv`:
  - One row per dataset (or `_reference_` for cross-dataset refs)
  - One column per output type (clear, shared column names)
  - Each cell = the Synapse path
- [ ] Keep `2026_05_07_state.tsv` checkboxes in sync with what's actually on Synapse.

## 7. Final documentation

- [ ] Create a Google Sheet capturing the same state (likely a friendlier view of the TSV from step 6).
- [ ] Link the jamboree planning doc to the Google Sheet.
- [ ] Make sure this repo as a whole as up to date as possibe. We will make some skills and mds for the purposes of allowing folks in the jamboree to explore the data and run analyses with agents