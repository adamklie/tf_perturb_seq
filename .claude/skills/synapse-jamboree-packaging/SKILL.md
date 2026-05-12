---
name: synapse-jamboree-packaging
description: Package TFP3 production datasets for the 2026 UTSW jamboree. Mirror canonical outputs (CRISPR pipeline, cNMF, energy distance, QC) from GCS or HPC to Synapse syn64423137/2026_UTSW/datasets/<dataset>/, conforming to JSON schemas in docs/jamborees/2026_UTSW/schemas/. Maintain the synapse_paths.tsv registry. Build working-group summary TSVs (wg1 QC / wg1 e-distance / wg1 TF-cross-lineage / wg3 / wg4 / wg5) for jamboree distribution. Triggers on keywords like jamboree, UTSW, Synapse, syn64423137, 2026_UTSW, mirror outputs, synapse_paths, working group, wg1, wg3, wg4, wg5, schema, packaging, upload to Synapse.
user_invocable: true
---

# Synapse Jamboree Packaging (2026 UTSW)

You are an interactive assistant for packaging TFP3 outputs for the 2026 UTSW jamboree. This skill mirrors per-dataset analysis bundles to Synapse, maintains the `synapse_paths.tsv` registry, and builds working-group summary TSVs.

## Constants

```
JAMB_ROOT:            docs/jamborees/2026_UTSW/
SYNAPSE_PARENT:       syn64423137
SYNAPSE_TARGET ROOT:  syn64423137/2026_UTSW/datasets/<dataset_id>/
LOCAL_REGISTRY:       docs/jamborees/2026_UTSW/synapse_paths.tsv
SCHEMAS:              docs/jamborees/2026_UTSW/schemas/{crispr_pipeline,cnmf,energy_distance,qc,...}.json
SCRIPTS:              docs/jamborees/2026_UTSW/scripts/
AUTH (Synapse):       SYNAPSE_AUTH_TOKEN in ~/.zshrc — invoke scripts via `zsh -ic '...'` to pick it up
AUTH (GCS):           gcloud auth as adamklie@igvf-pertub-seq-pipeline.iam.gserviceaccount.com
HPC STAGING:          /cellar/users/aklie/scratch/jamboree_pipeline_staging/<dataset>/   (NOT /tmp — only 20 GB on nrnb)
```

## What's packaged

For each of the **5 production datasets**:

| Dataset | Bundle types |
|---|---|
| Hon WTC11 Cardiomyocyte | CRISPR pipeline, cNMF, energy distance, QC |
| Huangfu HUES8 Definitive Endoderm | same |
| Huangfu HUES8 Embryonic Stem Cell | same |
| Gersbach WTC11 Hepatocyte | same |
| Engreitz WTC11 Endothelial | same |

Plus cross-dataset **reference data** (already uploaded): TF metadata, experimental metadata, IGVF GTF, guide library.

> _A 6th dataset (Gersbach HTv2 benchmark) is mirrored under `2026_UTSW/datasets/` as a small testbed only — **not** production, **not** enshrined in schemas. Per memory `[[htv2_benchmark_role]]`._

## Step 0: Identify action + dataset

Ask (or infer):

1. **What action?**
   - `mirror-crispr` — mirror CRISPR pipeline outputs from GCS → Synapse.
   - `mirror-cnmf` — mirror cNMF outputs from HPC → Synapse.
   - `mirror-edist` — mirror energy distance outputs from HPC → Synapse.
   - `mirror-qc` — mirror QC outputs from HPC → Synapse.
   - `build-wg1` / `build-wg3` / `build-wg4` / `build-wg5` — build working-group summary TSVs from Synapse-resident bundles.
   - `update-registry` — read current Synapse layout and update `synapse_paths.tsv`.
2. **Which dataset(s)?** One, several, or all 5 production.
3. **Dry-run first?** Most mirror scripts have `--dry-run` — recommended.

Then read the matching reference file.

| Topic | Reference file |
|---|---|
| Mirror recipes per bundle type | `references/01-mirror-recipes.md` |
| Schemas + canonical layouts | `references/02-schemas.md` |
| `synapse_paths.tsv` registry | `references/03-registry.md` |
| Working-group build scripts | `references/04-working-group-builds.md` |

## Step 1: Mirror outputs to Synapse

Each output type has its own mirror script. They all live in `docs/jamborees/2026_UTSW/scripts/` and follow the same pattern: stage from canonical source → upload to Synapse → record syn ID in `synapse_paths.tsv`.

```bash
# Authenticate (one-time per shell)
[[ -z "${SYNAPSE_AUTH_TOKEN:-}" ]] && { echo "SYNAPSE_AUTH_TOKEN not set; check ~/.zshrc"; exit 2; }
gcloud auth print-access-token >/dev/null || gcloud auth login

# CRISPR pipeline (GCS source)
zsh -ic 'python docs/jamborees/2026_UTSW/scripts/mirror_pipeline_outputs.py \
    --dataset <DATASET_ID> \
    --workdir /cellar/users/aklie/scratch/jamboree_pipeline_staging/<DATASET_ID> \
    --dry-run'   # remove --dry-run when ready

# cNMF (HPC source — run from nrnb)
zsh -ic 'python docs/jamborees/2026_UTSW/scripts/mirror_cnmf_outputs.py --dataset <DATASET_ID>'

# Energy distance (HPC source)
zsh -ic 'python docs/jamborees/2026_UTSW/scripts/mirror_edistance_outputs.py --dataset <DATASET_ID>'
```

For per-bundle invariants, recovery from partial uploads, and what to do when a dataset's source layout drifts from canonical, see `references/01-mirror-recipes.md`.

## Step 2: Validate against schemas

Each bundle has a JSON schema describing the expected directory tree and column contents. Spot-check after upload:

```bash
# Check Synapse contents
python -c "
import os, synapseclient
syn = synapseclient.Synapse(silent=True); syn.login(authToken=os.environ['SYNAPSE_AUTH_TOKEN'])
for f in syn.getChildren('<SYN_ID>'):
    print(f['name'], f['type'])
" 2>/dev/null

# Compare against the schema
cat docs/jamborees/2026_UTSW/schemas/crispr_pipeline.json | jq '.directories[].path'
```

A `validate_synapse_against_schema.py` script does not yet exist — eyeballing remains the workflow until a dataset trips on missing files.

For each bundle type's schema (top-level dirs, key files, columns), see `references/02-schemas.md`.

## Step 3: Update `synapse_paths.tsv`

Most mirror scripts auto-update the registry. To do it manually (e.g., after a partial mirror):

```bash
# View current state
column -t -s $'\t' docs/jamborees/2026_UTSW/synapse_paths.tsv

# Edit (vim / awk / pandas) to add the new (dataset, output_type) → syn_id row
```

The schema is documented in `references/03-registry.md`. After editing, commit — the registry is **tracked** (small, ground-truth for what's on Synapse).

## Step 4: Build working-group summaries

After at least 3 of 5 production datasets have canonical mirrors on Synapse, regenerate the working-group TSVs:

```bash
for wg in build_wg1_qc_summary build_wg1_edistance_summary build_wg1_significant_tfs \
          build_wg1_tf_cross_lineage build_wg1_trans_target_counts; do
  python docs/jamborees/2026_UTSW/scripts/${wg}.py
done

# WG3 / WG4 / WG5 — only after WG1 cross-dataset products exist
python docs/jamborees/2026_UTSW/scripts/build_wg3_tf_convergence_scorecard.py
python docs/jamborees/2026_UTSW/scripts/build_wg3_disease_tf_activity.py
python docs/jamborees/2026_UTSW/scripts/build_wg4_network_structure_by_lineage.py
python docs/jamborees/2026_UTSW/scripts/build_wg4_tf_gene_edges.py
python docs/jamborees/2026_UTSW/scripts/build_wg5_tf_family_scorecard.py
```

These read from `reference/*.tsv` (already-built cross-dataset summaries) and per-dataset Synapse mirrors. Outputs land under `working_groups/wg<N>_<topic>/`.

For per-script inputs, outputs, and dependencies, see `references/04-working-group-builds.md`.

## Status overview

Check the current packaging state any time:

```bash
column -t -s $'\t' docs/jamborees/2026_UTSW/synapse_paths.tsv
cat docs/jamborees/2026_UTSW/README.md | head -50    # 'Status at a glance' table
```

As of 2026-05-12 (per memory `[[2026_utsw_jamboree]]`):

- **CRISPR pipeline:** ✅ Hon CM (partial), Huangfu DE/ESC, Gersbach Hep (non-canonical); ☐ Engreitz Endothelial.
- **cNMF:** ✅ Huangfu DE/ESC; ⏳ Hon CM, Gersbach Hep, Engreitz.
- **Energy distance:** ✅ Huangfu DE/ESC (w/ calibration caveat), Hon CM, Gersbach Hep; ☐ Engreitz.
- **QC:** ✅ Huangfu DE/ESC, Hon CM; ⏳ Gersbach, Engreitz.

## Important notes

- **`zsh -ic` is critical** for Synapse-touching scripts. `SYNAPSE_AUTH_TOKEN` is in `~/.zshrc` (per memory `[[synapse_project]]`); a plain `bash python ...` won't pick it up.
- **HPC `/tmp` is 20 GB.** Pipeline bundles are ~60 GB each. Always pass `--workdir /cellar/users/aklie/scratch/jamboree_pipeline_staging/<id>` per memory `[[hpc_pipeline_mirror]]`.
- **`gcloud auth`** — use the IGVF service account, NOT your personal account. The IGVF service account has read access to `gs://igvf-pertub-seq-pipeline-data/`. Switch with `gcloud config configurations activate <name>` or `gcloud auth login adamklie@igvf-pertub-seq-pipeline.iam.gserviceaccount.com`.
- **cNMF on Synapse is FLAT** — `2026_UTSW/datasets/<DS>/cnmf/` (no `<cnmf_run_name>/` nesting). The local layout has `<DS>/<RUN>/cnmf/<cnmf_run_name>/` but the upload script flattens it. Inspect `mirror_cnmf_outputs.py` if confused.
- **Reference data is already uploaded** — TF metadata `syn74834227`, experimental `syn74834309`, GTF `syn74834518`, guide library `syn74834519`. Don't re-mirror unless updating.
- **The HTv2 testbed is NOT production** — see memory `[[htv2_benchmark_role]]`. Mirror to Synapse for staging if needed, but don't enshrine in schemas or working-group builds.
- **Schemas in `schemas/*.json`** are authoritative. If a mirror script disagrees, schemas win — file an issue against the script.
- **`synapse_paths.tsv` is the source of truth** for "what's on Synapse." If you edit Synapse manually (e.g., delete a stale folder), update this file in the same commit.
- For the full jamboree plan: [docs/jamborees/2026_UTSW/README.md](../../../docs/jamborees/2026_UTSW/README.md). For step-by-step prep: [TODO.md](../../../docs/jamborees/2026_UTSW/TODO.md).
