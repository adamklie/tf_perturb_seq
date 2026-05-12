# Outputs + local mirror

## GCS layout

After a completed run:

```
gs://igvf-pertub-seq-pipeline-data/<DATASET>/<DATA_DATE>/outs/<RUN_LABEL>/
├── pipeline_outputs/             # *** PRIMARY ***
│   ├── inference_mudata.h5mu     # final MuData (Stage 3 + Stage 5 consume this)
│   ├── perturbo_cis_per_element.tsv
│   ├── perturbo_cis_per_guide.tsv
│   ├── perturbo_trans_per_element.tsv
│   └── perturbo_trans_per_guide.tsv
├── pipeline_dashboard/           # *** PRIMARY ***
│   ├── dashboard.html            # human QC view
│   ├── inference_mudata.h5mu     # mirror of pipeline_outputs
│   └── additional_qc/            # gene / guide / intended_target / trans metrics
├── pipeline_info/                # *** TRACKED LOCALLY ***
│   ├── params_<timestamp>.json   # actual config used (RUN PROVENANCE)
│   └── versions.yml              # tool versions
├── mergedresults/                # final merged inference outputs (alt path to pipeline_outputs)
├── sceptre/                      # sceptre inference per-element/per-guide + chunk_manifest.tsv
├── inference/                    # additional inference outputs
├── evaluation/                   # built-in evaluation_output/, plots/
├── tf/benchmark_output/          # TF enrichment vs ChIP-seq
├── createmudata/                 # pre-inference combined MuData
├── mudata/                       # concat_mudata.h5mu (alternative)
├── preprocessanndata/            # filtered_anndata.h5ad + figures
├── anndata/                      # concatenated pre-QC anndata
├── filter/                       # per-MS hashing filtered h5ads
├── prepare/                      # per-MS pre-guide-assignment mudatas
├── guide/                        # per-MS post-guide-assignment mudatas
├── demultiplex/                  # per-MS HTO demux (HTO datasets only)
├── hashing/                      # concatenated HTO demux (HTO datasets only)
├── mappingscrna/<MS>_ks_transcripts_out/   # kallisto bus per measurement set
│   └── run_info.json             # n_processed, n_pseudoaligned, p_pseudoaligned
├── mappingguide/<MS>_ks_guide_out/         # per-MS guide alignment
├── mappinghashing/<MS>_ks_hashing_out/     # per-MS HTO alignment (HTO datasets only)
├── seqspecparser/                # parsed seqspecs
├── seqspeccheck/                 # seqspec position tables + plots
├── createguideref/               # guide index + mismatch ref
├── createhashingref/             # HTO index (HTO datasets only)
├── downloadreference/            # transcriptome .idx + t2g
└── downloadgtf/                  # GENCODE GTF
```

For per-stage detail, see [docs/analysis/crispr_pipeline/CRISPR_PIPELINE_OUTPUTS.md](../../../../docs/analysis/crispr_pipeline/CRISPR_PIPELINE_OUTPUTS.md).

## What to mirror locally

Mirror only the small, durable artifacts. Bulk outputs (h5mu, h5ad, kallisto matrices) stay on GCS — they're regenerable.

```
<BASE_DIR>/<RUN_LABEL>/crispr_pipeline/
├── pipeline_info/             # TRACKED — params + versions, ~50 KB
│   ├── params_*.json
│   └── versions.yml
├── pipeline_outputs/          # gitignored — bulk
│   ├── inference_mudata.h5mu
│   └── perturbo_*.tsv
├── pipeline_dashboard/        # gitignored — bulk + HTML
│   ├── dashboard.html
│   └── additional_qc/
└── anndata/                   # gitignored — h5mu/h5ad
```

`pipeline_info/` is the **only** tracked directory. Everything else is gitignored per `.gitignore`'s `crispr_pipeline/` rules — see `docs/data/DATA.md` §"Per-dataset layout" for the exact rules.

## Mirror recipe

```bash
DATASET=<...>
DATA_DATE=<YYYY_MM_DD>
RUN_LABEL=<...>
BASE_DIR=/path/to/repo/datasets/${DATASET}
GCS=gs://igvf-pertub-seq-pipeline-data/${DATASET}/${DATA_DATE}/outs/${RUN_LABEL}
LOCAL=${BASE_DIR}/${RUN_LABEL}/crispr_pipeline

# 1. pipeline_info (TRACKED — keep small files only)
mkdir -p ${LOCAL}/pipeline_info
gsutil -m cp -r ${GCS}/pipeline_info/* ${LOCAL}/pipeline_info/

# 2. pipeline_outputs (gitignored, but needed locally for Stage 3+)
mkdir -p ${LOCAL}/pipeline_outputs
gsutil -m cp ${GCS}/pipeline_outputs/inference_mudata.h5mu ${LOCAL}/pipeline_outputs/
gsutil -m cp ${GCS}/pipeline_outputs/perturbo_*.tsv ${LOCAL}/pipeline_outputs/

# 3. pipeline_dashboard (optional — useful for slides)
mkdir -p ${LOCAL}/pipeline_dashboard
gsutil -m cp -r ${GCS}/pipeline_dashboard/* ${LOCAL}/pipeline_dashboard/

# 4. anndata (optional — pre/post-QC anndata if needed for ad-hoc exploration)
# mkdir -p ${LOCAL}/anndata
# gsutil -m cp -r ${GCS}/preprocessanndata/*.h5ad ${LOCAL}/anndata/
```

Then commit only the `pipeline_info/` mirror:

```bash
cd <repo>
git add datasets/${DATASET}/${RUN_LABEL}/crispr_pipeline/pipeline_info/
git status   # verify nothing else is staged
```

There's a helper for mirroring (`scripts/sync_gcp_run.sh` at the repo root) — check there first before writing a one-off.

## `params_*.json` — the run-provenance file

The single most important file in the output. It records the **actual** values Nextflow used after profile/config merging, including all defaults you didn't set. Use it to:

- Verify two runs differ only by the param you intended to change.
- Reproduce a run months later (everything in `params {}` is here).
- Audit which container tag was used (`containers.base`, `containers.sceptre`, etc.).

Quick comparisons:

```bash
diff <(jq -S . runA/pipeline_info/params_*.json) <(jq -S . runB/pipeline_info/params_*.json)
```

Folder names alone are misleading. **Always read the params JSON to know what differs between two runs.**

## Handoff to Stage 3 (QC)

Stage 3 consumes `pipeline_outputs/inference_mudata.h5mu`. After the local mirror is in place:

```bash
ls ${LOCAL}/pipeline_outputs/inference_mudata.h5mu
```

Then invoke the QC skill (when built): `/qc-runner`. It expects `<RUN_LABEL>/crispr_pipeline/pipeline_outputs/inference_mudata.h5mu` and writes to `<RUN_LABEL>/qc/`.

## Common post-run checks

```bash
# Did all measurement sets map?
for d in $(gsutil ls ${GCS}/mappingscrna/); do
    name=$(basename "$d")
    p=$(gsutil cat "${d}run_info.json" 2>/dev/null | jq -r '.p_pseudoaligned')
    echo "$name  p_pseudoaligned=$p"
done

# Per-cell counts after QC?
gsutil ls ${GCS}/preprocessanndata/

# Did inference complete for both methods?
gsutil ls ${GCS}/mergedresults/inference_mudata.h5mu
```

Anything anomalous → re-check `dashboard.html` and Tower's resources tab; consider whether a re-run with adjusted QC params is warranted.
