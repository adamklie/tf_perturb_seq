# Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq

## Dataset

Production TF Perturb-seq dataset from the Huangfu lab. HUES8 hESCs differentiated to definitive endoderm, with CRISPRi perturbation of the full TF library (~2000 targets), pools ABCD. IGVF analysis set: [IGVFDS9951KTRR](https://data.igvf.org/analysis-sets/IGVFDS9951KTRR).

## Run Pipeline on GCP

```bash
cd $PROJECT_ROOT

# Generate samplesheet
bash datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/1_generate_per_sample_metadata.sh

# Transfer to GCP -- dry run first
DRY_RUN=true bash datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/2_upload_to_gcp.sh
bash datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/2_upload_to_gcp.sh

# Patch .tsv.gz / seqspec files to unzipped versions
bash datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/3_patch_gcp_files.sh

# Run pipeline on GCP
RUN_IN_BACKGROUND=true bash datasets/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/4_run_CRISPR_pipeline.sh
```

## Sample Metadata

| Version | File | Notes |
|---------|------|-------|
| Original | `sample_metadata_gcp_2026_04_09.csv` | Per-sample seqspecs from portal (gzipped) |
| Patched | `sample_metadata_gcp_2026_04_09_patched.csv` | Decompressed seqspec / barcode_onlist / guide_design files |

## Pipeline run reference

> **Naming convention — read this first.** Each entry below is keyed by our `RUN_LABEL` (set in `4_run_CRISPR_pipeline.sh`). `RUN_LABEL` drives the **output directory** (`.../outs/<RUN_LABEL>/`), the **config filename** (`<DATASET>_<RUN_LABEL>.config`), and the **log filename**. It is **not** the name you'll see on Seqera Tower — Nextflow auto-generates a separate run name per launch (e.g. `lonely_sammet`, `irreverent_liskov`). One `RUN_LABEL` can correspond to multiple Tower runs if we relaunch after fixes; all Tower names for a given `RUN_LABEL` are listed in its entry below.


### muddy_penguin — 2026-04-14 (spacer_tag +G)
- **Config:** `Huangfu_HUES8-..._muddy_penguin.config` — duplicated from `entertaining_hamster` with one change: `spacer_tag = "GAGTACATGGGGG"` (extra leading G, 5 G's instead of 4).
- **Sample metadata:** `sample_metadata_gcp_2026_04_09_patched.csv` (unchanged)
- **Output:** `gs://igvf-pertub-seq-pipeline-data/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/2026_04_09/outs/muddy_penguin`
- **Tower run names:** `distraught_borg` (initial launch), `friendly_morse` (resumed run — currently running). `-resume` was re-added for this relaunch to pick up from the prior cache.
- **Rationale (per Gary Yang, Slack 7:34 AM):** production data may carry an additional G in front of the guide sequences vs. the benchmark spacer. This run tests whether that extra G changes guide-assignment accuracy. Gary wants to validate accuracy on ES and DE before patching measurement/auxiliary sets with the corresponding seqspec files.

### entertaining_hamster — 2026-04-13 (initial run)
- **Tower run names:** `lonely_sammet` (first launch, failed on GCS auth), `irreverent_liskov` (relaunch after seqSpecCheck fix)

- **Config:** `Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq_entertaining_hamster.config`
- **Sample metadata:** `sample_metadata_gcp_2026_04_09_patched.csv`
- **Output:** `gs://igvf-pertub-seq-pipeline-data/Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/2026_04_09/outs/entertaining_hamster`
- First HUES8 definitive-endoderm run through the standardized CRISPR FG pipeline. Config adapted from the WTC11 benchmark (same 10x 3' v3 chemistry, SCEPTRE guide assignment + inference, MOI=high). Metadata CSVs keep their `2026_04_09` generation date; the run label is `entertaining_hamster`.

**Issues encountered:**
- _First launch (2026-04-13):_ failed at `seqSpecCheck` with `ValueError: No guide column found in metadata`. Root cause: the IGVF portal delivered the guide_design file as `IGVFFI8270UPKB.csv.gz` but the contents are tab-separated. `CRISPR_Pipeline/bin/seqSpecCheck.py:380` picks its pandas delimiter from the file extension (`,` for `.csv`, `\t` otherwise), so reading it as CSV collapsed everything into one column and no `spacer`/`guide` column was found. Fix: copied the GCS file to `IGVFFI8270UPKB.tsv`, updated every row of `sample_metadata_gcp_2026_04_09_patched.csv` to reference the `.tsv`, and patched `3_patch_gcp_files.sh` so future decompression writes `.tsv` directly.

**TODO — verify parameter passing:** Suspect the pipeline may be using `CRISPR_Pipeline/nextflow.config` defaults rather than (or in addition to) the `-c ..._entertaining_hamster.config` we pass. Need to confirm which params actually take effect (e.g. `is_10x3v3`, `spacer_tag`, `QC_pct_mito`, `INFERENCE_SCEPTRE_control_group`) — diff the two param blocks and inspect Nextflow's resolved params in the Tower run / `.nextflow.log`.
