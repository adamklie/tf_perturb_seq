# Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq

## Dataset

Hon-lab production TF Perturb-seq dataset. WTC11 hiPSCs differentiated to cardiomyocytes, full TF library (~2000 targets) via CRISPRi, **HTO-multiplexed**. IGVF analysis set [IGVFDS6332VCTO](https://data.igvf.org/analysis-sets/IGVFDS6332VCTO). See `README.md` for portal status / lane counts and the open question about measurement-set counts.

## Run Pipeline on GCP

```bash
cd $PROJECT_ROOT

# Generate samplesheet
bash datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/1_generate_per_sample_metadata.sh

# Transfer to GCP -- dry run first
DRY_RUN=true bash datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/2_upload_to_gcp.sh
bash datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/2_upload_to_gcp.sh

# Patch .gz files (seqspec, barcode_onlist, guide_design)
bash datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/3_patch_gcp_files.sh

# Run pipeline on GCP
RUN_IN_BACKGROUND=true bash datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/4_run_CRISPR_pipeline.sh
```

## Sample Metadata

| Version | File | Notes |
|---------|------|-------|
| Initial | `sample_metadata_gcp_2026_04_15.csv` | Per-sample seqspecs from portal (gzipped) |
| Patched | `sample_metadata_gcp_2026_04_15_patched.csv` | Decompressed seqspec / barcode_onlist / guide_design files; 112 rows (28 MS × {1 scRNA + 2 gRNA + 1 hash}); GCS paths verified 2026-04-29 |
| Patched v2 | `sample_metadata_gcp_2026_04_15_patched_v2.csv` | Same as `_patched` except the `seqspec` column points at `patch_v2/` yamls with i7/i5 sample-index reads stripped (see issue noted in `seqspec_v2` run below) |

## Pipeline run reference

> **Naming convention.** Each entry below is keyed by our `RUN_LABEL` (set in `4_run_CRISPR_pipeline.sh`). `RUN_LABEL` drives the **output directory** (`.../outs/<RUN_LABEL>/`), the **config filename** (`<DATASET>_<RUN_LABEL>.config`), and the **log filename**. It is **not** the name you'll see on Seqera Tower — Nextflow auto-generates a separate run name per launch (e.g. `lonely_sammet`, `irreverent_liskov`). One `RUN_LABEL` can correspond to multiple Tower runs if we relaunch after fixes; all Tower names for a given `RUN_LABEL` are listed in its entry below.

### seqspec_v2 — 2026-04-29 (relaunch with stripped seqspecs) — _PARTIAL FIX, ABANDONED_

**Status:** stripping the i7/i5 reads fixed RNA mapping (no failures observed) and would also have unblocked guide mapping (mappingGuide carries an awk fallback). It did NOT fix the **hashing modality**, which failed with `Error: empty sequence list ... -- Missing output file(s) *_ks_hashing_out/counts_unfiltered/adata.h5ad`. The right path forward is to source/build proper seqspecs (see "Open issue" below); abandoning this run.

**Guide failure (also unfixed by v2): doubled cell barcode.** After v2 stripped the Index 1 read, mappingGuide's chemistry came out as `0,0,16,1,106,122:0,16,26,1,96,106:1,0,0`. kallisto bus accepts it (both file indices exist), but the 1st chemistry component now declares cell barcode at **two locations** — R1 pos 0–16 (16bp) AND R2 pos 106–122 (16bp). kallisto concatenates them → 32bp barcodes, then `bustools correct` rejects against the 16bp 10x whitelist:
```
Error: barcode length and on-list length differ, barcodes = 32, on-list = 16
```
Same root cause as the hashing failure: the production yaml's `library_spec` describes a Sigma 5'-Perturb backbone where R2 reads back through the structure and hits the cell-barcode region a second time. `seqspec index -t kb` faithfully encodes both barcode positions, but the kb count workflow expects a single barcode location matching the whitelist length. Lucas's curated seqspec sidesteps this by truncating the library_spec to describe only R1 for cell ID, even when the physical chemistry has R2-readthrough — that simplification is what makes his yamls work end-to-end with kb count, not just the i7 omission.

**Hashing failure root cause (also unfixed by v2):**
- For the `tag` modality, `seqspec index -t kb` against the Hon production yamls produces a chemistry with an **empty 3rd component**: `0,0,16,1,10,25,1,57,73:0,16,26,1,47,57:` ← nothing after the second `:`. The IGVF portal yamls' `library_spec` for the HTO region (a 203bp `joined` block) doesn't expose a clean tag-region position that `seqspec index` can translate into the 3rd `file,start,end` triple.
- `kallisto bus -x` rejects an empty 3rd component with `Error: empty sequence list`.
- `mappingGuide` avoids this in practice via an awk fallback that fills `1,0,0` when the spacer position is empty (loose "scan all of R2 by k-mer" instruction — works for kite workflow but lossy vs. a precise window).
- `mappingHashing/main.nf` does **not** apply that same fallback — it passes the chemistry to `kb count` verbatim, so the empty 3rd component fails immediately.

**Why the benchmark dataset (`Hon_WTC11-benchmark_TF-Perturb-seq`) doesn't hit this:** It uses Lucas's manually curated seqspecs at `gs://igvf-pertub-seq-pipeline-data/scratch/bioinfolucas/new_03_seqspec/`, e.g. `gary_IGVFDS4761PYUO_hash.yaml`. Those yamls have:
- only R1 + R2 in `sequence_spec` (no Index 1 read)
- a hand-built `library_spec` that explicitly defines `barcode → umi → hashing(15bp) → common(10bp) → r2_primer`, so `seqspec index` produces a clean chemistry: `0,0,16:0,16,26:1,10,25` — a precise 15bp window for the hash tag in R2.

Per-modality chemistry comparison (Lucas curated vs. Hon production / patch_v2):
| modality | Lucas (works) | Hon production v2 (fails) |
|---|---|---|
| rna | `0,0,16:0,16,26:1,0,92` | `0,0,16:0,16,26:1,0,98` (close enough; works) |
| crispr | `0,0,16:0,16,26:` (awk fills `1,0,0`) | `0,0,16,1,106,122:0,16,26,1,96,106:` (awk fills) |
| tag | `0,0,16:0,16,26:1,10,25` ← clean | `0,0,16,1,10,25,1,57,73:0,16,26,1,47,57:` ← empty 3rd, kallisto error |

**TODO — dig into the `library_spec` structure of production vs. benchmark.** Before settling on a fix, we need to understand *why* the IGVF portal HTO yaml's `library_spec` is shaped the way it is (203bp `joined` region with no exposed tag-region position) versus Lucas's hand-built one (explicit `barcode → umi → hashing(15bp) → common → r2_primer` decomposition). Same physical assay should produce the same chemistry — so either (a) the portal seqspec is incomplete/incorrectly auto-generated, (b) Lucas's decomposition is right but bespoke and we need it for every dataset, or (c) one of them describes a different version of the kit. Compare the full `library_spec` blocks side-by-side for a benchmark hash yaml and a production hash yaml; resolve which is correct; lift the answer back to the IGVF portal team or to Lucas. The same comparison is worth doing for `crispr` (the production yaml's R2 reads back through the structure → barcode/UMI appear at both R1 and R2 positions; Lucas's doesn't), to confirm whether the production assay genuinely has that backbone or if the seqspec is mis-decoded.

**Open issue — seqspec quality / pipeline robustness.** Two complementary fixes worth tracking:
1. _Data side:_ obtain or build properly-curated Hon-CM seqspecs (1 per modality, à la Lucas's). Once we have a clean hash yaml with a defined tag-region position, the production run will work without any pipeline patch.
2. _Pipeline side:_ port the `mappingGuide` awk fallback into `mappingHashing/main.nf` (and consider also `mappingscRNA`) so an empty 3rd component falls back to `1,0,0`. Wouldn't fix the underlying seqspec issue but would make the pipeline robust to it.

### seqspec_v2 — original launch entry (kept for record)

- **Config:** `Hon_WTC11-..._seqspec_v2.config` (copy of `_initial_run.config`, no param changes)
- **Sample metadata:** `sample_metadata_gcp_2026_04_15_patched_v2.csv` — `seqspec` column points at `gs://.../2026_04_15/patch_v2/`, where all 84 yamls have had i7/i5 sample-index `!Read` blocks stripped.
- **Output:** `gs://igvf-pertub-seq-pipeline-data/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/2026_04_15/outs/seqspec_v2`
- **Tower run names:** _TBD_
- **Why this run exists:** `initial_run` failed completely on mapping because the Hon-supplied seqspecs included an extra `!Read` block for the i7 sample index. `parsing_guide_metadata.py` feeds every read_id into `seqspec index -t kb`, so the chemistry came out referencing 3 fastq files per lane (R1 + i7 + R2). The pipeline only stages 2 (R1 + R2), so kallisto bus rejected with `Number of files (4) does not match number of input files required by technology ... (3)`. Errors were swallowed by `errorStrategy = ignore` and the workflow finished with `succeededCount=10; failedCount=84; ignoredCount=84` and zero mapping outputs in OUTDIR.
- **Fix applied:**
  - Stripping rule: drop any `!Read` block whose `primer_id` is `index5`/`index7` OR whose `max_len < 12` (catches an 8bp i7 mislabeled as `truseq_read2` in all 28 RNA seqspecs).
  - Local seqspec 0.3.1 verification: chemistries before vs. after stripping for one yaml of each modality (file index `2` → `1`):
    - crispr `0,0,16,2,106,122:0,16,26,2,96,106:` → `0,0,16,1,106,122:0,16,26,1,96,106:`
    - tag `0,0,16,2,10,25,2,57,73:0,16,26,2,47,57:` → `0,0,16,1,10,25,1,57,73:0,16,26,1,47,57:`
    - rna `0,0,16:0,16,26:2,0,98` → `0,0,16:0,16,26:1,0,98`
  - Stripped yamls uploaded to `gs://.../2026_04_15/patch_v2/` (84 files); v2 metadata CSV generated; `3_patch_gcp_files.sh` updated so future runs strip on decompression.

### initial_run — 2026-04-29 (first run)
- **Config:** `Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq_initial_run.config` (renamed from `..._2026_04_15.config` to match the `RUN_LABEL` convention).
- **Sample metadata:** `sample_metadata_gcp_2026_04_15_patched.csv`
- **Output:** `gs://igvf-pertub-seq-pipeline-data/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/2026_04_15/outs/initial_run`
- **Tower run names:** `gigantic_edison` (initial launch, currently running)
- **Config differences vs HUES8 run (`entertaining_hamster`):**
  - `ENABLE_DATA_HASHING = true` (HTO multiplexing — `barcode_hashtag_map` column populated for hash rows)
  - `is_10x3v3 = false` (different 10x chemistry)
  - `reverse_complement_guides = true`
  - `spacer_tag = "TAGCTCTTAAAC"` (Hon-specific)
- **Open questions to track:**
  - Measurement-set count discrepancy (28 in metadata vs. 32/26 in README) — verify against analysis set IGVFDS6332VCTO before interpreting outputs.
  - Same param-passing concern as HUES8 — confirm pipeline is using `-c ..._initial_run.config` rather than `CRISPR_Pipeline/nextflow.config` defaults.
