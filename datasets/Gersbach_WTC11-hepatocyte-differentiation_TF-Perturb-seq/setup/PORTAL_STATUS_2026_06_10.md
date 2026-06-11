# Gersbach Hep — IGVF portal + GCS status (2026-06-10)

Checkpoint before re-running Stage 1 (portal → GCS sync). Tracking issue: [#28](https://github.com/adamklie/tf_perturb_seq/issues/28).

Portal queried 2026-06-10 via the IGVF API over the 47 GEX measurement sets
(`setup/samplesheets/hep_measurement_sets.txt`) and their 47 paired CRISPR
auxiliary sets (`setup/samplesheets/hep_ms_aux_pairs.tsv`), filtered to
`construct_library_sets.accession=IGVFDS3299AXST`. Counts are R1/R2 pairs
grouped by (sequencing_run, lane, flowcell_id, index), excluding
`deleted`/`revoked` files.

## Observations

### Read pairing
- **scRNA: 751 pairs** (one short of the expected 47 × 16 = 752).
- **gRNA: 752 pairs** (complete).
- **Exactly one unpaired read in the whole dataset**: `IGVFDS4761IGDX`
  (sub-pool `10XLane1-8_S10`), GEX measurement set, **run 1 / lane 3** —
  only **R2** is present:
  - R2 = `IGVFFI7067DKVR`, status `in progress`,
    `submitted_file_name = /work/rr151/TF_perturbSeqgenes_Helen/GEX/novaseq_1/B2_S10_L003_R2_001.fastq.gz`
  - No R1 mate exists for that run/lane in any status.
  - All other 15 lanes × 2 runs of S10 are paired; so are all 16 lanes of the
    other 46 sub-pools and all 47 gRNA aux sets.

### seqspecs
- **All 47 scRNA MS R1 files now carry seqspecs.**
- **10 of 47 gRNA aux sets carry seqspecs; 37 do not.**
- (Issue #28 recorded *zero* seqspecs across the dataset — these are being
  added progressively on the portal side. Non-blocking for us: the samplesheet
  uses fallback seqspec YAMLs, currently Hon CM's 10x v3 set.)

### onlist_files
- Still **null on all 47 measurement sets** (unchanged from issue #28 gap 3).
  Samplesheet uses fallback `IGVFFI9487JPEN` (`737K-august-2016.txt.gz`, 10x v2,
  matching Sara's `is_10x3v3=false` config).

### Resolved since issue #28
- Issue #28 gap 2 described S10 as missing a **gRNA** pair (16 scRNA / 15 gRNA).
  As of today S10 is **16 gRNA / 15 scRNA** — the gRNA side is now complete and
  the single remaining gap is on the scRNA side (the orphan R2 above).

## GCS state
- Bucket `gs://igvf-pertub-seq-pipeline-data/` has **no
  `Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/` prefix.** Nothing
  has been synced. (All other production datasets — Hon CM, Huangfu DE/ESC,
  benchmarks — are present.)
- No storage-transfer jobs exist since 2026-04-15 (`gcloud transfer operations
  list`), i.e. none for the June 3/4 attempts.
- The local `setup/samplesheets/sample_metadata_gcp_2026_06_0{3,4}.csv` contain
  `gs://` paths, but the referenced objects do not exist — the destination paths
  were written into the CSVs but the S3→GCS transfers were never executed/completed.
- `sample_metadata_gcp_2026_06_04.csv` matches today's portal state exactly
  (751 scRNA + 752 gRNA; S10 = 15 scRNA + 16 gRNA), so no re-query is needed —
  only the actual transfer (Step 2) + patch (Step 3).

## Ruled out
- **Not an access artifact**: the project service account
  (`adamklie@igvf-pertub-seq-pipeline.iam.gserviceaccount.com`) lists the bucket
  and every other production-dataset prefix; the Hep prefix is simply absent.
- **The missing scRNA R1 is not a status-filter artifact**: the orphan R2 is
  `in progress`, and no R1 exists for that run/lane under any status.

## Implication for the sync
- Re-running Step 2 today transfers all available files (751 scRNA + 752 gRNA
  pairs). S10 would run with 15/16 GEX lanes; the pipeline aggregates lanes, so
  the effect is one fewer lane of depth for that one sub-pool.
- When Ruhi uploads `B2_S10_L003_R1_001.fastq.gz`, that single pair can be
  re-synced and the samplesheet patched without redoing the rest.
