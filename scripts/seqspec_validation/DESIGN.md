# seqspec validation — setup step (DESIGN — DRAFT, not yet run)

Make seqspec validation a repeatable part of dataset setup. For every dataset under
`datasets/` that has seqspec YAML(s) + a samplesheet, run the `crispr_validator`
(Mode 1, `samplesheet`) against **one representative complete lane** (a unit that has
BOTH a scRNA and a gRNA row) to confirm the dataset's seqspec matches the actual
fastqs. Runs as a SLURM array on the NRNB `carter-compute` partition — **one array
task per dataset** — not locally (the guide exact-match scan is slow, so it must batch).

This doc + the three draft scripts in this directory are for review. **Nothing here has
been executed on the cluster and no jobs have been submitted.**

---

## 1. Inputs that already exist (verified)

### Dataset inventory (local repo, `datasets/`)

Only datasets with seqspec YAMLs under `setup/seqspec/` are validatable. Verified:

| Dataset | seqspec YAMLs (`setup/seqspec/`) | samplesheet | validatable? |
|---|---|---|---|
| `Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq` | `rna_seqspec.yml`, `guide_seqspec.yml` | `sample_metadata*.csv` (+ `_localseqspec_2026_06_10.csv`) | **yes** |
| `Hon_WTC11-benchmark_TF-Perturb-seq` | `rna_seqspec.yml`, `guide_seqspec.yml`, `hash_seqspec.yml` | `sample_metadata.csv` | **yes** |
| `Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq` | `rna_seqspec.yml`, `guide_seqspec.yml`, `hash_seqspec.yml` | `sample_metadata.csv` | **yes** |
| `Engreitz_WTC11-benchmark_TF-Perturb-seq` | none | yes | no (no seqspec) |
| `Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3` | none | yes | no (no seqspec) |
| `Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2` | none | yes | no (no seqspec) |
| `Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq` | none | yes | no (no seqspec) |
| `Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq` | none | yes | no (no seqspec) |
| `Huangfu_WTC11-benchmark_TF-Perturb-seq` | none | yes | no (no seqspec) |
| `technology-benchmark_WTC11_TF-Perturb-seq` | none | none | no |
| `ENCODE_WTC11_ChIP-seq` | none | none | no (not perturb-seq) |

The seqspec YAMLs under `setup/seqspec/` **are tracked in git**, so they are present at the
mapped repo path on NRNB (confirmed via `ls` on the cluster). The validator needs the seqspec
column to be an **absolute local path** — on NRNB that is
`/carter/users/aklie/projects/tf_perturb_seq/datasets/<ds>/setup/seqspec/{rna,guide}_seqspec.yml`.

### Samplesheet formats differ across datasets (important)

All sheets share the same 10 columns
(`R1_path,R2_path,file_modality,measurement_sets,sequencing_run,lane,seqspec,barcode_onlist,guide_design,barcode_hashtag_map`)
but the **cell values differ**:

| Dataset | R1/R2 | seqspec col | onlist / guide_design | lane | join key for scRNA↔gRNA |
|---|---|---|---|---|---|
| Hon benchmark | `gs://…` | `gs://…` | `gs://…` | empty | only one group (`run=1,lane=`); rows have different `measurement_sets` |
| Hon cardiomyocyte | `IGVFFI…` | `IGVFFI…` | `IGVFFI…` | populated, but **differs by modality** | **`measurement_sets`** (scRNA lane≠gRNA lane) |
| Hepatocyte (`_localseqspec_…`) | `IGVFFI…` | absolute local path | `IGVFFI…` | populated, same for both | `measurement_sets` (or `(run,lane)`) |

**Conclusion — the robust representative-lane key is `measurement_sets`.** In the
cardiomyocyte sheet scRNA and gRNA for the same biological sample sit on *different lanes*,
so grouping by `(sequencing_run, lane)` would NOT pair them. Grouping by `measurement_sets`
pairs them in every dataset (each `measurement_sets` accession carries one scRNA + one gRNA +
optionally one hash row). Hon-benchmark is the exception (rows differ in `measurement_sets`,
lane is blank) — there it falls back to "the single `(run,lane)` group" and just takes the
first scRNA + first gRNA row.

### The validator (verified working on NRNB)

`external/crispr_validator/seqspec_parser.py samplesheet` (gitignored; present on NRNB at
`/carter/users/aklie/projects/tf_perturb_seq/external/crispr_validator`). The working
invocation (reused verbatim from the local validation):

```
uv run --no-project --with pyyaml --with certifi --with seqspec python -u seqspec_parser.py samplesheet \
  --samplesheet <tiny_lane.csv> --analysis-root <out> --igvf-keypair <key.json> \
  --downloads-dir <cache> --barcode-sample-reads 3000 --feature-sample-reads 5000 --chunk-bytes 8000000
```

- `--group-by` default is `sequencing_run,lane`; we pass `--group-by measurement_sets` so the
  two rows in the tiny CSV form one group. (With a 2-row CSV the grouping is moot, but be
  explicit.)
- Output per group: `<analysis-root>/<group_label>/analysis_summary.json` (+ HTML reports).
  `analysis_summary.json["comparison_rows"]` is a list of `ComparisonRow` dicts — one per
  seqspec region — each with `modality`, `region`, `flag`
  (`perfect_match` / `close_enough` / `very_distant` / `missing_seqspec` / `missing_prediction`),
  `note`, `max_distance_bp`, `seqspec_*`/`prediction_*` interval/read/strand fields.

---

## 2. Per-dataset representative-lane selection logic

Implemented in `build_validation_samplesheets.py`. For each requested dataset:

1. Pick the **source samplesheet** (default: prefer `_localseqspec_*.csv`, else newest
   non-patched `sample_metadata*.csv`, else `sample_metadata.csv`; overridable per dataset).
2. Canonicalize `file_modality` (`scRNA sequencing`→`scRNA`, `gRNA sequencing`→`gRNA`, etc.).
3. Choose the representative unit:
   - Compute, per `measurement_sets` value, the set of modalities present.
   - Pick the **first** `measurement_sets` (sorted, deterministic) that has BOTH `scRNA` and
     `gRNA`. If none (e.g. Hon-benchmark with differing measurement_sets), fall back to: pick
     the first scRNA row and the first gRNA row regardless of measurement_sets and emit a
     warning that they may not be the same biological lane.
4. Emit a tiny 2-row CSV (`scRNA` + `gRNA`) with the same 10 columns, **rewriting**:
   - `seqspec` → the dataset's absolute seqspec path (rna for scRNA, guide for gRNA),
     expressed for the **target** (`--target nrnb|local`, default `nrnb`). NRNB base is
     `/carter/users/aklie/projects/tf_perturb_seq`.
   - `barcode_onlist` and `guide_design` → **prestaged absolute local paths** (see §3). This
     is the workaround for the validator's csv.gz/tsv.gz guide-design 404 bug.
   - `R1_path`/`R2_path` left as-is (IGVF accession or `gs://`) — the validator downloads a
     small chunk itself.
5. Write the CSV to `setup/seqspec_validation/<ds>__validation_lane.csv` under the dataset, and
   record a manifest row (dataset → csv path, chosen measurement_sets, seqspec paths,
   prestaged asset paths) used by the SLURM array driver.

## 3. Asset prestaging (onlist + guide_design) — the critical workaround

The validator infers `guide_design` IGVF accessions as `tabular-files/<acc>.tsv.gz`, but IGVF
guide-design files are frequently `file_format: csv` → `.csv.gz` → 404. Same risk for onlist.
**Workaround: pre-download both files and put absolute local paths in the CSV**, so the
validator treats them as local assets and skips its own download.

`build_validation_samplesheets.py` resolves each onlist/guide_design cell:
- **IGVF accession** (`IGVFFI…`): query the portal object
  (`https://api.data.igvf.org/tabular-files/<acc>/@@object?format=json`, HTTP basic auth from
  `IGVF_API_KEY`/`IGVF_SECRET_KEY` or `--keypair`), read its `href` + `file_format` to get the
  correct extension, then download `https://api.data.igvf.org<href>` to
  `setup/seqspec_validation/assets/<acc>.<ext>`.
- **`gs://…`**: `gsutil cp` to the assets dir.
- **absolute local path**: use as-is.

Prestaging writes assets under `datasets/<ds>/setup/seqspec_validation/assets/`. On NRNB the
CSV must reference the **NRNB** absolute asset path; so prestaging is intended to run **on
NRNB** (login node — gsutil + portal download are light) before the array submits. The build
script takes `--target nrnb|local` to decide which absolute prefix to bake into the CSV and is
designed to be runnable on the login node.

## 4. SLURM array layout (NRNB, carter-compute)

`submit_seqspec_validation.sh` — one array task per dataset.

- Driver reads a **manifest TSV** (`dataset, validation_csv, analysis_root, downloads_dir`,
  header row), one line per validatable dataset, produced by `build_validation_samplesheets.py
  --manifest <path>`.
- `SLURM_ARRAY_TASK_ID` → line `ID+1` (skip header), mirroring `scripts/qc_array.sh`.
- Headers match the repo convention (`-J`, `-c`, `--mem`, `--partition=carter-compute`, `-t`,
  `-o`/`-e` into a scratch log dir). NRNB needs `module load slurm/nrnb/23.02.7` to get
  `sbatch`/`sinfo` on PATH (login `.bashrc` does not `module` successfully); the submit
  instructions note this. `carter-compute` confirmed via `sinfo -s` (3 nodes, 14-day limit).
- Each task runs the verbatim `uv run … seqspec_parser.py samplesheet …` invocation against its
  dataset's tiny CSV, writing `analysis_summary.json` under
  `<analysis_root>=datasets/<ds>/setup/seqspec_validation/results/`.
- Modest resources (1–2 CPU, ~8G, ~30 min) — sampled reads are tiny.

## 5. Output location convention

```
datasets/<ds>/setup/seqspec_validation/
  <ds>__validation_lane.csv         # tiny 2-row representative CSV
  assets/<acc-or-name>.<ext>        # prestaged onlist + guide_design
  results/<group_label>/analysis_summary.json + *.html
```

Plus a run-level manifest + collected summary under
`scripts/seqspec_validation/runs/<date>/` (manifest.tsv, summary.tsv). `results/` is gitignored
by the existing `results/` rule, so only the CSV + summary TSV would be tracked if desired.

## 6. Pass/fail summary

`collect_results.py` walks every `datasets/*/setup/seqspec_validation/results/**/analysis_summary.json`
(or the manifest) and flattens `comparison_rows` into one TSV:

```
dataset  group_label  modality  region  flag  max_distance_bp  note
```

Plus a per-(dataset,modality) rollup: `worst_flag` (max over `FLAG_ORDER`:
perfect_match < close_enough < very_distant < missing_prediction < missing_seqspec) and a
boolean `pass` (= all regions perfect_match or close_enough). A dataset/modality with any
`very_distant`/`missing_*` is flagged for review.

Interpretation caveat (gotcha #3): the validator labels a region a "guide" only if its
`region_id`/`region_type` contains guide/grna/spacer/protospacer; a `region_type: cdna` guide
region is mislabeled rna and shows `missing_seqspec`. So a `missing_seqspec` on the gRNA
modality may be a labeling artifact, not a real mismatch — note in the summary, don't auto-fail.

## 7. Open questions

- **Hon-benchmark pairing**: rows have different `measurement_sets` and blank lane, so the
  scRNA + gRNA we pair are not provably the same biological lane. Acceptable for a seqspec
  structural check (the seqspec is per-modality, shared across lanes), but flag it.
- **NRNB checkout lag**: NRNB is a commit or two behind local — the hepatocyte
  `_localseqspec_2026_06_10.csv` is not on NRNB yet. Either `git pull` on NRNB first, or have
  the build script regenerate the lane CSV from `sample_metadata.csv` (which IS present). Build
  script defaults to regenerating, so it does not depend on the local-only sheet.
- **Credentials on NRNB**: confirm `IGVF_API_KEY`/`IGVF_SECRET_KEY` are exported (or a keypair
  JSON exists) on the login node for prestaging + R1/R2 chunk downloads.
- **hash modality**: out of scope (task says scRNA + gRNA). hash rows are dropped.
- Should this become a numbered `setup/scripts/5_validate_seqspec.sh` per dataset (matching the
  `1_…`–`4_…` convention) or stay centralized under `scripts/seqspec_validation/`? Drafted
  centralized; per-dataset wrapper is a thin follow-up if wanted.
