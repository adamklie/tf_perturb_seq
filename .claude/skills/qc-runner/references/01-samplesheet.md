# samples.tsv schema + building it

The QC array driver reads a 3-column TSV with a mandatory header row. One row per `inference_mudata.h5mu` to QC.

## Schema

```tsv
input	outdir	run_name
```

| Column | Meaning | Constraints |
|---|---|---|
| `input` | Absolute path to `inference_mudata.h5mu` on HPC | Must exist; must be readable from carter-compute |
| `outdir` | Absolute path to the run's `qc/` directory | Will be created if missing; conventionally `<BASE_DIR>/<RUN_LABEL>/qc` |
| `run_name` | Prefix used for output filenames | Recommend `<DATASET_SHORT>_<RUN_LABEL>`; alnum + dashes + underscores |

The header is literally the first line; `qc_array.sh` does `LINE_NUM=$((SLURM_ARRAY_TASK_ID + 1))` so task 1 reads line 2.

## Example

```tsv
input	outdir	run_name
/cellar/users/aklie/projects/tf_perturb_seq/datasets/Hon_WTC11-benchmark_TF-Perturb-seq/cleanser_800/crispr_pipeline/pipeline_outputs/inference_mudata.h5mu	/cellar/users/aklie/projects/tf_perturb_seq/datasets/Hon_WTC11-benchmark_TF-Perturb-seq/cleanser_800/qc	Hon_WTC11-benchmark_cleanser_800
/cellar/users/aklie/projects/tf_perturb_seq/datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/seqspec_v3/crispr_pipeline/pipeline_outputs/inference_mudata.h5mu	/cellar/users/aklie/projects/tf_perturb_seq/datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/seqspec_v3/qc	Hon_WTC11-cardio_seqspec_v3
```

## Building with the helper

```bash
python3 .claude/skills/qc-runner/scripts/build_samples_tsv.py \
  --repo-root /cellar/users/aklie/projects/tf_perturb_seq \
  --output    /cellar/users/aklie/projects/tf_perturb_seq/scratch/qc_array_logs/samples_$(date +%Y_%m_%d).tsv
```

Flags:

```
--repo-root <path>            Repo root (the parent of datasets/)
--output <path.tsv>           Output TSV (header + N rows)
--datasets <ds1> <ds2> ...    Restrict to specific dataset dirs (basename match)
--runs <run1> <run2> ...      Restrict to specific run labels (basename match)
--require-mudata              (default) Skip rows where inference_mudata.h5mu is missing
--include-missing             Emit rows even if the mudata is missing (use with --dry-run later to find what's not staged)
--qc-subdir <name>            Override 'qc' as the outdir basename (default: qc)
```

By default the builder scans `datasets/*/<run>/crispr_pipeline/pipeline_outputs/inference_mudata.h5mu` and emits one row per existing file. The `run_name` defaults to `<dataset_short>_<run_label>`, where `dataset_short` strips the `_TF-Perturb-seq` suffix.

## Mirroring `inference_mudata.h5mu` from GCS to HPC

Stage 2 outputs land on GCS. To run QC on HPC, mirror first:

```bash
DATASET=<...>
DATA_DATE=<YYYY_MM_DD>
RUN_LABEL=<...>
HPC_REPO=/cellar/users/aklie/projects/tf_perturb_seq
DEST=${HPC_REPO}/datasets/${DATASET}/${RUN_LABEL}/crispr_pipeline/pipeline_outputs
mkdir -p ${DEST}
gsutil -m cp \
  gs://igvf-pertub-seq-pipeline-data/${DATASET}/${DATA_DATE}/outs/${RUN_LABEL}/pipeline_outputs/inference_mudata.h5mu \
  ${DEST}/
```

You only need `inference_mudata.h5mu` for QC; perturbo TSVs are not consumed by the QC modules. (They're needed for DEG calibration — see `deg-calibration` skill.)

For per-jamboree mirror conventions, also see [docs/jamborees/2026_UTSW/](../../../../docs/jamborees/2026_UTSW/).

## Common builder usage

### All datasets, single sweep label

```bash
python3 build_samples_tsv.py --repo-root <REPO> --output sw.tsv --runs cleanser_800
```

### One dataset, all runs that have mudatas

```bash
python3 build_samples_tsv.py --repo-root <REPO> --output hon.tsv --datasets Hon_WTC11-benchmark_TF-Perturb-seq
```

### Audit: what mudatas are missing locally?

```bash
python3 build_samples_tsv.py --repo-root <REPO> --output /tmp/audit.tsv --include-missing
# then diff against ls of crispr_pipeline/pipeline_outputs/ across all (dataset, run) dirs
```

## Validating the TSV before submission

```bash
awk -F'\t' 'NR>1 {if (NF != 3) print "BAD ROW " NR ": " $0; if (!-f $1) print "MISSING INPUT " NR ": " $1}' samples.tsv
```

Catches:
- Wrong column count (split error, embedded tabs in paths).
- `input` paths that don't exist on the current FS.

A common mistake is editing the TSV in a spreadsheet that converts tabs to multiple spaces. Always edit in a real text editor.
