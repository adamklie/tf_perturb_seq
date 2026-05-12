# CLI option reference

Two entry points: the bash wrapper at `scripts/run_calibration.sh` (frozen, repo-root) and the underlying Python tool at `src/tf_perturb_seq/inference/calibrate.py`. The wrapper is convenient but exposes only a subset of options.

## `scripts/run_calibration.sh`

```
--project-root <path>     (required) Repo root; used to find calibrate.py and activate .venv
--trans-results <path>    (required) perturbo_trans_per_element_output.tsv.gz
--mudata <path>           (required) inference_mudata.h5mu
--outdir <path>           (required) Output directory (created if missing)
--prefix <str>            (required) Prefix for output filenames
--null-method t-fit|ecdf  (optional) Default: t-fit
--cis-results <path>      (optional) PARSED BUT DROPPED — never passed to Python.
                                     Leaving this flag in place is harmless but useless.
--dry-run                 (optional) Echo the command and exit
```

What the wrapper does:

1. Validates all `--*` paths.
2. `source ${PROJECT_ROOT}/.venv/bin/activate`.
3. Builds and runs the Python invocation.

**Quirk:** the wrapper accepts `--cis-results` for forward compatibility but doesn't pass it. Cis results are picked up by `calibrate.py` from the MuData / inferred paths instead.

## `src/tf_perturb_seq/inference/calibrate.py`

```
--trans-results <path>            (required) per-element trans TSV
--mudata <path>                   (required) inference_mudata.h5mu
--outdir / -o <path>              (required) output dir
--prefix <str>                    (required) output filename prefix
--null-method t-fit|ecdf          Default: t-fit. See references/02-null-methods.md
--cis-window <int>                Default: 100000. bp window for cis annotation
--non-targeting-label <str>       Default: "non_targeting". Value in mudata.mod["guide"].var["type"]
                                   that identifies NTC elements
--guide-mod-key <str>             Default: "guide". MuData modality key for guide data
--one-sided                       Flag. If set, computes one-sided (right-tail) p-values.
                                   Default is two-sided (test for any-direction effect).
```

All flags can be combined freely.

## When to use the Python tool directly

The wrapper covers ~80% of cases. Call `calibrate.py` directly when you need any of:

- Non-default `--cis-window` (e.g. 500 kb, 1 Mb)
- One-sided test (`--one-sided`)
- Non-standard NTC label in `guide.var["type"]` (`--non-targeting-label`)
- Non-standard MuData modality key (`--guide-mod-key`)

Recipe:

```bash
cd <REPO_ROOT>
source .venv/bin/activate
python src/tf_perturb_seq/inference/calibrate.py \
  --trans-results <...> \
  --mudata <...> \
  --outdir <...> \
  --prefix <...> \
  --null-method t-fit \
  --cis-window 250000 \
  --non-targeting-label NTC \
  --guide-mod-key gRNA \
  --one-sided
```

## Common invocations (TFP3)

### Benchmark dataset, defaults

```bash
DATASET=Hon_WTC11-benchmark_TF-Perturb-seq
RUN_LABEL=cleanser_800
BASE_DIR=<REPO>/datasets/${DATASET}
bash <REPO>/scripts/run_calibration.sh \
  --project-root <REPO> \
  --trans-results ${BASE_DIR}/${RUN_LABEL}/crispr_pipeline/pipeline_outputs/perturbo_trans_per_element_output.tsv.gz \
  --mudata        ${BASE_DIR}/${RUN_LABEL}/crispr_pipeline/pipeline_outputs/inference_mudata.h5mu \
  --outdir        ${BASE_DIR}/${RUN_LABEL}/calibration \
  --prefix        ${DATASET}_${RUN_LABEL}
```

### Production dataset, larger cis window via Python

```bash
DATASET=Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq
RUN_LABEL=seqspec_v3
BASE_DIR=<REPO>/datasets/${DATASET}
source <REPO>/.venv/bin/activate
python <REPO>/src/tf_perturb_seq/inference/calibrate.py \
  --trans-results ${BASE_DIR}/${RUN_LABEL}/crispr_pipeline/pipeline_outputs/perturbo_trans_per_element_output.tsv.gz \
  --mudata        ${BASE_DIR}/${RUN_LABEL}/crispr_pipeline/pipeline_outputs/inference_mudata.h5mu \
  --outdir        ${BASE_DIR}/${RUN_LABEL}/calibration \
  --prefix        ${DATASET}_${RUN_LABEL}_cis500kb \
  --cis-window    500000
```

### Sweep both null methods

```bash
for m in t-fit ecdf; do
  bash <REPO>/scripts/run_calibration.sh \
    --project-root <REPO> \
    --trans-results ${BASE_DIR}/${RUN_LABEL}/crispr_pipeline/pipeline_outputs/perturbo_trans_per_element_output.tsv.gz \
    --mudata        ${BASE_DIR}/${RUN_LABEL}/crispr_pipeline/pipeline_outputs/inference_mudata.h5mu \
    --outdir        ${BASE_DIR}/${RUN_LABEL}/calibration \
    --prefix        ${DATASET}_${RUN_LABEL}_${m} \
    --null-method   ${m}
done
```

## Errors

- **`ERROR: Trans results not found`** — Stage 2 mirror is incomplete; `gsutil -m cp` the file locally first.
- **`ERROR: MuData not found`** — same; mirror `inference_mudata.h5mu`.
- **`ERROR: calibrate.py not found`** — `--project-root` is wrong; should be the repo root containing `src/tf_perturb_seq/`.
- **`source: no such file: .venv/bin/activate`** — venv missing; run `uv sync` at the repo root.
- **`KeyError: 'guide'`** during run — wrong MuData modality key; pass `--guide-mod-key <key>`.
- **`KeyError: 'non_targeting'`** during run — NTC label mismatch; check `md.mod["guide"].var["type"].value_counts()` and pass `--non-targeting-label <actual>`.
