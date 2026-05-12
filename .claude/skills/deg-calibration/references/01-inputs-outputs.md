# Calibration inputs & outputs

## Inputs

### `perturbo_trans_per_element_output.tsv.gz`

PerTurbo's per-element trans-test TSV from Stage 2. Each row is one (element × gene) test.

Expected columns (subset; check actual header — pipeline evolves):

| Column | Meaning |
|---|---|
| `intended_target_name` | Element ID (gene symbol or NTC name). NTCs are identified by `--non-targeting-label`. |
| `gene_id` / `gene_name` | Tested gene |
| `log2_fc` | log2 fold change of gene expression in perturbed vs control cells |
| `log2_fc_std` | Standard error of `log2_fc` |
| `p_value` | PerTurbo posterior p-value (NOT calibrated) |
| `gene_chr` / `gene_start` / `gene_end` | For cis annotation |

### `perturbo_cis_per_element_output.tsv.gz`

Same schema, but restricted to gene–element pairs within the cis window upstream of Stage 2.

### `inference_mudata.h5mu`

MuData with at minimum:

- `mudata.mod["guide"].var` — guide metadata, including `type` column with NTC label.
- Gene metadata for cis-window annotation (chromosome, start, end).

## Outputs

All written to `<outdir>/<prefix>_*.tsv` (TSV, not gzipped — for easy `cut`/`awk`).

### `<prefix>_calibrated_trans_results.tsv`

Master table — every (targeting element × gene) test from the trans input, with calibration added.

| Column | Added by | Meaning |
|---|---|---|
| (all PerTurbo input columns) | upstream | passed through |
| `z_value` | calibration | `log2_fc / log2_fc_std`, computed before calibration |
| `empirical_pvalue` | calibration | Calibrated p-value from `--null-method` |
| `fdr` | calibration | BH-corrected q-value across the targeting tests only |
| `cis_distance` | calibration | bp from element midpoint to gene TSS (NaN for trans-only) |
| `is_cis` | calibration | bool: within `--cis-window` |
| `is_direct_target` | calibration | bool: gene matches the perturbed element's `intended_target_name` |

### `<prefix>_calibrated_direct_target_results.tsv`

Subset where `is_direct_target == True`. One row per element (TF) tested against its own target gene. Useful for asking "how well did each TF KD knock down its own mRNA?"

### `<prefix>_calibrated_cis_results.tsv`

Subset where `is_cis == True`. Cis-regulatory effects within `--cis-window`.

### `<prefix>_calibrated_trans_only_results.tsv`

Subset where `is_cis == False`. Trans effects (gene on a different chromosome or outside the cis window).

## What "calibration" actually does

For every (targeting element × gene) row:

1. Compute `z_value = log2_fc / log2_fc_std`.
2. Build the **empirical null distribution** from the z-values of NTC element rows in the same gene context.
3. Compute `empirical_pvalue` by comparing the targeting z to the null:
   - `ecdf`: rank-based, `p = (r + 1) / (B + 1)` where `r` = # null abs(z) ≥ abs(test) and `B` = # null z-values.
   - `t-fit`: fit a Student-t to the null abs(z), evaluate survival function at the test |z|.
4. After all tests are computed, apply Benjamini–Hochberg across the **targeting** discovery set.

NTC rows themselves get `empirical_pvalue = NaN` and `fdr = NaN` — they're the null reference, not part of the discovery set.

## Verifying NTC identification

```python
import mudata
md = mudata.read_h5mu("inference_mudata.h5mu")
print(md.mod["guide"].var["type"].value_counts())
# Expect: non_targeting  ...
#         targeting      ...
```

If the NTC label isn't `non_targeting`, pass `--non-targeting-label <label>` to `calibrate.py` directly (the bash wrapper doesn't expose this flag).

## Common shape sanity checks

```bash
# Row counts roughly match
zcat perturbo_trans_per_element_output.tsv.gz | wc -l    # ~ N rows
wc -l <prefix>_calibrated_trans_results.tsv               # ~ N - (#NTC tests)

# Most direct targets should be significant for a working TF library
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="fdr")c=i; next} c && $c<0.1' \
    <prefix>_calibrated_direct_target_results.tsv | wc -l
```

In a well-behaved TFP3 dataset (Hon WTC11 benchmark, etc.) the majority of TFs should show their own gene as a significant direct target at q<0.1. If <30% do, suspect a guide-assignment or labeling issue upstream.
