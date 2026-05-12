# Choosing `--null-method`

Two empirical-null methods are implemented. They give similar results when NTCs are plentiful; they diverge in the tails when NTCs are scarce.

## `t-fit` (default)

Fit a Student-t distribution (location fixed at 0) to the NTC z-values, then compute p-values from the fitted survival function.

**Use when:**

- You want smoother p-values in the tail (the t-fit can interpolate beyond the most extreme NTC observation).
- NTC count is modest (a few dozen to a few hundred elements) — t-fit extracts more information per NTC by assuming a parametric shape.
- You're comparing across runs where NTC counts differ — the parametric fit reduces sample-size sensitivity.

**Watch out for:**

- Heavy-tailed or skewed NTC distributions (e.g. one extreme outlier NTC) can pull the fit. Inspect the histogram of NTC z-values; if it's clearly bimodal or has a heavy positive tail, switch to `ecdf`.
- Default winsorization is `0.01` (1% at each tail) before fitting. If the fit looks bad, override (only callable from Python, not via the bash wrapper).

## `ecdf`

Rank-based empirical CDF. For each test, count the fraction of NTC statistics at least as extreme.

```
p = (r + 1) / (B + 1)
   where r = # NTC stats >= |test|, B = # NTC stats
```

The `+1`'s are a bias correction so that `p` can never be exactly 0.

**Use when:**

- NTCs are plentiful (hundreds+; production datasets with full TF library).
- You want fully distribution-free p-values (no parametric assumptions).
- Reviewers want the most "obviously correct" method.

**Watch out for:**

- **Floor on `p`**: minimum `p` is `1 / (B + 1)`. With 30 NTCs that's `~0.032` — you cannot get a p-value smaller than 0.032 no matter how strong the effect. For the TFP3 benchmark (~30 NTCs in a 50-gene library), this is the binding constraint and `t-fit` is preferred.
- "Step" pattern in the p-value histogram from coarse discretization.

## How to choose for TFP3

| Dataset class | NTC count | Recommendation |
|---|---|---|
| Benchmark (50-gene library, WTC11) | ~30 NTC elements | `t-fit` (default) — ecdf's p-value floor is too coarse |
| Production (full TF library, ~2000 elements) | Hundreds of NTC elements | Either; `ecdf` is more conservative and more defensible |
| Pilot / very small experiment | <20 NTC | `t-fit` mandatory; or consider whether calibration is even valid |

## Diagnostic: p-value histogram of NTC vs targeting

A common sanity check after calibration:

```python
import pandas as pd, matplotlib.pyplot as plt
df = pd.read_csv("<prefix>_calibrated_trans_results.tsv", sep="\t")
# NTCs aren't in this table (they're the null); plot the targeting p's
df["empirical_pvalue"].hist(bins=50)
plt.title("Calibrated p-values (targeting)")
```

**Expectations:**

- Under a well-calibrated null with no signal: uniform (flat) histogram.
- With signal: spike near 0, flat elsewhere.
- **Bad: spike near 0 *and* near 1, or strong concavity** — null is misspecified. Try the other method, or check NTC identification.
- **Bad: spike only near 1** — null is broader than the targeting distribution (over-conservative). Often a sign of NTC contamination by guides that have real effects.

## Switching methods mid-comparison

If you re-calibrate the same run with a different `--null-method`, **bump the `--prefix`** so the old outputs aren't overwritten:

```bash
bash scripts/run_calibration.sh ... --prefix <dataset>_<run>_tfit ...
bash scripts/run_calibration.sh ... --prefix <dataset>_<run>_ecdf --null-method ecdf ...
```

Then diff the two prefixes' `_calibrated_trans_results.tsv` to see where the methods agree/disagree:

```bash
join -t$'\t' -1 1 -2 1 \
  <(sort -t$'\t' -k1 <run>_tfit_calibrated_trans_results.tsv) \
  <(sort -t$'\t' -k1 <run>_ecdf_calibrated_trans_results.tsv) \
  | awk -F'\t' '{print $1, $X_tfit_p, $Y_ecdf_p}'  # X, Y = column indices
```
