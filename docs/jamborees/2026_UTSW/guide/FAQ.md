# FAQ

Common questions from people who've just opened a result file.

---

### Why are some of my p-values exactly zero?

A few possibilities:

1. **Permutation p-values bottom out at the resolution of the permutations.** If you ran 1,000 permutations and none beat the observed value, the smallest reportable p is ~ 1/1000 — often reported as 0 in compressed output.
2. **The test is anti-conservative.** If even the negative controls come out with p = 0, the null distribution is too tight relative to the real data. See the [calibration caveat](ENERGY_DISTANCE.md#the-calibration-caveat-read-this) for the energy-distance version of this. The fix is usually re-running with a different PCA basis (HVG-subset instead of all genes).
3. **The signal really is enormous.** For a positive control like *AARS* in a healthy run, p = 0 is the right answer.

If the negative controls also have p = 0, **don't gate on p alone**; use an effect-size cutoff (`distance_mean > NC max`).

---

### What's the difference between cis and trans?

- **Cis effect** = the knockdown changed *the gene it was supposed to knock down*. This is the QC: "did the perturbation work?"
- **Trans effect** = the knockdown changed *other* genes (downstream targets). This is the biology — the discovery layer.

Same data, different rows in the perturbo TSVs. Cis is the smaller table (one row per TF × its intended target); trans is the genome-wide one.

---

### How many TFs should be significant in a healthy dataset?

Depends on cell type and library, but rough ballparks for a ~2,000-TF library:

- **Cis-significant** (knockdown actually worked): 60–90% of guides should show ≥ 50% knockdown on their target. Lower than 50% = mediocre run; lower than 30% = a flag.
- **Trans-significant** (gave any signal at all): 50–200 TFs at energy-distance `distance_mean > NC max` is normal for a mature dataset. Big variance by lineage — cardiomyocytes are denser than ESCs in our data.

If you're seeing 1,000+ TFs called significant, the test is probably mis-calibrated. If you're seeing 0–10, check whether positive controls came out — if they did, you might just have a quiet cell type; if they didn't, something's wrong.

---

### Why do the absolute energy-distance numbers differ so much across datasets?

Because the distance metric depends on the PCA basis the cells were embedded in, which depends on the cells' actual variance structure. The Huangfu ESC distance scale is ~100× the Hon CM scale not because the perturbations are 100× stronger — but because ESC has more between-cell variance to begin with.

**Always interpret distances within a dataset, not across datasets.** For cross-dataset comparisons, use ranks or the "distance > NC max" boolean — never raw values.

---

### What's the difference between "negative control" and "non-targeting control" (NTC)?

- **NTC** = a guide whose sequence matches nowhere in the genome. The cells get the CRISPR machinery and no perturbation. These are the natural null distribution.
- **Negative control** = a guide that targets a benign region (a safe-harbor locus, an intergenic region). Those cells *do* have CRISPRi activity, but it shouldn't affect the cell state. Useful for catching off-target / general-CRISPR effects beyond what NTCs control for.

Both are present in our library. NTCs do the heavy lifting in p-value calibration; negative controls are extra null cells that should *also* not produce strong effects.

---

### I see "n_outlier_grnas_in_run" — what's an outlier guide?

The energy-distance pipeline flags individual guides as outliers if they disagree with the other guides targeting the same TF. Concretely:

- For a TF with 5 guides, we look at the energy distance each guide induces.
- If 4 of 5 land near each other and 1 is way off, that one is an outlier.
- The pipeline drops outlier guides before computing the per-TF distance, so the headline number reflects the consensus of the consistent guides.

`n_outlier_grnas_in_run` tells you how many of a TF's guides got dropped. If it's high (e.g. 3 of 5 guides flagged), be cautious — only 2 guides supported the call.

---

### What's HT-like chemistry?

In the 10x Genomics ecosystem, "HT" refers to the **high-throughput** versions of their kits (10x 3' HT, 10x 5' HT). They process more cells per lane than the standard kits, with slightly different chemistry / capture rates. Some of our datasets use HT-like chemistry; you'll see this in `experimental_metadata.tsv`. Functionally, it doesn't change how you read the results — it shifts cell counts up.

---

### What's a "knee" or "knee plot" and why does it matter?

In single-cell data, plotting log(cell rank) vs log(total UMIs per cell barcode) gives an L-shaped curve. The "knee" is the inflection between "real cells" (high UMI counts) and "empty droplets" (low). The pipeline's filter step uses the knee to draw a cell-calling threshold. A bad knee (no clear inflection, multiple inflections) means the cell-calling is unreliable — flag it.

---

### Why are there 25,312 genes in the gene-programs file but only 2,000 in the gene-universe file?

cNMF runs in two phases:
1. **Fit** on the top 2,000 HVGs (the "gene universe").
2. **Re-project** the program loadings onto *all* expressed genes — typically ~25,000 — to give every gene a loading score per program.

The 2,000-gene set is what cNMF saw during training; the 25,000-gene set is what cNMF *predicts* loadings for. The DACC spec calls the first one "gene universe" and the second one "gene programs" — the naming is confusing but technically correct.

---

### My dataset has `pipeline_outputs/` but no `pipeline_info/`. Should I worry?

Yes — `pipeline_info/` carries the Nextflow params + software versions that say *which version of the pipeline* produced these outputs. Without it, the outputs aren't reproducible (you can't tell if they came from `seqspec_v2` or `seqspec_v3`, which produce different trans-density counts). Ask whoever staged the outputs to re-mirror the full thing.

---

### I'm getting a "permission denied" error on Synapse

Three possible reasons:
1. **Not logged in.** Browser session expired; log back in at synapse.org.
2. **Token expired.** Personal access tokens have a lifetime; generate a new one in your Synapse settings.
3. **Pre-release file.** Some files are visible only to IGVF consortium members. Check with the dataset's lead.

---

### What's the right citation if I use a figure or number from these outputs in a talk/manuscript?

At minimum: cite the IGVF consortium, the CRISPR pipeline, perturbo (or sceptre, depending on the call), and the energy-distance method.

If the data are pre-publication, please loop in the relevant dataset lead (file a GitHub issue) before using them externally.

---

### Where do I report a bug in this doc?

File a GitHub issue at <https://github.com/adamklie/tf_perturb_seq/issues>. Including the exact sentence that confused you is more useful than "the doc is unclear."
