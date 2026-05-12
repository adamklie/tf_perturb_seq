# Module internals + known bandaids

What each of the three Python modules does and the non-obvious bits that bite when something fails.

## `mapping_gene.py`

**Reads:** `inference_mudata.h5mu` → `mdata.mod["gene"]`. Needs `.obs["batch"]`, `.obs["total_gene_umis"]`, `.obs["percent_mito"]`, and `.obs["num_expressed_genes"]`.

**Fallback:** if `num_expressed_genes` is missing, derives it from `np.expm1(obs["log1p_n_genes_by_counts"])`. This was added because some GCP pipeline runs don't emit `num_expressed_genes`.

**Per-batch breakdown:** groups by `obs["batch"]` and emits one row per batch plus an `overall` row.

**Plots:**

- Knee plot: ranked `total_gene_umis` on log-log.
- Histograms: 3-panel (UMI, genes, mito) with vertical median lines.
- Histograms by batch: same 3 panels, seaborn `histplot(element="step", common_norm=False)` colored by batch.
- Cells per batch: bar chart.

**Common edge case:** `obs["batch"]` is a string column. If a dataset has integer batch IDs, they're coerced to string for the groupby. If `obs["batch"]` is missing entirely, the module errors out with a `KeyError`.

## `mapping_guide.py`

**Reads:** `mdata.mod["guide"]`. Needs `.obs["batch"]`, `.obs["total_guide_umis"]`, `.var["label"]`, `.var["gene_name"]`, `.layers["guide_assignment"]` (binary or sparse).

**Pre-computation in `compute_guide_assignment_counts()`:**
- Binarizes `layers["guide_assignment"]` (sparse → dense bool).
- `n_guides_per_cell = (A > 0).sum(axis=1)` → stored on `obs`.
- `n_cells_per_guide = (A > 0).sum(axis=0)` → stored on `var`.

**Per-batch breakdown:** groups by `guide.obs["batch"]`.

### Engreitz bandaid (important)

If `mdata.mod["gene"].obs["batch"]` and `mdata.mod["guide"].obs["batch"]` have **zero overlap**, the module copies gene batch labels onto guide obs (cell-barcode-matched). This handles the Engreitz dataset where RNA and guide libraries are registered as different IGVF analysis set accessions, so each modality has its own batch IDs.

**Symptom that the bandaid fired:** in `mapping_guide_metrics.tsv`, the overall row has `n_batches=1` even though `mapping_gene_metrics.tsv` shows multiple batches. The bandaid logs to stderr — check the SLURM `.err` file for a line like `WARN: zero batch overlap between gene and guide; copying gene batches onto guide obs`.

**When to override:** never automatically. If the bandaid is wrong for a future dataset (e.g., the user intentionally wants per-modality batches), comment it out manually in `mapping_guide.py` for that one run. Don't generalize.

### `per_guide_capture` table

For each guide in `var`:
- `n_cells_detected = (A > 0)[:, g].sum()`
- `frac_cells_detected = n_cells_detected / n_cells`
- `total_umi`, `mean_umi`, `median_umi`, `std_umi`, `max_umi` computed from the **raw** (non-binarized) assignment matrix.

Sorted by `n_cells_detected` desc.

## `intended_target.py`

**Reads:**
- `mdata.mod["guide"].var` — needs `intended_target_name`, `gene_name`, `label`.
- `mdata.uns["trans_per_guide_results"]` — per-guide trans-test DataFrame with `guide_id`, `gene_id` (or `gene_name`), `log2_fc`, `p_value`.

**Workflow:**

1. **Build intended target map:** `guide.var["intended_target_name"]` (the TF the guide is supposed to perturb) → `guide.var["gene_name"]` (the gene tested in the trans table). For targeting guides these match; for NTCs `intended_target_name` is empty/NaN.

2. **Load trans table** from `uns`.

3. **Filter to intended-target rows:** rows where the guide's intended target gene matches the tested gene. NTC guides are excluded here (no intended target).

4. **Knockdown metrics:**
   - `log2_fc <= log2(0.4)` (≥60% knockdown).
   - `p_value < 0.05`.
   - Reported as counts + fractions in `_metrics.tsv`.

5. **Balanced evaluation table** (for AUROC/AUPRC):
   - **Positives:** targeting guides paired with their intended target.
   - **Negatives:** non-targeting guides paired with the **same set** of genes (downsampled to match positives count per gene).
   - **Score:** `1 - p_value` (higher = more knockdown evidence).

6. **AUROC / AUPRC** computed from scikit-learn on the balanced table.

**Plots:**

- Volcano: `log2_fc` vs `-log10(p_value)`, colored by `label`, top N most-significant points text-labeled.
- log2FC distribution: histogram colored by label, vertical lines at 0 and `log2(0.4)`.
- ROC/PR curves: 2-panel, AUROC/AUPRC annotated.

### Edge case: missing `intended_target_name`

If `guide.var["intended_target_name"]` is missing or all empty, the module raises a clear error. Check the IGVF pipeline output schema — some older runs use `target_name` instead. The schema fix is upstream (re-run the CRISPR pipeline with a fresh build) rather than patching the QC module.

### Edge case: zero NTC guides

If there are no `non_targeting` rows in `guide.var["label"]`, AUROC/AUPRC can't be computed (no negatives) and the module emits NaN in those columns. The knockdown counts still work.

## Shared assumptions

All three modules assume:

- MuData has both `gene` and `guide` modalities under those exact keys.
- `.obs["batch"]` exists on both modalities (the bandaid handles guide's case).
- The IGVF pipeline emitted `mdata.uns["trans_per_guide_results"]` (needed only by `intended_target`).

If any of these fails, the module errors out with a `KeyError` — not silently producing garbage. That's intentional.

## Where this came from

For the architectural rationale and history, see `src/tf_perturb_seq/qc/plan.md`. The modules were originally three separate scripts; they were standardized to a common CLI surface (`--input`, `--outdir`, `--prefix`) so the array driver could call them uniformly.
