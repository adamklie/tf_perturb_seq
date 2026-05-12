# Outputs + step 3 (cutoffs)

Each run writes to a single `OUTPUT_DIR`. All filenames are fixed by `config1_2.json` `output_file_name_list`.

## After steps 0-2.1

```
<OUTPUT_DIR>/
├── inference_mudata.h5mu           # Downloaded MuData (kept for resubmits)
├── preprocessed.h5ad               # Step 0: filtered + PCA'd AnnData (X_pca in obsm)
├── annotation_table.csv            # Step 0: cell → intended_target_name lookup
├── gRNA_dict.pickle                # Step 0: target → list of cell barcodes
├── pca_dataframe.pickle            # Step 0: cached PCA (skipped on re-runs unless OVERWRITE_PCA_DICT=true)
├── config1_2.json                  # Steps 0/1/2 config (regenerated each invocation)
├── config3.json                    # Step 3 cutoffs (regenerated each invocation)
├── targeting_outlier_table.csv     # Step 1: which targeting gRNAs flagged as outliers
├── non_targeting_outlier_table.csv # Step 1: which NTC gRNAs flagged as outliers
├── discordance_gRNA.csv            # Step 1: DISCO test results per gRNA
├── pval_edist_full.csv             # *** STEP 2: HEADLINE OUTPUT ***
└── *.png                           # Step 2.1: diagnostic figures
```

## After step 3

```
<OUTPUT_DIR>/
├── target_by_target_matrix.csv     # Pairwise e-distance between significant targets
└── edist_embedding_info.csv        # 2D t-SNE coords + cluster labels
```

## `pval_edist_full.csv` — the headline

One row per target (gene or NTC). Columns (subset; check actual header):

| Column | Meaning |
|---|---|
| `intended_target_name` | Gene symbol or NTC label |
| `edist_mean` | Mean energy distance across `num_of_bg` background draws |
| `edist_std` | SD across background draws |
| `pvalue` | One-sided p-value from permutation test |
| `n_cells` | Cells assigned to this target after filtering |
| Per-bg breakdown | One column per background draw (rare; usually you want the summary) |

Typical screening cutoffs: `pvalue < 0.05` AND `edist_mean > 0.5`. Production-scale datasets typically yield 200–500 significant targets at these thresholds; benchmark datasets fewer.

```bash
# Quick hit count at default cutoffs:
awk -F',' 'NR>1 && $P<0.05 && $E>0.5' pval_edist_full.csv | wc -l   # P/E = pvalue/edist column indices
```

## Picking cutoffs for step 3

Step 3 builds a target × target distance matrix only for "significant" targets. Step 3 is expensive (O(N²) where N = significant targets), so cutoffs matter.

Rule of thumb:
- **Strict** (`pvalue < 0.01` AND `edist > 1.0`): ~50–100 targets, small clean matrix, fast.
- **Default** (`pvalue < 0.05` AND `edist > 0.5`): ~200–500 targets for production datasets.
- **Loose** (`pvalue < 0.1` AND `edist > 0.3`): ~500–1000 targets, matrix gets large.

Look at the step 2.1 plots before committing: a "knee" in the edist vs pvalue scatter often picks itself.

## Running step 3

Once cutoffs are picked, edit `<OUTPUT_DIR>/config3.json`:

```json
{
  "cutoff": {
    "pvalue_cutoff": 0.05,
    "edist_cutoff": 0.5
  }
}
```

Then re-invoke with `--run-step3`. Idempotent — only step 3 runs:

```bash
bash /cellar/users/aklie/projects/tf_perturb_seq/scripts/run_energy_distance_pipeline.sh \
  --mudata-path <OUTPUT_DIR>/inference_mudata.h5mu \
  --output-dir  <OUTPUT_DIR> \
  --run-step3
```

**Caveat:** the runner regenerates `config3.json` from its heredoc every invocation. If you've edited cutoffs there, your edits get clobbered. Two workarounds:

1. Edit the heredoc directly (copy `run_energy_distance_pipeline.sh` to `datasets/<ds>/bin/` and edit there per the frozen-scripts rule).
2. Bypass the runner for step 3 only:

   ```bash
   apptainer exec --nv \
     --bind /cellar:/cellar \
     /cellar/users/aklie/opt/containers/edist_pipeline.sif \
     bash -c "PYTHONPATH=/cellar/users/aklie/projects/tf_perturb_seq/external/energy_dist_pipeline/bin \
       python /cellar/users/aklie/projects/tf_perturb_seq/external/energy_dist_pipeline/bin/3_e_distance_among_regions.py \
       <OUTPUT_DIR>/config1_2.json <OUTPUT_DIR>/config3.json"
   ```

## `target_by_target_matrix.csv` and `edist_embedding_info.csv`

After step 3:

- `target_by_target_matrix.csv` — square distance matrix indexed by `intended_target_name`. Use for hierarchical clustering / heatmaps.
- `edist_embedding_info.csv` — per-target 2D t-SNE coordinates plus any cluster labels the step computed. Use for scatter plots colored by gene family / pathway.

These two files are the inputs for the cross-dataset edistance summary (`src/tf_perturb_seq/energy_dist/cross_dataset_edistance_summary.py`).

## Validation

There's a helper for sanity-checking step 2 outputs:

```bash
source /cellar/users/aklie/projects/tf_perturb_seq/.venv/bin/activate
python /cellar/users/aklie/projects/tf_perturb_seq/src/tf_perturb_seq/energy_dist/validate_edistance_outputs.py \
  --outdir <OUTPUT_DIR>
```

Checks file existence, column schema, row counts, and a few sanity invariants (NTC edist should cluster near 0; targeting edist should have a long right tail).

## What to commit

Per `docs/data/DATA.md` and `.gitignore`:

- **Tracked:** `<run>/energy_distance/scripts/`, `<run>/energy_distance/configs/`.
- **Gitignored:** `image/`, `logs/`, `*.csv`, `*.tsv`, `*.h5mu`, `*.pickle`.

All actual outputs (CSVs, pickles) are gitignored — they live on HPC and optionally mirror to Synapse for sharing.
