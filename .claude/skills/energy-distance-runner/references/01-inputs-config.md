# Inputs + `config1_2.json` knobs

## MuData source (choose one)

| Flag | Source | Action |
|---|---|---|
| `--mudata-path <path>` | Local HPC file | Used as-is |
| `--gcs-mudata-path gs://...` | GCS | Downloaded to `<OUTPUT_DIR>/inference_mudata.h5mu` via `gcloud storage cp` |
| `--synapse-id synXXXXX` | Synapse | Downloaded via synapseclient; uses `SYNAPSE_AUTH_TOKEN` env var |

If `<OUTPUT_DIR>/inference_mudata.h5mu` already exists, the runner reuses it (idempotent). Delete it to force re-download.

**Synapse-only datasets (currently):** Hon CM uses `syn74522725` (2026_04_19_no_spacer pipeline_outputs). Once a full GCS bundle lands, switch to `--gcs-mudata-path`.

## Generated `config1_2.json` (steps 0, 1, 2, 2.1)

The runner writes this file in `<OUTPUT_DIR>/` each invocation. Default contents:

```json
{
  "output_file_name_list": { ... },
  "input_data": {
    "annotation_file": {
      "file_path": "<OUTPUT_DIR>/annotation_table.csv",
      "concatenate_key": "intended_target_name"
    },
    "h5ad_file": {
      "file_path": "<OUTPUT_DIR>/preprocessed.h5ad",
      "obsm_key": "X_pca"
    },
    "sgRNA_file": {
      "file_path": "<OUTPUT_DIR>/gRNA_dict.pickle"
    }
  },
  "gRNA_filtering": {
    "perform_targeting_filtering": true,
    "perform_nontargeting_filtering": true,
    "threshold_gRNA_num": 6,
    "combi_count": 4,
    "total_permute_disco": 1000,
    "combi_cell_num_max": 1000,
    "batch_num_basic": 120
  },
  "permutation_test": {
    "permute_per_bg": 1000,
    "num_of_bg": 20,
    "non_target_pick": 2000,
    "target_cell_num_max": 2000,
    "batch_num_basic": 200,
    "use_matched_bg": false
  },
  "aggregate": {
    "downsampling_maximum": 10000
  }
}
```

### Knobs that matter

| Block | Param | Default | When to change |
|---|---|---|---|
| `gRNA_filtering` | `threshold_gRNA_num` | 6 | Min cells per gRNA. Lower for sparse benchmark datasets (try 4 if many gRNAs are getting filtered out at 6) |
| `gRNA_filtering` | `combi_count` | 4 | Combinations of gRNAs to check; rarely changed |
| `gRNA_filtering` | `total_permute_disco` | 1000 | Permutations for DISCO outlier test. 100 for quick smoke-runs |
| `permutation_test` | `permute_per_bg` | 1000 | Permutations per background. Halve for speed, double for tighter p-values |
| `permutation_test` | `num_of_bg` | 20 | Number of NTC backgrounds sampled. More = tighter null estimate |
| `permutation_test` | `non_target_pick` | 2000 | Cells sampled per NTC background. Bound by available NTC cells |
| `permutation_test` | `target_cell_num_max` | 2000 | Cells per target. Cap above the largest target's cell count = no downsampling |
| `permutation_test` | `use_matched_bg` | false | Use batch-matched NTC backgrounds. True only if you suspect strong batch effects |
| `aggregate` | `downsampling_maximum` | 10000 | Cap for the aggregated PCA matrix |

### Editing the config

The runner regenerates `config1_2.json` from a hardcoded heredoc on every invocation. To use a custom config:

1. Run once to materialize the default, then edit it.
2. Add `OVERWRITE_PCA_DICT: false` is already in the default — set to `true` to force PCA recomputation.
3. Re-run with `--run-step3` — the heredoc still fires, so your edits get clobbered. **Workaround:** either patch the heredoc in `run_energy_distance_pipeline.sh` (copy to `datasets/<ds>/bin/` and edit per the frozen-scripts rule) or invoke the container directly bypassing the runner:

   ```bash
   apptainer exec --nv --bind /cellar:/cellar <CONTAINER_PATH> \
     bash -c "PYTHONPATH=<PIPELINE_BIN> python <PIPELINE_BIN>/2_e_distance_nontargeting.py <YOUR_EDITED_CONFIG>"
   ```

## Required Stage 2 columns in `inference_mudata.h5mu`

The preprocess script (`preprocess_mudata_local.py`) expects:

- `mdata.mod["gene"]` with `.obsm["X_pca"]` (computed by the CRISPR pipeline). If missing, preprocess will run PCA itself.
- `mdata.mod["guide"].var["intended_target_name"]` — used as the concatenation key.
- `mdata.mod["guide"].layers["guide_assignment"]` — to build `gRNA_dict.pickle`.

Missing any of these → preprocess fails with a clear `KeyError`. Re-run Stage 2 with the right config — don't patch around it.

## Cross-references

- Full upstream docs: https://github.com/Chikara-Takeuchi/energy_dist_pipeline
- Wrapper docs: https://github.com/Chikara-Takeuchi/energy_dist_TFperturb
- TFP3 wrapper notes: `docs/analysis/energy_dist/ENERGY_DISTANCE.md`
