# Benchmark dir cleanup — 2026-04-28

Moves are reversible. All paths below are relative to `datasets/technology-benchmark_WTC11_TF-Perturb-seq/` unless noted.

Archive root: `scratch/2026_04_28_benchmark_cleanup/` (at repo root).

## What moved and why

### Tier 1 — clearly one-off

| From | To | Why |
|---|---|---|
| `bin/slurm_logs/` (entire dir, ~947K, 100+ files) | `scratch/2026_04_28_benchmark_cleanup/slurm_logs/from_benchmark_bin/` | Historical SLURM `.out`/`.err`/`.log` from past job runs |
| `_static/` (entire dir, 4 PNGs) | `scratch/2026_04_28_benchmark_cleanup/_static/from_benchmark_root/` | One-off PNG dump for an old presentation |
| `results/2026_03_25_benchmark_update (1).pdf` | `scratch/2026_04_28_benchmark_cleanup/meeting_artifacts/` | March 25 meeting deck |
| `results/speaker_notes.md` | `scratch/2026_04_28_benchmark_cleanup/meeting_artifacts/` | Notes for that meeting |
| `results/cross_tech_comparison/benchmark_images.zip` | `scratch/2026_04_28_benchmark_cleanup/meeting_artifacts/` | Bundle for sharing meeting figures |
| `results/cross_tech_comparison/benchmark_report.html` | `scratch/2026_04_28_benchmark_cleanup/meeting_artifacts/` | Self-contained meeting HTML |
| `results/cross_tech_comparison/benchmark_results.zip` | `scratch/2026_04_28_benchmark_cleanup/meeting_artifacts/` | Bundle for sharing meeting results |
| `results/cross_tech_comparison/benchmark_results_lite.zip` | `scratch/2026_04_28_benchmark_cleanup/meeting_artifacts/` | Lite version of above |

### Tier 2 — likely one-off but kept addressable

| From | To | Why |
|---|---|---|
| `bin/2_qc/scripts/extract_per_guide_capture.py` | `scratch/2026_04_28_benchmark_cleanup/benchmark_legacy_scripts/` | Marked "(legacy)" in CLAUDE.md |
| `bin/2_qc/scripts/run_extract_per_guide_capture.sh` | `scratch/2026_04_28_benchmark_cleanup/benchmark_legacy_scripts/` | Wrapper for the legacy script |
| `manifests/all_new_runs_inference_input.tsv` | `scratch/2026_04_28_benchmark_cleanup/benchmark_one_off_manifests/` | Batch inference input from a past iteration |
| `manifests/all_new_runs_pathway_input.tsv` | `scratch/2026_04_28_benchmark_cleanup/benchmark_one_off_manifests/` | Batch pathway input from a past iteration |
| `manifests/pending_calibration_input.tsv` | `scratch/2026_04_28_benchmark_cleanup/benchmark_one_off_manifests/` | Past calibration input |
| `manifests/rerun_qc_input.tsv` | `scratch/2026_04_28_benchmark_cleanup/benchmark_one_off_manifests/` | Past rerun input |
| `results/cross_tech_comparison/{filtering_flow*,filtering_funnel*,dataset_characteristics.png,dataset_summary_metrics.tsv,per_cell_metrics.png,per_lane_*,per_ms_read_counts.tsv,positive_control_*}` (17 files) | `scratch/2026_04_28_benchmark_cleanup/benchmark_top_level_aggregates/` | One-off cross-dataset aggregate figures from `generate_upstream_figures.py` for the March 25 update. Regenerable via that script |

### Removed-then-restored (the `_local.tsv` files are load-bearing)

I initially deleted three `_qc_paths_local.tsv` files thinking they were redundant copies of `_qc_paths.tsv` (the diff was empty). They were restored by copy because **the QC notebooks (`bin/2_qc/{1,2,3}_*.ipynb`) read `_qc_paths_local.tsv`** — the `_qc_paths.tsv` form is only what the sync script writes. For dashboard-available runs the two are byte-identical (just a copy); for fallback runs they differ.

Restored:
- `manifests/cleanser_500_mito_15pc_qc_paths_local.tsv` (copy of `_qc_paths.tsv`)
- `manifests/cleanser_800_mito_15pc_qc_paths_local.tsv` (copy of `_qc_paths.tsv`)
- `manifests/cleanser_knee2_mito_15pc_qc_paths_local.tsv` (copy of `_qc_paths.tsv`)

Untouched (already present, fallback form):
- `manifests/cleanser_unified_qc_paths_local.tsv`
- `manifests/sceptre_v11_qc_paths_local.tsv`

**Cleanup of the cleanup**: the sync script (`bin/1_get_data/scripts/sync_benchmark.sh`) writes only `_qc_paths.tsv`. After every sync, copy to `_qc_paths_local.tsv` until either (a) the notebooks are updated to fall back to `_qc_paths.tsv`, or (b) the sync script writes both forms. Worth fixing properly later.

## What stayed (the "core")

```
bin/
  0_cross_run_comparison.ipynb       # cross-run comparison (e.g. cleanser_unified vs sceptre_v11)
  6_cross_parameter_comparison.ipynb # cross-param comparison (the active sweep notebook)
  generate_report.py
  run_all_notebooks.sh
  1_get_data/                        # GCP sync orchestrator
  2_qc/                              # QC notebooks (1-3) + plotting scripts
  3_inference/                       # calibration + pathways + trans + TF enrichment
  4_integrate/                       # cross-dataset integration
  5_upload/                          # Synapse upload
manifests/                           # per-run input TSVs (kept)
results/
  cross_tech_comparison/
    cleanser_500_mito_15pc/          # per-run QC outputs
    cleanser_800_mito_15pc/
    cleanser_knee2_mito_15pc/
    cleanser_unified/
    sceptre_v11/
    cross_parameter/                 # cross-config comparison outputs
    cross_run/                       # cross-run comparison outputs
    generate_upstream_figures.py     # regenerator for the moved-out aggregates
    pipeline_stage_audit.md
    TAKEAWAYS.md
  integration/                       # 38GB of integrated h5ad outputs (left in place)
README.md, CLAUDE.md, SHARE_README.md
```

## Archive 2 (later same day, 2026-04-28) — superseded runs

Per Adam: cleanser_unified and sceptre_v11 are runs we've moved past. If we want a new comparison we'd run new ones. Plus the 38GB `results/integration/` directory contains integrated h5ad outputs derived from those runs, so it's archived alongside.

| From | To | Size |
|---|---|---|
| `datasets/<5 datasets>/runs/cleanser_unified/` | `scratch/2026_04_28_benchmark_cleanup/archived_runs/cleanser_unified/<dataset>/` | 6.8 GB |
| `datasets/<5 datasets>/runs/sceptre_v11/` | `scratch/2026_04_28_benchmark_cleanup/archived_runs/sceptre_v11/<dataset>/` | 5.5 GB |
| `results/cross_tech_comparison/cleanser_unified/` | `scratch/2026_04_28_benchmark_cleanup/archived_runs/cleanser_unified/_results/` | (in 6.8 GB above) |
| `results/cross_tech_comparison/sceptre_v11/` | `scratch/2026_04_28_benchmark_cleanup/archived_runs/sceptre_v11/_results/` | (in 5.5 GB above) |
| `manifests/cleanser_unified_qc_*.tsv` | `scratch/2026_04_28_benchmark_cleanup/archived_runs/cleanser_unified/` | small |
| `manifests/sceptre_v11_qc_*.tsv` | `scratch/2026_04_28_benchmark_cleanup/archived_runs/sceptre_v11/` | small |
| `results/integration/` (38 GB — h5ad outputs derived from the two archived runs) | `scratch/2026_04_28_benchmark_cleanup/archived_runs/integration_38GB/` | 38 GB |

Also updated:
- `datasets/technology-benchmark_WTC11_TF-Perturb-seq/bin/1_get_data/gcp_inference_mudata_manifest.csv` — removed 5 `Benchmark_cleanser_unified` rows
- `bin/6_cross_parameter_comparison.ipynb` — removed `cleanser_unified` from RUNS dict; markdown header now lists 5 active filter configs

The two scrublet runs were synced today and added to the manifest as new run labels:
- `scrublet_off_cleanser_800_mito_15pc` (← `Benchmark_scrublet_cleanser_800_mito_15pc`)
- `scrublet_on_sceptre_800_mito_15pc` (← `REAL_SCRUBLETS_RENAME_OTHER_Benchmark_scrublet_sceptre_800_mito_15pc`)

## Open structural questions (deferred)

- `bin/0_cross_run_comparison.ipynb` vs `bin/6_cross_parameter_comparison.ipynb` — both exist at top of `bin/`. Leaving both in place; if `0_*` becomes legacy, move it then
- `results/cross_tech_comparison/` is a confusing name now that we have multiple cross-axis comparisons (cross-tech, cross-param, cross-run, cross-version). Consider renaming to `results/comparisons/` with subdirs `per_run/`, `cross_param/`, `cross_run/`, etc. Defer until P1 lands so we don't break the active notebook
- `results/integration/` is 38GB — leave in place for now, decide separately whether to archive to a non-repo path

## Re-pulling moved items

```bash
mv /cellar/users/aklie/projects/tf_perturb_seq/scratch/2026_04_28_benchmark_cleanup/<archive_subdir>/<file> \
   /cellar/users/aklie/projects/tf_perturb_seq/datasets/technology-benchmark_WTC11_TF-Perturb-seq/<original_path>
```

The archive subdirs (`slurm_logs/`, `_static/`, `meeting_artifacts/`, `benchmark_legacy_scripts/`, `benchmark_one_off_manifests/`, `benchmark_top_level_aggregates/`) name the original location.
