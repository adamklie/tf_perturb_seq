# Scripts

Python scripts that produce or mirror the contents of `data/`. Most write to Synapse under `syn64423137/2026_UTSW/`; one writes a local portal snapshot.

> **Warning**: the mirror and upload scripts write to Synapse. Confirm the target before running. All Synapse-touching scripts require `SYNAPSE_AUTH_TOKEN` in the environment (set in `~/.zshrc` and re-source the shell before invoking).

Run all scripts from the jamboree root (`docs/jamborees/2026_UTSW/`).

## `query_igvf_portal.py`

- **Purpose**: snapshot IGVF portal state for the TF Perturb-seq Project (one JSON + TSV per record type, plus a manifest), and refresh the `latest` symlink.
- **Inputs**: IGVF portal REST API (`data.igvf.org/search/?type=<type>&collections=TF+Perturb-seq+Project&collections!=Benchmark`).
- **Output**: local snapshot under `docs/data/portal_snapshots/<utc-iso>/`. No Synapse writes.
- **Invocation**:
  ```
  uv run python data/scripts/query_igvf_portal.py
  uv run python data/scripts/query_igvf_portal.py --types MeasurementSet,AnalysisSet
  ```

## `mirror_pipeline_outputs.py`

- **Purpose**: mirror a dataset's CRISPR pipeline outputs (`pipeline_dashboard/`, `pipeline_info/`, `pipeline_outputs/`) from GCS to Synapse.
- **Inputs**: GCS (path read from `reference/experimental_metadata.tsv` -> `gcs_output_path`). Requires `gcloud auth` set to the IGVF service account.
- **Output**: Synapse folder `2026_UTSW/datasets/<dataset_id>/crispr_pipeline/`.
- **Invocation**:
  ```
  uv run python data/scripts/mirror_pipeline_outputs.py --dataset <dataset_id>
  uv run python data/scripts/mirror_pipeline_outputs.py --dataset <dataset_id> --dry-run
  ```

## `mirror_pipeline_outputs_hpc.py`

- **Purpose**: same as `mirror_pipeline_outputs.py`, but reads from a local HPC directory instead of GCS. Use when the run lives on the HPC rather than GCS.
- **Inputs**: HPC `--source-dir` containing `pipeline_dashboard/`, `pipeline_info/`, `pipeline_outputs/` (and optional `calibration/`, `tf/` with `--include-aux`).
- **Output**: Synapse folder `2026_UTSW/datasets/<dataset_id>/crispr_pipeline/`.
- **Invocation**:
  ```
  uv run python data/scripts/mirror_pipeline_outputs_hpc.py \
      --dataset <dataset_id> --source-dir <DIR> [--include-aux] [--readme-path <FILE>]
  ```

## `mirror_cnmf_outputs.py`

- **Purpose**: mirror a dataset's cNMF run outputs from the HPC to Synapse, following the inclusion rule in [`../schemas/cnmf.json`](../schemas/cnmf.json) (selected-k full data plus sweep-as-provenance).
- **Inputs**: HPC `--source-dir` (a PerturbNMF/cNMF result directory) and `--selected-k`.
- **Output**: Synapse folder `2026_UTSW/datasets/<dataset_id>/cnmf/`.
- **Invocation** (run on the HPC):
  ```
  uv run python data/scripts/mirror_cnmf_outputs.py \
      --dataset <dataset_id> --source-dir <DIR> --selected-k <K>
  ```

## `mirror_edistance_outputs.py`

- **Purpose**: mirror a dataset's energy-distance pipeline deliverables from the HPC to Synapse. Intermediates and the source MuData are skipped (see the script docstring for the whitelist).
- **Inputs**: HPC `--source-dir` (an energy-distance run output directory).
- **Output**: Synapse folder `2026_UTSW/datasets/<dataset_id>/energy_distance/`.
- **Invocation** (run on the HPC):
  ```
  uv run python data/scripts/mirror_edistance_outputs.py \
      --dataset <dataset_id> --source-dir <DIR>
  ```

## `upload_to_synapse.py`

- **Purpose**: idempotently upload a single local file to Synapse and return the resulting `syn` ID. Used by the mirror scripts and for ad-hoc uploads.
- **Inputs**: local `--src` file and a `--remote-path` under the Synapse parent (`syn64423137` by default).
- **Output**: Synapse file under the requested path. Existing folders are reused; existing files are versioned in place.
- **Invocation**:
  ```
  uv run python data/scripts/upload_to_synapse.py \
      --src <LOCAL_FILE> \
      --remote-path 2026_UTSW/<remote/path/to/file> \
      --description "<description>"
  ```
