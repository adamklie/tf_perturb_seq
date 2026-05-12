# Mirror recipes per bundle type

Each bundle type has its own mirror script. They all share a pattern: stage from canonical source → upload to Synapse → record `syn_id` in `synapse_paths.tsv`.

## CRISPR pipeline (`mirror_pipeline_outputs.py`)

**Source:** GCS — `gs://igvf-pertub-seq-pipeline-data/<DATASET>/<DATA_DATE>/outs/<RUN_LABEL>/` (`gcs_output_path` column in `reference/experimental_metadata.tsv`).

**Destination:** `syn64423137/2026_UTSW/datasets/<DATASET_ID>/crispr_pipeline/`

**Mirrored sub-bundles (3):**
- `pipeline_dashboard/` (~40 GB) — dashboard.html, inference_mudata.h5mu (dup), additional_qc/, evaluation_output/, figures/, guide_seqSpec_plots/, svg/
- `pipeline_info/` (~10 KB) — params_<timestamp>.json, nf_core_pipeline_software_versions.yml
- `pipeline_outputs/` (~20 GB) — inference_mudata.h5mu, perturbo_*_per_*_output.tsv.gz, sceptre/

**Recipe:**

```bash
# On HPC (recommended — bandwidth-bound; ~60 GB → ~60 min from GCS to Synapse)
zsh -ic '
python docs/jamborees/2026_UTSW/scripts/mirror_pipeline_outputs.py \
    --dataset Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq \
    --workdir /cellar/users/aklie/scratch/jamboree_pipeline_staging/Hon_CM \
    --dry-run
'
```

**Flags:**

| Flag | Default | Notes |
|---|---|---|
| `--dataset <id>` | required | Dataset ID (folder name under `datasets/`) |
| `--workdir <path>` | `/tmp/jamboree_pipeline_<id>` | **Use `/cellar/.../scratch/` on HPC** — `/tmp` is only 20 GB |
| `--dry-run` | off | List actions only |
| `--skip-download` | off | Reuse existing workdir (resume after partial upload) |

**Known quirks:**

- The GCS path is read from `reference/experimental_metadata.tsv` `gcs_output_path` — if a dataset's GCS run hasn't landed in that file, the mirror errors out clearly. Update the metadata first.
- `inference_mudata.h5mu` appears in both `pipeline_dashboard/` and `pipeline_outputs/` (intentional duplication per the schema). The mirror uploads both copies.
- Hon CM as of 2026-05-09: only `pipeline_dashboard/` + `pipeline_outputs/` are on Synapse (`syn74520421`); `pipeline_info/` is missing because Hon's Synapse upload didn't include it. Working with Weizhou to get the full bundle.
- Gersbach Hepatocyte (`syn70518849`): non-canonical layout (Sara's prior structure). Migration to canonical pending.

## cNMF (`mirror_cnmf_outputs.py`)

**Source:** HPC — `datasets/<DS>/<RUN>/cnmf/<cnmf_run_name>/Result/`

**Destination:** `syn64423137/2026_UTSW/datasets/<DATASET_ID>/cnmf/` (**FLAT** — no `<cnmf_run_name>/` nesting)

**Mirrored:**
- Per-K h5mu (default: selected-K only; ~5–6 GB at K=200)
- `Evaluation/<K>_2_0/` (all per-K evaluation TSVs)
- `Plot/k_selection/K-selection_panel_2.0.png` + `.svg`
- `Plot/Perturb_gene_<K>_2_0/` (per-TF PDFs)
- `Interpretation/Summary_table/<K>_2_0/cNMF_<K>_2_0.xlsx`
- `README.md` (per-run summary)

**Skipped by default:**
- Per-K h5mus other than selected K (saves ~25 GB per dataset; total ~30+ GB → ~7 GB)
- `Inference/adata/` raw outputs (large; selected K's h5mu has everything)
- `Data/` (input, reproducible)

**Recipe:**

```bash
# On HPC (cNMF outputs are HPC-resident)
zsh -ic '
python <REPO>/datasets/<DS>/<RUN>/cnmf/<cnmf_run_name>/Script/upload_to_synapse.py \
    --dataset DE \
    --upload-only-k200-h5mu \
    --dry-run
'
```

The cNMF mirror is wired as a **per-run** script under `<cnmf_run_name>/Script/upload_to_synapse.py` (not a centralized one in `docs/jamborees/.../scripts/`). This reflects the per-run nature of cNMF — each run has different selected K, different bug workarounds, etc.

**Flags:**

| Flag | Notes |
|---|---|
| `--dataset <DE\|ESC\|CM\|Hep>` | Selects the synapse target dataset code |
| `--upload-only-k<N>-h5mu` | Skip non-selected K h5mus (recommended; ~7–8 GB total upload) |
| `--include-all-k-h5mu` | Upload every K (~30+ GB; for archive runs only) |
| `--include-inference-adata` | Upload `Inference/adata/` (very large; rarely needed) |
| `--dry-run` | List actions only |

## Energy distance (`mirror_edistance_outputs.py`)

**Source:** HPC — `datasets/<DS>/<RUN>/energy_distance/` (canonical) or `datasets/<DS>/results/energy_distance/<RUN>/` (legacy).

**Destination:** `syn64423137/2026_UTSW/datasets/<DATASET_ID>/energy_distance/`

**Mirrored:**
- `pval_edist_full.csv` — headline output
- `target_by_target_matrix.csv` (if step 3 run)
- `edist_embedding_info.csv` (if step 3 run)
- `targeting_outlier_table.csv`, `non_targeting_outlier_table.csv`, `discordance_gRNA.csv`
- `config1_2.json`, `config3.json` (provenance)
- Step 2.1 PNGs

**Recipe:**

```bash
zsh -ic '
python docs/jamborees/2026_UTSW/scripts/mirror_edistance_outputs.py --dataset <DATASET_ID>
'
```

**Known caveats:**

- p-value calibration: Huangfu DE and Huangfu ESC are on Synapse but with the original p-values (not BH-adjusted) — see jamboree TODO.md for the recalibration task. The `pval_edist_full.csv` has a `pvalue` column; downstream WG1 e-distance summary applies BH on the fly.

## QC (`mirror_qc_outputs.py` — TBD)

Not yet a centralized script as of 2026-05-12. Per-dataset QC outputs are mirrored ad-hoc via `synapseclient` until a wrapper lands.

**Source:** HPC — `datasets/<DS>/<RUN>/qc/{mapping_gene,mapping_guide,intended_target}/`

**Destination:** `syn64423137/2026_UTSW/datasets/<DATASET_ID>/qc/`

**Mirrored:**
- All `*_metrics.tsv` (small; the cross-dataset summaries depend on these)
- All `*.png` (slideware)
- Plots and per-guide capture TSV

**Ad-hoc recipe (until a wrapper exists):**

```python
import os, synapseclient
syn = synapseclient.Synapse(silent=True)
syn.login(authToken=os.environ['SYNAPSE_AUTH_TOKEN'])

DATASET = 'Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq'
RUN = 'seqspec_v3'
LOCAL = f'/cellar/users/aklie/projects/tf_perturb_seq/datasets/{DATASET}/{RUN}/qc'
SYN_PARENT = '<syn_id_for_dataset_qc_folder>'  # create under .../datasets/<DATASET_ID>/qc/

# walk local and mirror, using `synapseutils.syncToSynapse` or manual File()s
```

The cross-dataset QC summary (`build_wg1_qc_summary.py`) reads from the **local repo** rather than Synapse — so the immediate jamboree need is to have these mirrored locally on HPC, not necessarily on Synapse.

## Common patterns

### Authentication failures

```
synapseclient.core.exceptions.SynapseAuthenticationError: 401
```

`SYNAPSE_AUTH_TOKEN` is unset or expired. Re-source `~/.zshrc` and run via `zsh -ic '...'`.

```
ServiceException: Anonymous users may not access this object
```

The PAT is for a different user / doesn't have access. Confirm the token belongs to a user with read+write on `syn64423137`.

### Recovering from partial uploads

- `--skip-download` reuses an existing `--workdir`, skipping the GCS pull. The upload step is idempotent on Synapse — re-uploading the same file overwrites.
- For partial uploads, delete the affected sub-folder on Synapse (via UI or `syn.delete(syn_id)`) and re-run. Synapse versions files on overwrite, but folder layout changes can leave orphans.

### Disk space

Bundles are large. On HPC with 200 GB usable in `/cellar/users/aklie/scratch`, you can stage at most ~3 datasets in parallel. Serialize mirrors to avoid filling disk.
