# Nextflow `.config` spec

Per-run Nextflow config at `setup/configs/<DATASET>_<RUN_LABEL>.config`. Loaded via `nextflow run main.nf -c <config> ...`.

The config has four sections: `params {}`, `profiles {}`, `process {}` (resources per task), and singularity/tower/report blocks. Per-dataset edits are almost all in `params {}`.

## `params {}` — knobs you'll actually edit

### Chemistry / protocol flags

| Param | Type | Typical values | When to change |
|---|---|---|---|
| `ENABLE_DATA_HASHING` | bool | `true` (Hon HTO), `false` (others) | HTO/cell-hashing datasets only |
| `ENABLE_SCRUBLET` | bool | `true` / `false` | One axis of the benchmark sweep |
| `use_igvf_reference` | bool | `true` (default) | Set `false` only if pointing at a non-IGVF reference build |
| `is_10x3v3` | bool | `true` (3' v3) / `false` (5' HT v2, GEM-X, CC Perturb-seq) | Chemistry-specific |
| `reverse_complement_guides` | bool | `true` / `false` | Depends on guide construct orientation; ask the data producer |
| `spacer_tag` | string | `""` (most), `"TAGCTCTTAAAC"` (Hon cardio) | Set when guide construct has a spacer; empty otherwise |
| `DUAL_GUIDE` | bool | `false` | `true` for dual-guide libraries (rare in TFP3) |
| `REFERENCE_transcriptome` | string | `'human'` | Always human for TFP3 |
| `REFERENCE_gtf_download_path` | URL | GENCODE v46 default | Override only for special builds |
| `REFERENCE_gtf_local_path` | path | leave as `/path/to/...` placeholder | Used when you've pre-staged a GTF in `singularity-cache/` |

### QC

| Param | Type | Default | Notes |
|---|---|---|---|
| `QC_min_genes_per_cell` | int | `500` | Per-cell gene count floor |
| `QC_min_cells_per_gene` | float | `0.05` | Fraction (not count) of cells expressing a gene |
| `QC_pct_mito` | int | `15` (TFP3) / `20` (pipeline default) | Per-memory `[[project_qc_pct_mito]]` TFP3 uses 15% across the benchmark sweep |
| `QC_barcode_filter` | string | `'knee'`, `'cleanser'`, `'cleanser_500'`, `'cleanser_800'`, `'cleanser_extremes_200'`, `'cleanser_extremes_2000'`, `'cleanser_knee2'` | Main axis of the benchmark parameter sweep |

### Guide assignment

| Param | Type | Values | Notes |
|---|---|---|---|
| `Multiplicity_of_infection` | string | `'high'` / `'low'` | Affects guide-assignment defaults |
| `GUIDE_ASSIGNMENT_method` | string | `'sceptre'` / `'cleanser'` | Sweep axis |
| `GUIDE_ASSIGNMENT_capture_method` | string | `'CROP-seq'` | Almost always |
| `GUIDE_ASSIGNMENT_cleanser_probability_threshold` | float | `1` (no threshold) | Used when method=cleanser |
| `GUIDE_ASSIGNMENT_SCEPTRE_probability_threshold` | string/float | `'default'` or numeric | Used when method=sceptre |
| `GUIDE_ASSIGNMENT_SCEPTRE_n_em_rep` | string/int | `'default'` | Sceptre EM replicates |

### Inference

| Param | Type | Notes |
|---|---|---|
| `INFERENCE_method` | string | `'default'` (= sceptre + perturbo) |
| `INFERENCE_target_guide_pairing_strategy` | string | `'default'` |
| `INFERENCE_predefined_pairs_to_test` | path | Override for cis-only / custom pair lists |
| `INFERENCE_max_target_distance_bp` | int | `1000000` (1 Mb cis window) |
| `INFERENCE_SCEPTRE_side` | string | `'both'` / `'left'` / `'right'` |
| `INFERENCE_SCEPTRE_grna_integration_strategy` | string | `'union'` |
| `INFERENCE_SCEPTRE_resampling_approximation` | string | `'skew_normal'` |
| `INFERENCE_SCEPTRE_control_group` | string | `'default'` |
| `INFERENCE_SCEPTRE_resampling_mechanism` | string | `'default'` |
| `INFERENCE_SCEPTRE_formula_object` | string | `'default'` (or R formula string) |

### Network / dashboard

| Param | Notes |
|---|---|
| `NETWORK_custom_central_nodes` | TFs to highlight in dashboard network plots; `'undefined'` for default |
| `NETWORK_central_nodes_num` | Count of central nodes shown |

### Containers

```groovy
containers {
    base     = 'sjiang9/conda-docker:0.3'
    cleanser = 'ghcr.io/gersbachlab-bioinformatics/cleanser:1.2.1'
    sceptre  = 'sjiang9/sceptre-igvf:0.1'
    perturbo = 'ghcr.io/pinellolab/perturbo:sha-f3dc8ca'
    aria2    = 'biasofpriene/aria2c'
}
```

Pin these per run — bumping container tags between runs makes comparisons invalid.

### Google Cloud

```groovy
google_bucket  = 'gs://igvf-pertub-seq-pipeline-data'
google_project = 'igvf-pertub-seq-pipeline'
google_region  = 'us-central1'
```

### Resource caps (rarely change)

```groovy
max_cpus = 128
max_memory = 256.GB
```

## `profiles {}`

Three profiles: `local`, `slurm`, `google`. Always run with `-profile google` for TFP3. The `google` profile sets:

- `executor = 'google-batch'`
- `errorStrategy = retry on [137,143,50001,50002,50003,50006]; ignore on 0; terminate otherwise`
- `maxRetries = 3`
- `batch.spot = true` (cost; expect preemptions)
- `batch.maxSpotAttempts = 5`
- `batch.bootDiskSize = 100.GB`
- `workDir = gs://<bucket>/work`
- `singularity.cacheDir = gs://<bucket>/singularity-cache`

## `process { withName: ... }` blocks

These pin CPU/memory/machineType per task name. Defaults are usually fine. Notable:

- `mappingGuide|mappingHashing|mappingscRNA` — `n2-highmem-32`, 100GB/8 CPU base, scales with attempt.
- `guide_assignment_cleanser` — `n2-highmem-32`, 300GB, `batch.spot=false` (long-running, expensive to preempt), `batch.scratch=true`, `batch.disk=500GB`.
- `guide_assignment_sceptre|inference_sceptre` — `n2-highmem-64`, 200GB.
- `inference_perturbo|inference_perturbo_trans` — `a2-ultragpu-1g` (GPU), 40GB.
- `downloadReference` — `n1-highmem-32`, 100GB.

Edit these only if you're hitting OOM or want to switch GPU classes; check Tower for actual memory used first.

## Per-dataset cheat sheet

| Dataset | Chemistry | `ENABLE_DATA_HASHING` | `is_10x3v3` | `spacer_tag` |
|---|---|---|---|---|
| Hon_WTC11-benchmark | 10x 5' HT v2 + HTO | `true` | `false` | `""` |
| Hon_WTC11-cardio | 10x 5' HT v2 + HTO | `true` | `false` | `"TAGCTCTTAAAC"` |
| Huangfu_WTC11-benchmark | 10x 3' v3 | `false` | `true` | `""` |
| Huangfu_HUES8-DE/ESC | 10x 3' v3 | `false` | `true` | `""` |
| Engreitz_WTC11-benchmark | CC Perturb-seq | `false` | `false` | `""` |
| Gersbach_WTC11-benchmark_GEM-Xv3 | 10x GEM-X 3' | `false` | `true` | `""` |
| Gersbach_WTC11-benchmark_HTv2 | 10x HT v2 | `false` | `false` | `""` |
| Gersbach_WTC11-hepatocyte | 10x HT v2 | `false` | `false` | `""` |

When in doubt, copy the latest config from a sibling dataset with the same lab/chemistry and adjust.

## Cross-references

- IGVF CRISPR_Pipeline source of truth for params: https://github.com/IGVF/CRISPR_Pipeline
- Existing per-dataset configs: `datasets/*/setup/configs/*.config`
