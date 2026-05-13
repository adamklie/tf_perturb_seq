# Energy-distance p-value calibration concern

**Status**: ⚠ Huangfu DE + ESC runs complete on Synapse, p-values are anti-conservative (every negative-control target gets `pval_mean = 0`). Hon CM run also complete and shows healthy calibration. HTv2 benchmark also healthy. Need a fix recommendation before re-running Huangfu.

**Owner of fix**: us (apply whatever preprocess change Chikara recommends, re-run on HPC) ± Chikara Takeuchi (upstream guidance on what produces calibrated p-values for the affected datasets).

## TL;DR

We ran the energy-distance pipeline on five Gersbach/Hon/Huangfu runs. Three are calibrated, two are not:

- **Calibrated**: Gersbach HTv2 benchmark (both Chikara's reference and our re-run), Hon WTC11 cardiomyocyte production
- **Anti-conservative**: Huangfu HUES8 definitive-endoderm production, Huangfu HUES8 embryonic-stem-cell production

Same wrapper, same container, same parameters, same library, same scale of cell counts and target counts across the production runs — but the two Huangfu runs come out broken and the Hon CM run does not.

## Observations

### Distribution numbers per run

Numbers from each run's `pval_edist_full.csv` (schema in [`schemas/energy_distance.json`](../../schemas/energy_distance.json)):

| Run | Targets | Cells/target (median) | Targeting `distance_mean` median | NC `distance_mean` median | NC vs targeting ratio | `pval_mean` median | NCs with `pval_mean=0` |
|---|---:|---:|---:|---:|---:|---:|---:|
| HTv2 verified reference (Chikara, [`syn74381183`](https://www.synapse.org/Synapse:syn74381183)) | 64 | 2000 | 1.78 | _no NC class in run_ | — | (healthy) | _n/a_ |
| HTv2 our rerun (`cleanser_800_mito_15pc`, 2026-05-10) | 65 | 2000 | 19.46 | _no NC class in run_ | — | 0.0316 | _n/a_ |
| **Hon WTC11 cardiomyocyte** ([`syn74897350`](https://www.synapse.org/Synapse:syn74897350)) | **2030** | up to 2000 (range 5-2000) | **1.30** | (calibrated) | (calibrated) | **0.339** ✅ | (calibrated) |
| Huangfu HUES8 DE ([`syn74883327`](https://www.synapse.org/Synapse:syn74883327)) | 2267 | 2000 | 93.21 | 91.70 | 1.02× | **0** | **100 / 100** |
| Huangfu HUES8 ESC ([`syn74883475`](https://www.synapse.org/Synapse:syn74883475)) | 2267 | 2000 | 527.33 | 528.25 | 1.00× | **0** | **100 / 100** |

Side observation worth flagging: in both Huangfu runs, positive controls have *lower* median distance than NCs (DE: pos-ctrl 77 vs NC 92; ESC: pos-ctrl 500 vs NC 528).

### Plots (Huangfu DE + ESC vs HTv2 reference)

#### `pval_mean` distribution by target type
![pval_mean by type](01_pval_mean_by_type.png)

NC bars pile at 0 in both Huangfu runs.

#### `distance_mean` distribution by target type
![distance_mean by type](02_distance_mean_by_type.png)

NC and targeting `distance_mean` distributions overlap near-perfectly in the Huangfu runs.

#### Volcano: `distance_mean` vs `−log10(pval_mean)`
![volcano by type](03_volcano_by_type.png)

NCs (blue) sit at the top of the volcano next to targeting (red) in the Huangfu runs.

#### Cross-run distance scale comparison
![distance scale comparison](04_distance_scale_comparison.png)

Huangfu distances are 1-2 orders of magnitude larger than HTv2 reference. Hon CM distances are in the same regime as HTv2.

## What we have ruled out

### Pipeline / wrapper / preprocess deviation
We re-ran the e-distance pipeline on Gersbach HTv2 benchmark with our wrapper (`cleanser_800_mito_15pc`, 65 targets) and compared bit-for-bit against Chikara's HTv2 run on Synapse [`syn74895081`](https://www.synapse.org/Synapse:syn74895081):

```
65/65 targets overlap, 0 unique to either side
Pearson distance correlation: 1.0000
Pearson pval correlation:     1.0000
Per-target abs diff:          median 0, max 0
```

### Wrong background source
Step 2 (`2_e_distance_nontargeting.py`, line 55) reads non-targeting outliers explicitly: `clear_nt_sgRNA_list = nontargeting_outlier_df[nontargeting_outlier_df["pval_outlier"]>0.05].index.tolist()`. The permutation null is built from `non-targeting`-typed gRNAs only. Negative controls are tested as targets, not used as background.

### Guide-metadata labeling
Both libraries use OR (olfactory receptor) gene-targeting gRNAs as the de-facto negative-control class — just labeled differently:

| Library | NC labeling | OR-targeting gRNAs | OR pval_mean median | OR frac p<0.05 |
|---|---|---|---|---|
| HTv2 (416 gRNAs) | Class doesn't exist | 54, labeled `type=targeting` | 0.0654 | 44% |
| Huangfu DE (14k gRNAs) | Explicit `type=negative control` | 592 of 598 NCs match `^OR[digit]` | 0 | 100% |

The labeling difference doesn't affect the test mechanically (cells go through the same comparison either way). HTv2's OR-targeting subset (44% sig) sits below HTv2 real targeting (63% sig), showing calibration works there even without an explicit NC class.

### Scale alone
Earlier we suspected the breakdown was driven by the larger non-targeting pool in Huangfu (600 NT gRNAs vs HTv2's 30) creating a tighter permutation null. Hon CM has the same 600 NT gRNAs, the same 14k gRNA library, the same ~2k unique target groups, ~270k cells — i.e. matched scale to Huangfu — and is calibrated. So scale alone is not the cause.

## Open question

Three calibrated runs (HTv2 Chikara, HTv2 our rerun, Hon CM) and two broken runs (Huangfu DE, Huangfu ESC) using the same pipeline, same wrapper, same params, and (between Hon CM and Huangfu) the same library. **What's different about the Huangfu data that breaks the test?**

Possibilities to consider — without picking a working hypothesis:
- Cell-state heterogeneity across batches / differentiation timepoints (Huangfu is HUES8-derived definitive endoderm + embryonic stem cells; Hon is WTC11 cardiomyocyte, terminally differentiated)
- Cell-cycle / lineage gradient variance dominating the PCA at this scale
- A QC / filtering threshold mismatch upstream of the e-distance preprocess
- Something specific to how Huangfu's `inference_mudata.h5mu` was constructed (kallisto / sceptre version, MOI, batch effects)
- Something else Chikara has seen at production scale

## Setup details

**Container**: `docker.io/takechikara/energy_distance_env:latest` (apptainer .sif at `/cellar/users/aklie/opt/containers/edist_pipeline.sif`).

**Pipeline**: [`external/energy_dist_pipeline`](../../../../../external/energy_dist_pipeline/) submodule pinned at `5821450a1eacfddb3c83be8c69fd593eaf76a61c` (post-downsampling-fix main).

**Wrapper**: [`scripts/run_energy_distance_pipeline.sh`](../../../../../scripts/run_energy_distance_pipeline.sh) (based on Chikara's [`external/energy_dist_TFperturb`](../../../../../external/energy_dist_TFperturb/) template).

**Preprocess**: matches upstream `preprocess_mudata.py` step-for-step, in [`src/tf_perturb_seq/crispr_pipeline/preprocess_mudata_local.py`](../../../../../src/tf_perturb_seq/crispr_pipeline/preprocess_mudata_local.py):
1. Take `gene` modality from MuData
2. `sc.pp.filter_genes(min_counts=1)`
3. `sc.pp.normalize_total`
4. `sc.pp.log1p`
5. `sc.pp.scale`
6. `sc.tl.pca(n_comps=50)`

The only structural deviations from upstream:
- We read MuData from local disk (Huangfu MuData isn't on Synapse) rather than via `--synapse-id`
- The annotation table writes the promoter-format string `<ENSG>|<chr>:<start>-<end>` into `intended_target_name` (rather than into a separate `intended_target_promoter` column), with config `concatenate_key="intended_target_name"`. Same effective grouping as upstream.

**Step-1/2/2.1 config** — identical to upstream `external/energy_dist_pipeline/config.json` defaults:
- gRNA filtering: `threshold_gRNA_num=6`, `combi_count=4`, `total_permute_disco=1000`, `combi_cell_num_max=1000`, `batch_num_basic=120`
- Permutation: `permute_per_bg=1000`, `num_of_bg=20`, `non_target_pick=2000`, `target_cell_num_max=2000`, `batch_num_basic=200`, `use_matched_bg=false`

**Inputs (Huangfu DE/ESC, Hon CM)**: 12,934 targeting + 600 non-targeting + 19 positive controls + 598 negative controls (all from IGVF library `IGVFFI8270UPKB`, pools A-D). After grouping by `intended_target_promoter`, ~2k unique target regions per dataset.

## Acceptance criteria for the fix

After applying whatever fix Chikara recommends and re-running:

- [ ] NC `pval_mean` distribution roughly uniform on [0, 1] (not piled at 0) for Huangfu DE and ESC
- [ ] NC `distance_mean` median **less than** targeting `distance_mean` median by a meaningful margin
- [ ] Volcano plot shows NCs separating from targeting (NCs near origin, targeting spread out)
- [ ] Hon CM and HTv2 results unchanged (or, if the fix changes them, still calibrated)
- [ ] All 4 layers of [`src/tf_perturb_seq/edistance/validate_edistance_outputs.py`](../../../../../src/tf_perturb_seq/edistance/validate_edistance_outputs.py) still PASS

## Pointers

| Object | Path |
|---|---|
| Validator (4 layers: file presence, CSV schema, value-range sanity, schema-identity vs HTv2 reference) | [`src/tf_perturb_seq/edistance/validate_edistance_outputs.py`](../../../../../src/tf_perturb_seq/edistance/validate_edistance_outputs.py) |
| Pipeline runner | [`scripts/run_energy_distance_pipeline.sh`](../../../../../scripts/run_energy_distance_pipeline.sh) |
| Preprocess (the file we'd edit for any preprocess-side fix) | [`src/tf_perturb_seq/crispr_pipeline/preprocess_mudata_local.py`](../../../../../src/tf_perturb_seq/crispr_pipeline/preprocess_mudata_local.py) |
| Per-dataset entrypoints | `datasets/<dataset>/5_run_energy_distance.sh` |
| Mirror script | [`scripts/mirror_edistance_outputs.py`](../../scripts/mirror_edistance_outputs.py) |
| Schema | [`schemas/energy_distance.json`](../../schemas/energy_distance.json) |
| Analysis-level walkthrough | [`docs/analysis/ENERGY_DISTANCE.md`](../../../../analysis/ENERGY_DISTANCE.md), [`docs/analysis/ENERGY_DISTANCE_OUTPUTS.md`](../../../../analysis/ENERGY_DISTANCE_OUTPUTS.md) |
| Upstream pipeline | [`external/energy_dist_pipeline/`](../../../../../external/energy_dist_pipeline/) (pinned to `5821450`) |
| Upstream wrapper | [`external/energy_dist_TFperturb/`](../../../../../external/energy_dist_TFperturb/) (Chikara's template) |

## Synapse links

- Huangfu DE bundle: [`syn74883327`](https://www.synapse.org/Synapse:syn74883327)
- Huangfu ESC bundle: [`syn74883475`](https://www.synapse.org/Synapse:syn74883475)
- Hon CM bundle (calibrated reference): [`syn74897350`](https://www.synapse.org/Synapse:syn74897350)
- HTv2 reference (Chikara, newer, used for bit-perfect reproducibility check): [`syn74895081`](https://www.synapse.org/Synapse:syn74895081)
- HTv2 reference (Chikara, older, used for schema cross-check): [`syn74381167`](https://www.synapse.org/Synapse:syn74381167)
