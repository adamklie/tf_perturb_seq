# Progress Tracking

Last updated: 2026-03-25

Target: **UTSW Jamboree, May 13-16, 2026**

## Benchmark Datasets (5 total)

All 5 share the WTC11 iPSC line and ~50-gene guide library.

| Dataset | CRISPR Pipeline | QC | Energy Distance | cNMF | Cross-tech Integration |
|---------|:-:|:-:|:-:|:-:|:-:|
| Hon benchmark | done | done | done (elated_almeida run) | done | done |
| Huangfu benchmark | done | done | **not started** | script exists, **not run** | done |
| Engreitz benchmark | done | done (batch ID issue) | **not started** | **no script** | done |
| Gersbach GEM-Xv3 | done | done | **not started** | script exists, **not run** | done |
| Gersbach HTv2 | done | done | **not started** | script exists, **not run** | done |

**Pipeline runs available:** All 5 have both `cleanser_unified` and `sceptre_v11` runs with `inference_mudata.h5mu`.

**Cross-tech benchmark** (`technology-benchmark_WTC11_TF-Perturb-seq/`):
- QC comparison report: done (`results/cross_tech_comparison/benchmark_report.html`)
- Integration (scVI): done for both cleanser and sceptre (`results/integration/`)
- Sceptre calibration: blocked (missing `log2_fc_std` outputs)

### Known issues
1. **Pipeline version mismatch**: CLEANSER and Sceptre used slightly different pipeline versions — will need eventual rerun for consistency
2. **Engreitz outlier**: Clear outlier in many QC metrics; believed to be a processing issue, not underlying data
3. **Trans calibration**: Calibrating trans results gives wildly different DEG calls for the same perturbation across datasets — root cause unknown, needs close investigation

### What's missing for benchmark "completion"
- **Energy distance**: Only Hon has this (old `elated_almeida` run, to be archived). 4 datasets need it. **Owner: Chikara** — has gotten it working, slow going.
- **cNMF**: Only Hon has results (old run, to be archived). Huangfu/Gersbach_GEM-Xv3/Gersbach_HTv2 have `6_run_cnmf.sh` but haven't been executed. Engreitz doesn't even have the script. **Owner: Alexandra Mo** — need a reproducible worked example that translates across MuDatas.

---

## Production Datasets (5 total, need >= 2 through full pipeline)

| Dataset | CRISPR Pipeline (standard) | H5MU | cNMF | Energy Distance | Notes |
|---------|:-:|:-:|:-:|:-:|-------|
| Hon cardiomyocyte | done (Hon lab internal) | no (h5ad only) | **not started** | **not started** | 137 GB processed data in HonLabInternal/; NOT run through standard CRISPR FG pipeline |
| Hon cardiomyocyte (pilot) | **not started** | no | no | no | Only has `1_generate_per_sample_metadata.sh` and a comparison notebook |
| Huangfu endoderm | done (external) | yes (3.8 GB) | **not started** | **not started** | Processed externally, not through standard pipeline scripts |
| Huangfu stemcell | done (external) | yes (2.4 GB) | **not started** | **not started** | Processed externally, has additional DE analysis |
| Gersbach hepatocyte | **config only** | no | no | no | Has nextflow.config + samplesheets from Oct 2025; no execution |

### Notes
- Hon cardiomyocyte and both Huangfu datasets have processed outputs from **lab-internal pipelines**, not the standardized CRISPR FG pipeline. The Huangfu datasets have MuData files ready for downstream analysis.
- Gersbach hepatocyte: **Sara Geraghty** is troubleshooting the CRISPR pipeline. Config exists from Oct 2025 but not yet executed. These are large datasets (~14k guides x ~1M cells) run on individual HPCs.
- Hon cardiomyocyte: **Weizhou** is troubleshooting the CRISPR pipeline. Also large-scale.
- **No production dataset has cNMF or energy distance results yet.**

### Strategy
Push the **two Huangfu datasets** (endoderm + stemcell) through the standardized CRISPR FG pipeline. The existing internal processing outputs are stale and non-uniform — they will be archived to `scratch/`. Pipeline scripts (steps 1-4) have been scaffolded; the blocking step is finding the IGVF portal accession IDs.

### Next steps
1. Find portal accession IDs for both Huangfu HUES8 datasets
2. Run `1_generate_per_sample_metadata.sh` to understand portal structure
3. Upload to GCP and configure pipeline
4. Run CRISPR pipeline (assess GCP cost vs. local SLURM)
5. Once pipeline outputs exist, run E-distance and cNMF

---

## Deadline Status

| Deadline | Target | Actual Status | Assessment |
|----------|--------|---------------|------------|
| **Mar 12** | Benchmark through CRISPR + cNMF + E-distance | CRISPR + QC done for all 5. cNMF + E-distance done for Hon only. | **MISSED** for 4/5 datasets on cNMF + E-distance |
| **Mar 26** | >= 2 production through CRISPR pipeline | 3 datasets have processed outputs (Hon cardio, Huangfu x2) but via internal pipelines, not standardized CRISPR FG pipeline | **AMBIGUOUS** — depends on whether internal processing counts |
| **Apr 9** | Iterate on additional production datasets | — | upcoming |
| **Apr 23** | >= 2 production through cNMF + E-distance | No production dataset has cNMF or E-distance | **not started** |
| **May 7** | Prep for jamboree | — | upcoming |
