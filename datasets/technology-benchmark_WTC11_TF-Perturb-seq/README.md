# WTC11 TF-Perturb-seq Technology Benchmark

Master tracking document for running the CRISPR Pipeline on WTC11 benchmark datasets.

## Dataset Overview

| Dataset | Lab | Technology | Hashing | MOI | Pipeline Status |
|---------|-----|------------|---------|-----|-----------------|
| [Hon](../Hon_WTC11-benchmark_TF-Perturb-seq/) | Hon | 10x 5' | ✅ Yes | High | ✅ **Complete** |
| [Huangfu](../Huangfu_WTC11-benchmark_TF-Perturb-seq/) | Huangfu | 10x 3' v3 | ❌ No | High | ⚠️ Partial |
| [Gersbach GEM-Xv3](../Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/) | Gersbach | 10x GEM-X 3' | ❌ No | Low | ❌ Not started |
| [Gersbach HTv2](../Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/) | Gersbach | 10x HT v2 | ❌ No | Low | ❌ Not started |

---

## Critical Path to Running Pipeline

### 1. Gersbach Datasets (GEM-Xv3 & HTv2) - Blocked

**Blocker**: Seqspecs and guide metadata not available locally

| Task | Priority | Status |
|------|----------|--------|
| Download seqspecs from IGVF portal | P0 | ❌ |
| Confirm seqspecs are pipeline-compatible | P0 | ❌ |
| Download guide metadata (`IGVFFI5765HMZH`) | P0 | ❌ |
| Unzip barcode onlists | P1 | ❌ |
| Upload files to GCS | P1 | ❌ |
| Update CSV paths to GCS | P1 | ❌ |
| Run pipeline | P2 | ❌ |

**Key Question**: Are portal seqspecs usable as-is, or do they need `.formatted.fixed` / `.simple_formatted` processing?

### 2. Huangfu Dataset - Partially Ready

**Blocker**: Files in wrong GCS location (`beer_test/`)

| Task | Priority | Status |
|------|----------|--------|
| Seqspecs uploaded to new location | - | ✅ Done |
| Move barcode onlist to production location | P0 | ❌ |
| Move guide metadata to production location | P0 | ❌ |
| Update CSV paths | P1 | ❌ |
| Clarify guide metadata version (v1 vs v3_rc) | P1 | ❌ |
| Run pipeline | P2 | ❌ |

**Key Question**: Why does Huangfu use different guide metadata (`v1.tsv`) and barcode onlist (`IGVFFI4695IKAL`)?

### 3. Hon Dataset - Complete (Maintenance)

| Task | Priority | Status |
|------|----------|--------|
| Pipeline run | - | ✅ Done |
| Update guide metadata to latest spec | P2 | ❌ |
| Rerun with updated metadata | P2 | ❌ |
| Get CellRanger/PySpade outputs from Mpathi | P3 | ❌ |

---

## Shared Resources

### Guide Metadata
| Accession | Version | Used By |
|-----------|---------|---------|
| `IGVFFI5765HMZH` | v3_rc | Hon, Gersbach GEM-Xv3, Gersbach HTv2 |
| `benchmark_guide_metadata_v1.tsv` | v1 | Huangfu |

### Barcode Onlists
| Accession | Used By |
|-----------|---------|
| `IGVFFI9487JPEN` | Hon, Gersbach HTv2 |
| `IGVFFI1697XAXI` | Gersbach GEM-Xv3 |
| `IGVFFI4695IKAL` | Huangfu |

---

## GCS Directory Structure (Target)

```
gs://igvf-pertub-seq-pipeline-data/
├── Hon_WTC11-benchmark_TF-Perturb-seq/           ✅ Done
│   └── 2026_01_18/
│       ├── rna_seqspec.yml
│       ├── guide_seqspec.yml
│       ├── hash_seqspec.yml
│       └── [data files]
├── Huangfu_WTC11-benchmark_TF-Perturb-seq/       ⚠️ Partial
│   ├── seqspecs/                                 ✅ Done
│   ├── onlists/                                  ❌ TODO
│   └── guide_metadata/                           ❌ TODO
├── Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/  ❌ TODO
│   ├── seqspecs/
│   ├── onlists/
│   └── guide_metadata/
└── Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/    ❌ TODO
    ├── seqspecs/
    ├── onlists/
    └── guide_metadata/
```

---

## Configuration Differences

| Parameter | Hon | Huangfu | Gersbach |
|-----------|-----|---------|----------|
| `GUIDE_ASSIGNMENT_method` | sceptre | sceptre | cleanser |
| `Multiplicity_of_infection` | high | high | low |
| `ENABLE_DATA_HASHING` | true | false | false |
| `is_10x3v3` | false | **true** | false |
| `reverse_complement_guides` | true | - | true (in params.config) |

---

## Immediate Next Steps

1. **Download and test Gersbach seqspecs** from IGVF portal
   - Determine if `.formatted.fixed` processing is needed
   
2. **Relocate Huangfu files** from `beer_test/` to production location
   - Clarify guide metadata version discrepancy

3. **Standardize configurations** across datasets
   - Decide on SCEPTRE vs CLEANSER for guide assignment
   - Reconcile MOI settings

---

## Links

- **IGVF Data Portal**: https://data.igvf.org/
- **GCS Bucket**: `gs://igvf-pertub-seq-pipeline-data`
- **Pipeline Repo**: (add link)

---

*Last updated: 2026-01-23*
