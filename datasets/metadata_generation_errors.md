# Metadata Generation Error Catalog

**Date:** 2026-01-20

**Script:** `scripts/generate_per_sample.py`

---

## Summary

All datasets failed due to **missing IGVF API credentials**.

| Dataset | Accession | Status | Error |
|---------|-----------|--------|-------|
| Huangfu_WTC11-benchmark_TF-Perturb-seq | IGVFDS4761PYUO | FAILED | Missing IGVF credentials |
| Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3 | IGVFDS6673ZFFG | FAILED | Missing IGVF credentials |
| Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2 | IGVFDS6237URFJ | FAILED | Missing IGVF credentials |
| Engreitz_WTC11-benchmark_TF-Perturb-seq | IGVFDS5057HJKP | FAILED | Missing IGVF credentials |
| Hon_WTC11-cardiomyocyte-differentiation_Pilot-TF-Perturb-seq | IGVFDS4389OUWU | FAILED | Missing IGVF credentials |

### Skipped (already working)
- Hon_WTC11-benchmark_TF-Perturb-seq
- Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq

---

## Root Cause

The script requires IGVF Portal API credentials to fetch analysis set metadata.

**Required:** Either:
1. Environment variables: `IGVF_API_KEY` and `IGVF_SECRET_KEY`
2. Or a keypair JSON file passed via `--keypair` flag

---

## Resolution

**Option 1:** Set environment variables before running:
```bash
export IGVF_API_KEY="your_key"
export IGVF_SECRET_KEY="your_secret"
```

**Option 2:** Provide keypair JSON file:
```bash
python3 scripts/generate_per_sample.py --keypair /path/to/keypair.json ...
```

**Note:** Must also activate the venv first:
```bash
source /Users/adamklie/Desktop/projects/tf_perturb_seq/.venv/bin/activate
```

---

## Full Error Output (Example: Huangfu)

```
==========================================
Generate Per-Sample Metadata
==========================================
Dataset:    Huangfu_WTC11-benchmark_TF-Perturb-seq
Accession:  IGVFDS4761PYUO
Output:     /Users/adamklie/Desktop/projects/tf_perturb_seq/datasets/Huangfu_WTC11-benchmark_TF-Perturb-seq/sample_metadata.csv

Generating per-sample CSV from analysis set 'IGVFDS4761PYUO'...
[ERROR] No credentials provided. Set IGVF_API_KEY and IGVF_SECRET_KEY or provide a keypair JSON.
Traceback (most recent call last):
  File "/Users/adamklie/Desktop/projects/tf_perturb_seq/scripts/generate_per_sample.py", line 290, in <module>
    main()
  File "/Users/adamklie/Desktop/projects/tf_perturb_seq/scripts/generate_per_sample.py", line 277, in main
    auth = get_auth(args.keypair)
  File "/Users/adamklie/Desktop/projects/tf_perturb_seq/scripts/generate_per_sample.py", line 74, in get_auth
    raise RuntimeError('No credentials provided. Set IGVF_API_KEY and IGVF_SECRET_KEY or provide a keypair JSON.')
RuntimeError: No credentials provided. Set IGVF_API_KEY and IGVF_SECRET_KEY or provide a keypair JSON.
```

---

## Previous Error (Resolved)

Initially all scripts failed with `ModuleNotFoundError: No module named 'requests'`. This was resolved by activating the venv at `.venv/`.
