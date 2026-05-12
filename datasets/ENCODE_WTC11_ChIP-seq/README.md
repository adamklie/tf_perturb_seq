# ENCODE WTC11 TF ChIP-seq

ENCODE TF ChIP-seq experiments in WTC11 (human iPSC) cells. Used to identify TF binding sites that overlap with perturbation targets in the TF Perturb-seq project.

## Files

| File | Description |
|------|-------------|
| `experiment_report.tsv` | ENCODE experiment metadata (90 TF ChIP-seq experiments). First line is the query URL, second line is the header. |
| `bed_manifest.tsv` | File manifest for bed narrowPeak files (IDR thresholded peaks) with download URLs. |
| `bigwig_manifest.tsv` | File manifest for bigWig signal files (fold change over control) with download URLs. |
| `bed_files.txt` | Download URLs for bed narrowPeak files. |
| `bigwig_files.txt` | Download URLs for bigWig files. First line is the query URL. |
| `benchmark_tfs.tsv` | Curated list of TFs for perturbation benchmarking with rationale and expected effect sizes. |
| `check_overlap.ipynb` | Checks overlap between benchmark TFs and ENCODE ChIP-seq targets. Found 5 overlapping TFs: AFF4, HMGA2, SMAD3, SMAD4, TCF7. |
| `downloads/` | Destination for downloaded files. |

## Downloading data

```bash
# Bed narrowPeak files
xargs -n 1 curl -O -J -L < bed_files.txt

# BigWig files (skip header line)
tail -n +2 bigwig_files.txt | xargs -n 1 curl -O -J -L
```

## Source

[ENCODE Portal](https://www.encodeproject.org/) -- TF ChIP-seq, WTC11, Homo sapiens, released.
