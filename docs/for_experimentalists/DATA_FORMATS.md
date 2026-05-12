# Data formats — what is this file?

Quick reference for the file types you'll encounter. None of this is exotic; once you know the format, you can open it in Python, R, or sometimes Excel.

---

## TSV / CSV / `.tsv.gz` / `.csv.gz`

Plain-text tables, one row per record, columns separated by tabs (TSV) or commas (CSV). The `.gz` versions are gzip-compressed; any modern tool handles them transparently.

**Open without code:** Drag-drop into Excel/Numbers (it'll prompt for delimiter). Or in a terminal: `zcat file.tsv.gz | head` to peek at the first 10 lines.

**Open in Python:** `pandas.read_csv("file.tsv.gz", sep="\t")` (pandas detects the gzip automatically).

**Open in R:** `read.delim("file.tsv.gz")` (or `readr::read_tsv` if you prefer the tidyverse).

We use TSV (tabs) for everything in this project because gene symbols and program names can contain commas.

---

## `.h5ad` — AnnData

A single-cell expression matrix, stored in HDF5. One `.h5ad` file contains:

- The cell × gene count matrix (`.X`).
- Per-cell metadata in `.obs` (cell type, guide identity, QC flags, …).
- Per-gene metadata in `.var` (gene symbol, ENSG, biotype, …).
- Embeddings / dim-reductions in `.obsm` (UMAP, PCA, …).
- Misc unstructured metadata in `.uns`.

**Open without code:** Not really practical. The h5ad is binary; you need an AnnData-aware tool.

**Open in Python:**
```python
import anndata
adata = anndata.read_h5ad("file.h5ad")
print(adata)              # summary
print(adata.obs.head())   # per-cell metadata
print(adata.var.head())   # per-gene metadata
```

**Open in R:** Use the `anndataR` package or convert to Seurat via `zellkonverter`/`SeuratDisk`.

---

## `.h5mu` — MuData

A multi-modal extension of AnnData. One `.h5mu` file holds **multiple AnnData objects together** under one key per modality. For our pipelines:

- `mdata.mod["gene"]` — the RNA AnnData (cells × genes).
- `mdata.mod["guide"]` — the guide AnnData (cells × guides; cell-by-guide assignments).
- `mdata.mod["hashing"]` — for HTO-multiplexed experiments only (cell × HTO barcodes).

The cells line up across modalities (same cell barcodes, same order).

**Open in Python:**
```python
import mudata
mdata = mudata.read_h5mu("inference_mudata.h5mu")
print(mdata)                          # summary of all modalities
print(mdata.mod["gene"].shape)        # (n_cells, n_genes)
print(mdata.mod["guide"].shape)       # (n_cells, n_guides)
print(mdata.obs.head())               # shared cell metadata
```

---

## `.pkl` — Python pickle

A serialized Python object. The IGVF portal occasionally hosts raw feature-barcode matrices as pickles. You'll only encounter these as upstream inputs; we don't ship them as deliverables.

**Open in Python:**
```python
import pickle
with open("file.pkl", "rb") as f:
    obj = pickle.load(f)
```

⚠ **Security note**: never `pickle.load` a file you don't trust — pickles can execute arbitrary code.

---

## `.hdf5` / `.h5` — Generic HDF5

A hierarchical binary container. The IGVF portal stores feature-barcode matrices (from `cellranger`) as `.hdf5`. Inside, there's a `/matrix` group with sparse data.

**Open in Python:** Use `scanpy.read_10x_h5("file.hdf5")` — returns an AnnData.

---

## BED files

A 3+ column tab-delimited format for genomic coordinates: `chrom`, `start`, `end`, …  Used for "element universes" (which genomic regions a guide library targets).

**Open in Excel** for inspection; in pipelines, use `bedtools` or `pyranges`.

---

## GTF — Gene annotation

Tab-delimited gene-coordinate file from GENCODE/Ensembl. Tells you where each gene's exons/transcripts live in the genome. We use a single canonical IGVF GTF across the project; you'll rarely need to crack it open.

**Open in Python:** `pyranges.read_gtf("file.gtf.gz")` or `pandas.read_csv` with `sep="\t"` and comment-line filtering.

---

## Where files live

### IGVF data portal
- URL: <https://data.igvf.org/>
- Each file has a stable accession like `IGVFFI3617IJOW`.
- Each analysis set (a grouped collection of files belonging to one experiment) has an accession like `IGVFDS6332VCTO`.
- To download, click the file's page → "Download" button. Public files are open; some pre-release files require an IGVF consortium account.

### GCS (Google Cloud Storage)
- Bucket path: `gs://igvf-pertub-seq-pipeline-data/<dataset>/<date>/outs/<run>/`.
- For Nextflow CRISPR-pipeline outputs (intermediate + final).
- Download: `gsutil -m cp -r gs://… .` (you'll need a Google account on the project).

### Synapse
- URL: <https://www.synapse.org/>
- For project-level mirroring of finalized outputs (the place to look once everything has been packaged).
- Each file/folder has an ID like `syn74834952`.
- Download: log in via browser, then click "Download Options," or use `synapseclient` in Python:
  ```python
  import synapseclient
  syn = synapseclient.Synapse()
  syn.login(authToken="<your-token>", silent=True)
  path = syn.get("syn74834952").path
  ```

### HPC (local lab compute)
- For active analyses by the computational team. You normally won't need to access this — outputs you'd care about are mirrored to Synapse.

---

## Sizes you'll see

| Format | Typical size | Why |
|---|---|---|
| Per-target TSV (energy distance, regulators-per-program) | KB – MB | One row per TF, ~20–50 columns. |
| Per-perturbation × per-gene TSV (perturbo trans) | 1–10 GB | Millions of rows (n_perturbations × n_genes); reasonable to filter before downloading. |
| `inference_mudata.h5mu` | 10–30 GB | All cells × all genes + guide matrix + dim-reductions. |
| Pipeline dashboard | 30–60 GB | Includes raw count matrices + figures. |
| cNMF bundle (selected k) | 5–7 GB | Integrated MuData + loadings + cell usages. |

If a file isn't loading on your laptop, check the size first — it might be a server-side file you should subset before downloading.
