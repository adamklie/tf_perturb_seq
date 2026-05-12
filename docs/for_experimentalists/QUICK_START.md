# Quick start — from "I have a Synapse link" to "I can see results"

This is the 10-minute on-ramp. After this you can poke at a result file without writing serious code. For deeper analysis you'll want to read the per-output `UNDERSTAND_*.md` files.

---

## Step 1 — Get access

### Synapse
1. Make an account at <https://www.synapse.org/> (it's free).
2. Sign the IGVF data-use terms if asked.
3. Ask whoever pointed you at the data for the right Synapse parent folder (typically `syn64423137` for our project).
4. (Optional, for command-line) Get a personal access token: Profile → Settings → "Personal Access Tokens" → "Create New Token."

### IGVF data portal
1. Public files: just navigate to <https://data.igvf.org/> and search by accession.
2. Pre-release files (consortium-only): get an account via the IGVF coordinating center.

---

## Step 2 — Peek at one TSV without code

Most cross-dataset summaries live in `reference/` folders and are small TSVs.

1. Click the Synapse link, hit "Download Options" → "Download File."
2. Open the `.tsv` (or `.tsv.gz`) in Excel/Numbers/Google Sheets — it will ask for a delimiter; choose **Tab**.
3. Skim the first few rows; column headers are usually self-describing.
4. To narrow down: filter / sort by a column you care about (e.g. `sig_distance_gt_NC_max == True`).

That's it — for many "give me a number" questions, this is enough.

---

## Step 3 — Look at a dataset's `inference_mudata.h5mu`

The MuData is the canonical CRISPR-pipeline output. It contains every cell, every gene, every guide.

### Option A — In Python (recommended)

Make a virtual environment first:

```bash
uv venv && source .venv/bin/activate
uv pip install mudata anndata scanpy pandas synapseclient
```

Then:

```python
import os
import synapseclient
import mudata

syn = synapseclient.Synapse()
syn.login(authToken=os.environ["SYNAPSE_AUTH_TOKEN"], silent=True)

# Replace with the actual file accession (find it under the dataset's pipeline_dashboard/ on Synapse)
mu_path = syn.get("syn-id-of-inference_mudata.h5mu").path
mdata = mudata.read_h5mu(mu_path)
print(mdata)
print(mdata.mod["gene"].shape)              # n_cells × n_genes
print(mdata.mod["gene"].obs.head())         # per-cell metadata (guide id, etc.)
print(mdata.mod["gene"].var.head())         # per-gene metadata

# Cells assigned to a specific TF guide
sub = mdata.mod["gene"][mdata.mod["gene"].obs["intended_target_name"] == "SOX17"]
print(sub.shape)
```

### Option B — Don't touch the MuData

For most "see the result" questions, the **perturbo TSVs** alongside the MuData are what you actually want (one row per perturbation × gene). See [`UNDERSTAND_CRISPR_OUTPUTS.md`](UNDERSTAND_CRISPR_OUTPUTS.md).

---

## Step 4 — Look at the QC dashboard

Every CRISPR-pipeline bundle has a `pipeline_dashboard/dashboard.html`. Download it (and the `figures/` folder alongside it) and open in any browser. You'll see:

- Per-measurement-set cell counts.
- Guide-recovery rates per TF.
- Knockdown-efficiency distributions.
- UMAPs of cells colored by lineage / cell state / guide identity.
- Pipeline run summary (which steps ran, software versions).

This is the fastest path to "does this dataset look OK." See [`UNDERSTAND_CRISPR_OUTPUTS.md`](UNDERSTAND_CRISPR_OUTPUTS.md) for what to look for.

---

## Step 5 — Look at "which TFs did something"

If the dataset has **energy-distance** results, that's the headline statistic. Find `pval_edist_full.csv` under the dataset's `energy_distance/` folder on Synapse.

```python
import pandas as pd
ed = pd.read_csv(syn.get("syn-id-of-pval_edist_full.csv").path, index_col=0)

# Show the most-perturbed TFs (largest distance_mean)
top = ed[ed["type"] == "targeting"].sort_values("distance_mean", ascending=False).head(20)
print(top[["distance_mean", "pval_mean"]])
```

See [`UNDERSTAND_ENERGY_DISTANCE.md`](UNDERSTAND_ENERGY_DISTANCE.md) — including the important calibration caveat.

---

## Step 6 — Look at gene programs (cNMF)

cNMF outputs are the place to find "what biological programs did your TF disrupt." See [`UNDERSTAND_CNMF_OUTPUTS.md`](UNDERSTAND_CNMF_OUTPUTS.md).

---

## When something doesn't open

- **"Not enough memory"** — your file is the multi-GB MuData. Filter it server-side (HPC / Colab) or only load the cells/columns you need.
- **"Invalid file format"** — confirm the extension matches the content. `.csv.gz` is gzipped CSV; needs `gzip` decoding (Python/R handle automatically; raw `cat` won't).
- **"Permission denied" on Synapse** — you may not be logged in, or the file is still embargoed. Check with whoever shared the link.

For everything else: open a GitHub issue (link in [`docs/REFERENCES.md`](../REFERENCES.md)) or ping the project Slack channel.

---

## What to read next

- [`UNDERSTAND_CRISPR_OUTPUTS.md`](UNDERSTAND_CRISPR_OUTPUTS.md) — for the CRISPR pipeline outputs.
- [`UNDERSTAND_ENERGY_DISTANCE.md`](UNDERSTAND_ENERGY_DISTANCE.md) — for energy-distance results.
- [`UNDERSTAND_CNMF_OUTPUTS.md`](UNDERSTAND_CNMF_OUTPUTS.md) — for cNMF gene programs.
- [`COMPLETE_DATASET_CONTENTS.md`](COMPLETE_DATASET_CONTENTS.md) — to see whether your dataset is fully baked.
