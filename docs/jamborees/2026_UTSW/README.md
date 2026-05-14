# 2026 UTSW Jamboree

CRISPRi Perturb-seq of ~2,000 TFs across five human cell lineages. We've packaged the data, mirrored everything to Synapse, written a participant guide, and set up a GitHub + Synapse folder per working group. This page tells you where to go for what you need and how a typical WG runs at the jamboree.
---

## Logistics

| | |
|---|---|
| **Dates** | 2026-05-13 to 2026-05-16 |
| **Venue** | *[FILL IN]* |
| **Schedule** | *[FILL IN — link Google Doc]* |
| **Hybrid / Zoom** | *[FILL IN]* |
| **Slack** | `#tf-perturb-enhancer-tiger-team` |
| **GitHub** | [adamklie/tf_perturb_seq](https://github.com/adamklie/tf_perturb_seq) |
| **Issues / blockers** | [GitHub Issues](https://github.com/adamklie/tf_perturb_seq/issues) |
| **Synapse project root** | [`syn64423137/2026_UTSW/`](https://www.synapse.org/Synapse:syn64423137) |
| **Day-of contact** | Adam Klie (`aklie@ucsd.edu`) |

## How the jamboree works

Each working group runs through the same three steps:

1. **Brainstorm** — sketch the analyses the group wants to accomplish: example figures, summary tables, pseudocode. Capture this in a notebook, doc, or notes in your WG folder.
2. **Execute** — run the analyses; commit code to your WG's GitHub folder under [`working_groups/wg<N>_*/`](working_groups/). Notebooks, scripts, and supporting docs all live there.
3. **Share** — upload reusable objects (figures, tables, intermediate data) to your WG's Synapse folder, and record the syn ID in your WG README so the next person can find it. The [Synapse IDs at a glance](#synapse-ids-at-a-glance) table below has all the starting points.

## I want to…

| …do this | …go here |
|---|---|
| Understand what the outputs mean | [`guide/`](guide/) — glossary, FAQ, per-output interpretation |
| See what's available per dataset | [`data/README.md`](data/README.md) |
| Pick a working group | [`working_groups/README.md`](working_groups/README.md) |
| Set up my environment + Synapse | [Set up your environment](#set-up-your-environment) below |
| Drill into a single dataset | [`data/<dataset>/`](data/) (per-dataset READMEs) |
| Find DACC-spec deliverables | [`reference/tf_universe.tsv`](reference/tf_universe.tsv) + [`reference/element_universe.bed`](reference/element_universe.bed) for the library; per-dataset gene_universe files live on Synapse |
| Look up a Synapse ID | [Synapse IDs at a glance](#synapse-ids-at-a-glance) below |
| Report a blocker / open question | [GitHub Issues](https://github.com/adamklie/tf_perturb_seq/issues) |

## Datasets at a glance

| Dataset | CRISPR pipeline | cNMF | Energy distance |
|---|:---:|:---:|:---:|
| Hon WTC11 Cardiomyocyte | ready | - | ready |
| Huangfu HUES8 Definitive Endoderm | ready | ready | ready |
| Huangfu HUES8 Embryonic Stem Cell | ready | ready | ready |
| Gersbach WTC11 Hepatocyte | ready | - | ready |
| Engreitz WTC11 Endothelial | - | - | - |

Full per-dataset cards: [`data/README.md`](data/README.md) (and each dataset's own `README.md`).

## Working groups

| WG | Folder | Topic | Lead questions |
|---|---|---|---|
| 1 | [`wg1_data_qc/`](working_groups/wg1_data_qc/) | Topic 1 / Fig 1 | Guide detection & repression; transcriptome-wide significance; cross-lineage shared TFs |
| 2 | [`wg2_gene_programs/`](working_groups/wg2_gene_programs/) | Topic 2.1 / Fig 2 | Cross-lineage program similarity; regulators per program |
| 3 | [`wg3_disease_gwas/`](working_groups/wg3_disease_gwas/) | Topic 2.2 / Fig 2 | Disease-gene TFs; convergent vs divergent activity |
| 4 | [`wg4_grn_inference/`](working_groups/wg4_grn_inference/) | Topic 2.3 / Fig 3 | TF→gene edges; cross-lineage rewiring |
| 5 | [`wg5_tf_family_case_studies/`](working_groups/wg5_tf_family_case_studies/) | Topic 2.4 / Fig 4 | TF family activity scorecard; per-family deep-dives |
| 6 (opt) | [`wg6_predictive_modeling/`](working_groups/wg6_predictive_modeling/) | Topic 3 / Fig 5 | Predictive model task spec |

Each WG folder has its own README with scope, lead questions, dataset status, and the WG's Synapse output folder.

## Synapse IDs at a glance

Everything lives under the [`syn64423137`](https://www.synapse.org/Synapse:syn64423137) project root.

**Project + jamboree folders**

| Object | Syn ID |
|---|---|
| Project root (`tf_perturb_seq`) | [`syn64423137`](https://www.synapse.org/Synapse:syn64423137) |
| 2026 UTSW jamboree folder | [`syn74834225`](https://www.synapse.org/Synapse:syn74834225) |
| Working groups root | [`syn74954078`](https://www.synapse.org/Synapse:syn74954078) |

**Working-group output folders** (upload your WG's deliverables here)

| Working group | Syn ID |
|---|---|
| WG1 — data QC | [`syn74954079`](https://www.synapse.org/Synapse:syn74954079) |
| WG2 — gene programs | [`syn74954081`](https://www.synapse.org/Synapse:syn74954081) |
| WG3 — disease and GWAS | [`syn74954083`](https://www.synapse.org/Synapse:syn74954083) |
| WG4 — GRN inference | [`syn74954084`](https://www.synapse.org/Synapse:syn74954084) |
| WG5 — TF family case studies | [`syn74954085`](https://www.synapse.org/Synapse:syn74954085) |
| WG6 — predictive modeling | [`syn74954086`](https://www.synapse.org/Synapse:syn74954086) |

**Reference data**

| Object | Syn ID |
|---|---|
| TF metadata (`tf_metadata.tsv`) | [`syn74834227`](https://www.synapse.org/Synapse:syn74834227) |
| Experimental metadata (`experimental_metadata.tsv`) | [`syn74834309`](https://www.synapse.org/Synapse:syn74834309) |
| Guide library (`IGVFFI8270UPKB.csv.gz`) | [`syn74834519`](https://www.synapse.org/Synapse:syn74834519) |
TODO

**Per-dataset CRISPR pipeline outputs**

| Dataset | Syn ID |
|---|---|
| Huangfu HUES8 Definitive Endoderm | [`syn74834952`](https://www.synapse.org/Synapse:syn74834952) |
| Huangfu HUES8 Embryonic Stem Cell | [`syn74835010`](https://www.synapse.org/Synapse:syn74835010) |
| Hon WTC11 Cardiomyocyte | *[FILL IN]* |
| Gersbach WTC11 Hepatocyte | *[FILL IN]* |
| Engreitz WTC11 Endothelial | *[FILL IN]* |

Per-dataset cards (with cNMF + energy-distance Synapse pointers as they land) are under [`data/<dataset>/README.md`](data/).

---

## Set up your environment

### 0. Install uv

[`uv`](https://github.com/astral-sh/uv) is the Python package manager this project uses. Install it with `pip`, `brew`, or the official installer (see the [uv install docs](https://docs.astral.sh/uv/getting-started/installation/)):

```bash
pip install uv
```

### 1. Clone the repo + Python env

```bash
git clone https://github.com/adamklie/tf_perturb_seq.git
cd tf_perturb_seq
uv venv && source .venv/bin/activate
uv pip install mudata anndata scanpy pandas synapseclient
```

### 2. Synapse access

```bash
# Get a token at https://www.synapse.org/Profile:v/settings ("Personal Access Tokens")
export SYNAPSE_AUTH_TOKEN="<your-token>"
# Persist in ~/.zshrc or ~/.bashrc
```

```python
import os, synapseclient
syn = synapseclient.Synapse()
syn.login(authToken=os.environ["SYNAPSE_AUTH_TOKEN"], silent=True)
```

### 3. Load reference data

```python
import pandas as pd

# TF metadata (1,983 unique TF target genes × 16 cols)
tf_meta = pd.read_csv(syn.get("syn74834227").path, sep="\t")

# Experimental metadata (5 datasets × 26 cols)
exp_meta = pd.read_csv(syn.get("syn74834309").path, sep="\t")

# Guide library (14,150 guides × 18 cols, IGVF pools A-D; gzipped TSV despite .csv.gz)
guides = pd.read_csv(syn.get("syn74834519").path, sep="\t")
```

### 4. Load a dataset's CRISPR pipeline output

```python
import mudata as md
# Huangfu DE example — Synapse syn74834952; find inference_mudata.h5mu under pipeline_dashboard/
inf_id = "syn-id-of-inference_mudata.h5mu"  # syn ID of the .h5mu file inside pipeline_dashboard/
mdata = md.read_h5mu(syn.get(inf_id).path)
print(mdata.mod["gene"].shape)              # n_cells × n_genes
print(mdata.mod["guide"].shape)             # n_cells × n_guides
```

See [Synapse IDs at a glance](#synapse-ids-at-a-glance) for the per-dataset Synapse IDs and what's landed so far. Deeper walkthroughs: [`guide/CRISPR.md`](guide/CRISPR.md), [`guide/ENERGY_DISTANCE.md`](guide/ENERGY_DISTANCE.md), [`guide/CNMF.md`](guide/CNMF.md).
