# 2026 UTSW Jamboree — start here

CRISPRi Perturb-seq of ~2,000 TFs across five human cell lineages. We've packaged the data, mirrored everything to Synapse, written a participant guide, and rolled up cross-dataset summaries per working group. This page tells you where to go for what you need.

> **Online docs**: published participant guide → [tf-perturb-seq.readthedocs.io](https://tf-perturb-seq.readthedocs.io/) *(URL active once the first RTD build completes)*.

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

## I want to…

| …do this | …go here |
|---|---|
| Understand what the outputs mean | [`guide/`](guide/) — glossary, FAQ, per-output interpretation |
| See what's available per dataset | [`data/DATASET_STATUS.md`](data/DATASET_STATUS.md) |
| Pick a working group | [`working_groups/README.md`](working_groups/README.md) |
| See the scientific scope | [`working_groups/TOPICS.md`](working_groups/TOPICS.md) + [`working_groups/WORKING_GROUPS.md`](working_groups/WORKING_GROUPS.md) |
| Set up my environment + Synapse | [Set up your environment](#set-up-your-environment) below |
| Use Claude Code with this repo | [`CLAUDE_QUICKSTART.md`](CLAUDE_QUICKSTART.md) |
| Find Synapse IDs for every output | [`data/synapse_paths.tsv`](data/synapse_paths.tsv) |
| Drill into a single dataset | [`data/<dataset>/`](data/) (per-dataset READMEs + per-analysis subdirs) |
| Find DACC-spec deliverables | [`data/<dataset>/dacc/`](data/) for gene_universe; [`reference/tf_universe.tsv`](reference/tf_universe.tsv) + [`reference/element_universe.bed`](reference/element_universe.bed) for the library |
| Report a blocker / open question | [GitHub Issues](https://github.com/adamklie/tf_perturb_seq/issues) |
| Know who works on what | [`docs/TEAM.md`](../../TEAM.md) |
| See cross-dataset summaries (slide-deck-sized) | [`reference/cross_dataset_pipeline_summary.tsv`](reference/cross_dataset_pipeline_summary.tsv) + [`reference/cross_dataset_edistance_summary.tsv`](reference/cross_dataset_edistance_summary.tsv) |

## Datasets at a glance

| Dataset | CRISPR pipeline | cNMF | Energy distance |
|---|:---:|:---:|:---:|
| Hon WTC11 Cardiomyocyte | ⚠ | ☐ | ✅ |
| Huangfu HUES8 Definitive Endoderm | ✅ | ✅ | ⚠ |
| Huangfu HUES8 Embryonic Stem Cell | ✅ | ✅ | ⚠ |
| Gersbach WTC11 Hepatocyte | ⚠ | ☐ | ☐ |
| Engreitz WTC11 Endothelial | ☐ | ☐ | ☐ |

✅ on Synapse, ready • ⚠ on Synapse with a known caveat • ☐ blocked. Full per-dataset cards: [`data/DATASET_STATUS.md`](data/DATASET_STATUS.md).

## Working groups

| WG | Folder | Topic | Lead questions |
|---|---|---|---|
| 1 | [`wg1_data_qc/`](working_groups/wg1_data_qc/) | Topic 1 / Fig 1 | Guide detection & repression; transcriptome-wide significance; cross-lineage shared TFs |
| 2 | [`wg2_gene_programs/`](working_groups/wg2_gene_programs/) | Topic 2.1 / Fig 2 | Cross-lineage program similarity; regulators per program |
| 3 | [`wg3_disease_gwas/`](working_groups/wg3_disease_gwas/) | Topic 2.2 / Fig 2 | Disease-gene TFs; convergent vs divergent activity |
| 4 | [`wg4_grn_inference/`](working_groups/wg4_grn_inference/) | Topic 2.3 / Fig 3 | TF→gene edges; cross-lineage rewiring |
| 5 | [`wg5_tf_family_case_studies/`](working_groups/wg5_tf_family_case_studies/) | Topic 2.4 / Fig 4 | TF family activity scorecard; per-family deep-dives |
| 6 (opt) | [`wg6_predictive_modeling/`](working_groups/wg6_predictive_modeling/) | Topic 3 / Fig 5 | Predictive model task spec |

Each WG folder has its own README, pre-computed roll-up TSVs, and an `examples.py` showing how to load + filter the data.

---

## Set up your environment

### 1. Synapse access

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

### 2. Clone the repo + Python env

```bash
git clone https://github.com/adamklie/tf_perturb_seq.git
cd tf_perturb_seq
uv venv && source .venv/bin/activate
uv pip install mudata anndata scanpy pandas synapseclient
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
# Huangfu DE example — Synapse syn74834952 (find inference_mudata.h5mu under pipeline_dashboard/)
mu_path = syn.get("syn-id-of-inference_mudata.h5mu").path
mdata = md.read_h5mu(mu_path)
print(mdata.mod["gene"].shape)              # n_cells × n_genes
print(mdata.mod["guide"].shape)             # n_cells × n_guides
```

Deeper walkthroughs: [`guide/UNDERSTAND_CRISPR_OUTPUTS.md`](guide/UNDERSTAND_CRISPR_OUTPUTS.md), [`guide/UNDERSTAND_ENERGY_DISTANCE.md`](guide/UNDERSTAND_ENERGY_DISTANCE.md), [`guide/UNDERSTAND_CNMF_OUTPUTS.md`](guide/UNDERSTAND_CNMF_OUTPUTS.md).

---

## Where to ask for help

- **Slack**: `#tf-perturb-enhancer-tiger-team`
- **GitHub Issues**: <https://github.com/adamklie/tf_perturb_seq/issues>
- **Online docs**: [tf-perturb-seq.readthedocs.io](https://tf-perturb-seq.readthedocs.io/)
- **Who works on what**: [`docs/TEAM.md`](../../TEAM.md)
- **Day-of**: Adam Klie (`aklie@ucsd.edu`)
