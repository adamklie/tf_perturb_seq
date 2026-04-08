# Team

## GPS Tiger Team

**Meeting:** Wednesdays 10a PST / 12p CST / 1p EST
**Slack:** #tf-perturb-enhancer-tiger-team

## Analysis leads and assignments

| Person | Role | Sub-aim | Assigned Datasets |
|--------|------|---------|-------------------|
| Adam Klie | Project coordinator, repo maintainer, QC pipeline | 1 (Integrative) | Huangfu_HUES8-embryonic-stemcell, Gersbach hepatocyte |
| Sara Geraghty | Project coordinator, data processing, QC, energy distance analysis | 1 (Integrative) | Gersbach hepatocyte |
| Sid (Sidharth) Raghavan | Integrative analysis, jamboree coordination | 1 (Integrative) | — |
| Chikara Takeuchi | Energy distance pipeline, gene programs | 1 (Integrative) | Hon neuron |
| Weizhou Qian | PerturbNet, perturbation prediction | 2 (Prediction) | Hon cardiomyocyte |
| Gary Yang | gkmSVM, sequence models, multiome | 3 (Sequence models) | — |
| Alexandra Mo | cNMF-GPU implementation | — | — |
| Mpathi Nzima | Guide metadata, portal uploads, pipeline runs |
| Deins Torre |

## DACC
| Ian Whaling |


## CRISPR pipeline team

| Person | Role |
|--------|------|
| Lucas Ferreira | Pipeline developer, TF enrichment automation |
| Logan Blaine | Perturbo confidence/uncertainty estimates |

| Sushama Sivakumar | Validated TF-target pairs curation |

## Data producing labs

| Lab | PI | Cell Models | Benchmark Tech |
|-----|----|-----------| --------------|
| Hon | Gary Hon | WTC11 iPSC, cardiomyocytes, neurons, trophoblasts | 10x 5' HT v2 w/ HTO |
| Huangfu | Danwei Huangfu | HUES8 hESC, definitive endoderm | 10x 3' v3 |
| Engreitz | Jesse Engreitz | WTC11 iPSC, endothelial | CC Perturb-seq |
| Gersbach | Charles Gersbach | WTC11 iPSC, hepatocytes | 10x GEM-X 3', 10x HT v2 |

## Analysis sub-aims

### Sub-aim 1: Basic Integrative Analysis (Adam, Sid, Chikara)
- Cross-lineage energy distance comparison
- DE gene network construction and comparison
- Identify conserved vs. lineage-specific TF effects

### Sub-aim 2: Predicting Perturbation Effects (Weizhou)
- Enhance PerturbNet for TF Perturb-seq
- Counterfactual predictions for unseen perturbations
- Output: H5AD with perturbation predictions

### Sub-aim 3: Integrating Sequence Models (Revant, Gary)
- ChromBPNet on multiome data
- E2G linking (cis-regulatory modules)
- Integrate with Perturb-seq GRNs
- Variant interpretation and in silico perturbation
