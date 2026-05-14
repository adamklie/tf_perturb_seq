# Understanding cNMF outputs

cNMF (consensus non-negative matrix factorization) is the project's tool for discovering **gene programs** — groups of genes that vary together across cells. The outputs let you ask: *"What biological modules are active in my cell type, and which TFs disrupted which modules?"* Companion: [`CRISPR.md`](CRISPR.md), [`ENERGY_DISTANCE.md`](ENERGY_DISTANCE.md).

For the technical pipeline reference, see [`../../../analysis/cnmf/cNMF_OUTPUTS.md`](../../../analysis/cnmf/cNMF_OUTPUTS.md). This doc is the **interpretation** layer.

---

## What's a gene program?

Imagine sorting all of the genes in your dataset into 25–250 **groups** by which-genes-go-up-and-down-together-across-cells. Each group is a "program." A program is defined by:

- **Loading scores** per gene: how strongly each gene belongs to that program (some genes are "core" to a program; others are weak members).
- **Usage scores** per cell: how active that program is in each cell.

Categories of programs you might find (cell-type-dependent — vet before assuming):

- A cell-cycle / proliferation program.
- A mitochondrial / oxidative-phosphorylation program.
- A ribosome / translation program.
- A stress-response or interferon-response program.
- A lineage-identity program tied to the cell type at hand.
- A maturation or differentiation-state program.

For Perturb-seq, programs become a powerful intermediate layer: instead of asking "which of 30,000 genes did my TF move," you can ask "which of 50 programs did my TF move." Much more interpretable.

---

## The k-selection decision (and why it matters)

cNMF doesn't know how many programs ("k") to fit. You run it across a sweep (typical: k = 30, 50, 60, 80, 100, 200, 250, 300) and pick the value that balances stability (programs are reproducible across runs) with resolution (more programs = finer biology). The outputs include a `k_selection.png` plot showing stability vs. reconstruction error across the sweep, plus a `README.txt` explaining why a particular k was chosen.

You'll see the selected-k subdirectory prominently and a "sweep-as-provenance" record of every k value (so you can revisit the decision without re-running cNMF).

---

## Files you'll see

A complete cNMF run (selected-k plus sweep-as-provenance) lives under `datasets/<dataset>/cnmf/<run_name>/`.

### Selected-k files (use these for analysis)

| File | Contents |
|---|---|
| `Inference/cNMF_<k>_2_0.h5mu` | The integrated MuData: original cells × genes + a new "cNMF" modality of cells × programs. |
| `Inference/Inference.gene_spectra_score.k_<k>.dt_2_0.txt` | **Gene loadings**: genes × programs, z-scored. The headline "which genes are in which program" table. |
| `Inference/Inference.gene_spectra_tpm.k_<k>.dt_2_0.txt` | Same loadings, TPM-normalized instead of z-scored (useful for ranking genes within a program). |
| `Inference/Inference.usages.k_<k>.dt_2_0.consensus.txt` | **Cell usages**: cells × programs. How much each cell uses each program. |
| `Evaluation/<k>_2_0/<k>_perturbation_association_results_all.txt` | **Regulators per program**: which TF perturbations significantly shifted each program's usage. Columns: target, program, log2FC of program usage, q-value. |
| `Evaluation/<k>_2_0/<k>_geneset_enrichment.txt` | MSigDB enrichments per program. |
| `Evaluation/<k>_2_0/<k>_GO_term_enrichment.txt` | GO term enrichments per program. |

### Sweep-as-provenance files (keep, but you won't normally consult them)

- `Inference/Inference.k_selection.png` plus `Inference.k_selection_stats.df.npz` — the k-selection figure and stats.
- `Inference/Inference.clustering.k_<K>.dt_2_0.png` for every K — visual stability snapshots.
- A `README.txt` documenting the selected k and rationale.

---

## How to use this

### "What programs are active in my cells?"

```python
import pandas as pd
usages = pd.read_csv("Inference.usages.k_50.dt_2_0.consensus.txt", sep="\t", index_col=0)
print(usages.shape)              # (n_cells, n_programs)
print(usages.head())

# Cells with strongest activation of program 7
top_cells_p7 = usages.nlargest(20, "Usage_7")
```

### "What does a given program mean biologically?"

```python
spectra = pd.read_csv("Inference.gene_spectra_score.k_50.dt_2_0.txt", sep="\t", index_col=0)
print(spectra["Program_7"].sort_values(ascending=False).head(20))   # top-loaded genes

# And the enrichment table:
enrich = pd.read_csv("50_geneset_enrichment.txt", sep="\t")
print(enrich[enrich["program"] == "Program_7"].head(10))
```

### "Which TFs perturbed a given program?"

```python
regs = pd.read_csv("50_perturbation_association_results_all.txt", sep="\t")
my_p7 = regs[(regs["program"] == "Program_7") & (regs["qval"] < 0.05)]
print(my_p7.sort_values("log2fc", ascending=True).head(10))   # TFs that knocked the program down
```

### "Compare programs across two datasets"

For cross-lineage analysis, take the `gene_spectra_score` from each dataset and compute the all-vs-all cosine similarity between program loadings. Programs with high similarity across lineages are candidate **lineage-shared** modules; programs with no similar partner are candidate **lineage-specific** modules. The actual thresholds and conclusions are for the analyst to vet against the underlying biology.

```python
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity

A = pd.read_csv("<dataset_A>/gene_spectra_score.k_50.dt_2_0.txt", sep="\t", index_col=0)
B = pd.read_csv("<dataset_B>/gene_spectra_score.k_50.dt_2_0.txt", sep="\t", index_col=0)

shared = A.index.intersection(B.index)
sim = cosine_similarity(A.loc[shared].T, B.loc[shared].T)
```

---

## What to look for

These are methodological pointers — what the outputs *let you ask*, not claims about what you'll find. Vet every interpretation against the underlying gene lists and enrichments.

1. **"Core" recurrent programs** (cell cycle, mitochondria, ribosome, stress-response are common candidates across cell types). Annotating and naming these is a useful first deliverable.
2. **Lineage-identity programs** — programs whose top-loaded genes are cell-type marker genes. The TFs that regulate these programs are candidates for lineage-relevant regulators.
3. **Programs disrupted by a TF of interest** — pull the regulators-per-program file, sort by significance, see which programs the TF moved.
4. **Cross-lineage shared programs** — programs that match across cell types (high cosine similarity) and have overlapping regulator sets are candidates for general regulatory modules.
5. **TFs that disrupt multiple programs** — candidate pleiotropic regulators in the perturbation set.

---

## Issues

- **k is a judgment call.** A different selected k will give you partially different programs. The sweep-as-provenance files let you revisit the decision; if a program your TF disrupted only shows up at a different k, it's worth re-running with that k.
- **Programs are not pathways.** A program enriched for "OXPHOS" doesn't mean only OXPHOS genes are in it — it means OXPHOS is the dominant signal in the enrichment test. Always look at the top-loaded genes, not just the enrichment terms.
- **TPM vs z-scored loadings** answer different questions. Use z-scored for cross-program comparisons within one dataset; use TPM for "what are the actual most-abundant genes in this program."

---

## Who to ask

- **cNMF method / k-selection** — the cNMF / PerturbNMF maintainers (Engreitz Lab); also see [`../../../analysis/cnmf/PerturbNMF.md`](../../../analysis/cnmf/PerturbNMF.md).
- **Program annotation / biology** — meeting-room decision; check the run's `README.txt` for the k-selection meeting notes.
- **Project-specific context** — file a GitHub issue tagging the dataset's lead.
