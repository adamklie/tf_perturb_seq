"""WG2 examples — gene-program cross-lineage similarity + regulators per program.

**Most sections are GATED** on production cNMF runs landing. Today (2026-05-12):
HTv2 testbed verified; Huangfu DE + ESC production launches pending; Hon CM held
on collaborator deliverables from Weizhou. See `working_groups/wg2_gene_programs/README.md` for status.

What the sections will do once data is mirrored:

  §1  load each dataset's `gene_spectra_score.k_<sel>.dt_2_0.txt`
  §2  build a cross-dataset cosine-similarity matrix of program loadings
  §3  classify programs (lineage_shared / lineage_specific / cell_state) from §2
  §4  extract top-N genes per program (table)
  §5  load `<sel>_perturbation_association_results_*.txt` per dataset → regulators

Synapse paths for the mirrored cNMF folders are recorded by the mirror scripts
under `data/scripts/`; each dataset's Synapse folder also appears in its
`data/<dataset>/README.md`.
"""
from __future__ import annotations

from pathlib import Path

JAMBOREE_ROOT = Path(__file__).resolve().parents[3]

DATASETS = {
    "HonCM": "Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq",
    "HuangfuDE": "Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq",
    "HuangfuESC": "Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq",
    "GersbachHep": "Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq",
    "EngreitzEndo": "Engreitz_WTC11-endothelial-cells_TF-Perturb-seq",
}


# §1 — Load gene_spectra_score per dataset (GATED) ----------------------------
print("\n§1 Load gene_spectra_score (GATED — needs cNMF outputs mirrored)")
print("    Once mirrored under datasets/<id>/cnmf/<run>/, do something like:")
print("""
    import pandas as pd
    from pathlib import Path
    JAMBOREE_ROOT = Path(__file__).resolve().parents[3]
    selected_k = {"HuangfuDE": 25, "HuangfuESC": 25, "HonCM": None, ...}
    spectra = {}
    for short, k in selected_k.items():
        if k is None:
            continue
        full = DATASETS[short]
        run_dir = next((JAMBOREE_ROOT / "datasets" / full / "cnmf").glob("*"))
        path = run_dir / "Result" / f"gene_spectra_score.k_{k}.dt_2_0.txt"
        spectra[short] = pd.read_csv(path, sep="\\t", index_col=0)
""")


# §2 — Cross-dataset cosine similarity (GATED) -------------------------------
print("\n§2 Cosine similarity of program loadings across datasets")
print("""
    from sklearn.metrics.pairwise import cosine_similarity
    import numpy as np
    blocks = []
    labels = []
    for short, df in spectra.items():
        for prog in df.columns:
            blocks.append(df[prog].values)
            labels.append(f"{short}::{prog}")
    M = cosine_similarity(np.stack(blocks))
    sim = pd.DataFrame(M, index=labels, columns=labels)
    sim.to_csv("program_similarity_matrix.tsv", sep="\\t")
""")


# §3 — Classify programs (GATED) ---------------------------------------------
print("\n§3 Lineage_shared / lineage_specific / cell_state classification")
print("""
    Heuristic: a program with max off-diagonal cosine ≥ 0.5 against ≥1 other
    lineage is lineage_shared; ≥ 0.5 only within-dataset is cell_state;
    < 0.5 across-board is lineage_specific. Tune thresholds with the group.
""")


# §4 — Top-N genes per program (GATED) ---------------------------------------
print("\n§4 Top-N genes per program — per-dataset companion artifact")
print("""
    rows = []
    for short, df in spectra.items():
        for prog in df.columns:
            top = df[prog].sort_values(ascending=False).head(20)
            for rank, (gene, score) in enumerate(top.items(), 1):
                rows.append({"dataset": short, "program": prog, "rank": rank,
                             "gene": gene, "score": score})
    pd.DataFrame(rows).to_csv("wg2_top_genes_per_program.tsv", sep="\\t", index=False)
""")


# §5 — Regulators per program (GATED) ----------------------------------------
print("\n§5 Per-dataset regulators per program (TF perturbations associated)")
print("""
    # Source: <sel>_perturbation_association_results_*.txt under
    # datasets/<id>/cnmf/<run>/Eval/<sel>_<dt>/
    # See schemas/cnmf.json for the full file inventory.
""")

print("\nAll sections gated. Once a cNMF mirror lands, replace the print blocks "
      "with the indented code above (also queued at examples/build_wg2_*.py).")

