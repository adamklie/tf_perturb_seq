"""
CM_ED_SC_overlap.py
-------------------
Cross-dataset comparison of expressed genes and regulators (guide-target genes)
across 3 IGVF perturb-seq datasets:

    CM = cardiomyocyte                  (IGVF_Hon_CM)
    DE = definitive endoderm            (IGVF_Huangfu_definitive-endoderm)
    SC = embryonic stem cell            (IGVF_Huangfu_embryonic-stemcell)

Answers six questions in one pass:
    1) common expressed genes
    2) tissue-specific expressed genes
    3) common regulators
    4) tissue-specific regulators
    5) common expressed regulators           (regulators that are also expressed in their tissue)
    6) tissue-specific expressed regulators

Outputs (PDF + PNG Venn diagrams, plus 2 TSVs):
    Meta_program/Result/CM_ED_SC_analysis/
"""
# %% Imports
from __future__ import annotations
from pathlib import Path

import re

import muon as mu
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3


# %% Config
DATASETS: dict[str, str] = {
    "CM": "/oak/stanford/groups/engreitz/Users/ymo/Project/IGVF_Hon_CM/Data/inference_mudata.h5mu",
    "DE": "/oak/stanford/groups/engreitz/Users/ymo/Project/IGVF_Huangfu_definitive-endoderm/Result/Adam_run/Data/inference_mudata.h5mu",
    "SC": "/oak/stanford/groups/engreitz/Users/ymo/Project/IGVF_Huangfu_embryonic-stemcell/Result/Adam_run/Data/inference_mudata.h5mu",
}

OUT_DIR = Path(
    "/oak/stanford/groups/engreitz/Users/ymo/Testing_Scripts/Meta_program/Result/CM_ED_SC_analysis"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)


# %% Load only var tables from each h5mu (backed mode -> no X load)
def load_vars(path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    md = mu.read(path, backed="r")
    gene_var = md["gene"].var.copy()
    guide_var = md["guide"].var.copy()
    md.file.close()
    return gene_var, guide_var


gene_vars: dict[str, pd.DataFrame] = {}
guide_vars: dict[str, pd.DataFrame] = {}
for name, path in DATASETS.items():
    print(f"[load] {name}: {path}")
    gene_vars[name], guide_vars[name] = load_vars(path)
    print(f"       gene.var: {gene_vars[name].shape}   guide.var: {guide_vars[name].shape}")


# %% Build the six set universes
# Expressed genes — Ensembl IDs are directly comparable across datasets
expressed_ensembl: dict[str, set[str]] = {
    k: set(v.index.astype(str)) for k, v in gene_vars.items()
}
# Expressed genes in symbol space — needed to intersect with regulator symbols
expressed_symbol: dict[str, set[str]] = {
    k: set(v["symbol"].dropna().astype(str)) for k, v in gene_vars.items()
}

# Regulators — targeting guides only, deduplicated by gene_name (gene symbol).
# (intended_target_name in these files is Ensembl ID, gene_name is the symbol.)
def regulator_set(guide_var: pd.DataFrame) -> set[str]:
    mask = guide_var["targeting"].astype(str).str.lower().isin(["true", "1"])
    vals = guide_var.loc[mask, "gene_name"].dropna().astype(str)
    # Drop empties and any leftover control labels
    bad = {"", "nan", "none"}
    return {
        v for v in vals.unique()
        if v.lower() not in bad and not v.lower().startswith(("non-targeting", "safe-harbor", "safe_harbor"))
    }

regulators: dict[str, set[str]] = {k: regulator_set(v) for k, v in guide_vars.items()}

# Expressed regulators — regulator symbols that are also expressed in that tissue
expr_regulators: dict[str, set[str]] = {
    k: regulators[k] & expressed_symbol[k] for k in DATASETS
}


# %% Venn diagrams
def plot_venn3(sets: dict[str, set[str]], title: str, stem: Path) -> None:
    fig, ax = plt.subplots(figsize=(6, 6))
    venn3(
        subsets=[sets["CM"], sets["DE"], sets["SC"]],
        set_labels=("CM", "DE", "SC"),
        ax=ax,
    )
    ax.set_title(title)
    for ext in (".pdf", ".png"):
        fig.savefig(stem.with_suffix(ext), dpi=200, bbox_inches="tight")
    plt.close(fig)


plot_venn3(expressed_ensembl, "Expressed genes (Ensembl IDs)",        OUT_DIR / "venn_expressed_genes")
plot_venn3(regulators,        "Regulators (targeting guide symbols)", OUT_DIR / "venn_regulators")
plot_venn3(expr_regulators,   "Expressed regulators (symbol)",        OUT_DIR / "venn_expressed_regulators")


# %% Summary tables
def venn_regions(d: dict[str, set[str]]) -> dict[str, set[str]]:
    CM, DE, SC = d["CM"], d["DE"], d["SC"]
    return {
        "common_all_three": CM & DE & SC,
        "CM_only":          CM - DE - SC,
        "DE_only":          DE - CM - SC,
        "SC_only":          SC - CM - DE,
        "CM_DE_not_SC":     (CM & DE) - SC,
        "CM_SC_not_DE":     (CM & SC) - DE,
        "DE_SC_not_CM":     (DE & SC) - CM,
    }


universes: dict[str, dict[str, set[str]]] = {
    "expressed_genes":       expressed_ensembl,
    "regulators":            regulators,
    "expressed_regulators":  expr_regulators,
}

# Combined Ensembl -> symbol lookup (used to add symbols for expressed_genes regions)
ens_to_sym: dict[str, str] = {}
for k, gv in gene_vars.items():
    sym_col = gv["symbol"].dropna().astype(str)
    for ens, sym in sym_col.items():
        ens_to_sym.setdefault(str(ens), sym)


# Per-region gene lists
REGIONS_DIR = OUT_DIR / "regions"
REGIONS_DIR.mkdir(exist_ok=True)

rows: list[dict[str, object]] = []
for cat, sets in universes.items():
    for region, items in venn_regions(sets).items():
        rows.append({"category": cat, "region": region, "n": len(items)})
        items_sorted = sorted(items)
        out_tsv = REGIONS_DIR / f"{cat}__{region}.tsv"
        if cat == "expressed_genes":
            # universe is Ensembl IDs -> also write gene symbol
            df = pd.DataFrame({
                "ensembl_id": items_sorted,
                "symbol":     [ens_to_sym.get(e, "") for e in items_sorted],
            })
        else:
            # regulators / expressed_regulators are already gene symbols
            df = pd.DataFrame({"symbol": items_sorted})
        df.to_csv(out_tsv, sep="\t", index=False)

summary_long = pd.DataFrame(rows)
summary_long.to_csv(OUT_DIR / "summary_counts.tsv", sep="\t", index=False)

# Wide per-dataset set-size table
size_df = pd.DataFrame(
    {
        "n_expressed_genes":        {k: len(v) for k, v in expressed_ensembl.items()},
        "n_expressed_symbols":      {k: len(v) for k, v in expressed_symbol.items()},
        "n_regulators":             {k: len(v) for k, v in regulators.items()},
        "n_expressed_regulators":   {k: len(v) for k, v in expr_regulators.items()},
    }
)
size_df.index.name = "dataset"
size_df.to_csv(OUT_DIR / "set_sizes.tsv", sep="\t")


# %% Bundle everything into one Excel workbook
SHEET_PREFIX = {
    "expressed_genes":      "eg",
    "regulators":           "reg",
    "expressed_regulators": "er",
}
xlsx_path = OUT_DIR / "CM_ED_SC_overlap.xlsx"
with pd.ExcelWriter(xlsx_path, engine="openpyxl") as xw:
    summary_long.to_excel(xw, sheet_name="summary_counts", index=False)
    size_df.to_excel(xw, sheet_name="set_sizes")
    for cat, sets in universes.items():
        for region, items in venn_regions(sets).items():
            items_sorted = sorted(items)
            if cat == "expressed_genes":
                df = pd.DataFrame({
                    "ensembl_id": items_sorted,
                    "symbol":     [ens_to_sym.get(e, "") for e in items_sorted],
                })
            else:
                df = pd.DataFrame({"symbol": items_sorted})
            sheet = f"{SHEET_PREFIX[cat]}__{region}"[:31]  # Excel sheet name limit
            df.to_excel(xw, sheet_name=sheet, index=False)


# %% Print headline summary
print("\n=== Set sizes per dataset ===")
print(size_df.to_string())

print("\n=== Triple-intersection sizes ===")
for cat, sets in universes.items():
    common = sets["CM"] & sets["DE"] & sets["SC"]
    print(f"  {cat:>22s}: |CM ∩ DE ∩ SC| = {len(common):,}")

print(f"\nOutputs written to: {OUT_DIR}")
print(f"  Excel workbook: {xlsx_path.name}")
print(f"  Per-region TSVs: {REGIONS_DIR.relative_to(OUT_DIR)}/")

# %%
