"""Summarize the IGVF production guide library (IGVFFI8270UPKB).

Reads the guide table and writes a counts TSV plus a composition figure.
Despite the .csv extension, the file is tab-separated.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

REPO_ROOT = Path(__file__).resolve().parents[6]
GUIDE_TSV = REPO_ROOT / "docs/jamborees/2026_UTSW/reference/IGVFFI8270UPKB.csv"
OUT_DIR = Path(__file__).resolve().parent

sys.path.insert(0, str(REPO_ROOT / "config"))
from loader import load_colors  # noqa: E402

PALETTE = load_colors("guide_library")
TYPE_ORDER = PALETTE["type_order"]
TYPE_PALETTE = PALETTE["type_colors"]
STRAND_PALETTE = PALETTE["strand_colors"]


def load_guides() -> pd.DataFrame:
    df = pd.read_csv(GUIDE_TSV, sep="\t", dtype=str)
    df["type"] = pd.Categorical(df["type"], categories=TYPE_ORDER, ordered=True)
    return df


def write_summary_tsv(df: pd.DataFrame) -> None:
    rows = []
    for guide_type, sub in df.groupby("type", observed=True):
        rows.append(
            {
                "type": guide_type,
                "n_guides": len(sub),
                "n_unique_intended_targets": sub["intended_target_name"].nunique(dropna=True),
                "n_unique_gene_names": sub["gene_name"].replace("", pd.NA).nunique(dropna=True),
            }
        )
    summary = pd.DataFrame(rows)
    summary.loc[len(summary)] = {
        "type": "ALL",
        "n_guides": len(df),
        "n_unique_intended_targets": df["intended_target_name"].nunique(dropna=True),
        "n_unique_gene_names": df["gene_name"].replace("", pd.NA).nunique(dropna=True),
    }
    summary.to_csv(OUT_DIR / "library_counts.tsv", sep="\t", index=False)


def plot_composition(df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    counts = df["type"].value_counts().reindex(TYPE_ORDER)
    sns.barplot(
        x=counts.index,
        y=counts.values,
        hue=counts.index,
        palette=TYPE_PALETTE,
        legend=False,
        ax=axes[0, 0],
    )
    for i, v in enumerate(counts.values):
        axes[0, 0].text(i, v, f"{v:,}", ha="center", va="bottom", fontsize=9)
    axes[0, 0].set_title("Guides per type")
    axes[0, 0].set_ylabel("n guides")
    axes[0, 0].set_xlabel("")
    axes[0, 0].tick_params(axis="x", rotation=20)

    targeting = df[df["type"] == "targeting"]
    guides_per_target = targeting.groupby("intended_target_name").size()
    bins = list(range(1, int(guides_per_target.max()) + 2))
    axes[0, 1].hist(guides_per_target.values, bins=bins, color=TYPE_PALETTE["targeting"], edgecolor="white")
    axes[0, 1].set_title(
        f"Guides per target gene (targeting only)\n"
        f"{guides_per_target.size:,} unique targets, median = {int(guides_per_target.median())}"
    )
    axes[0, 1].set_xlabel("guides per target")
    axes[0, 1].set_ylabel("n targets")

    strand_by_type = df.groupby(["type", "strand"], observed=True).size().unstack(fill_value=0)
    strand_by_type = strand_by_type.reindex(TYPE_ORDER)
    strand_colors = [STRAND_PALETTE.get(s, "#cccccc") for s in strand_by_type.columns]
    strand_by_type.plot(kind="bar", stacked=True, ax=axes[1, 0], color=strand_colors)
    axes[1, 0].set_title("Strand by guide type")
    axes[1, 0].set_ylabel("n guides")
    axes[1, 0].set_xlabel("")
    axes[1, 0].tick_params(axis="x", rotation=20)
    axes[1, 0].legend(title="strand")

    chrom_order = [f"chr{i}" for i in list(range(1, 23)) + ["X", "Y"]]
    chrom_counts = df[df["guide_chr"].isin(chrom_order)]["guide_chr"].value_counts().reindex(chrom_order, fill_value=0)
    axes[1, 1].bar(range(len(chrom_counts)), chrom_counts.values, color="#555555")
    axes[1, 1].set_xticks(range(len(chrom_counts)))
    axes[1, 1].set_xticklabels(chrom_counts.index, rotation=90, fontsize=8)
    axes[1, 1].set_title("Guides per chromosome (all types)")
    axes[1, 1].set_ylabel("n guides")

    fig.suptitle(f"IGVFFI8270UPKB guide library — {len(df):,} guides", fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "library_composition.png", dpi=150)
    plt.close(fig)


def main() -> None:
    df = load_guides()
    write_summary_tsv(df)
    plot_composition(df)
    print(f"wrote {OUT_DIR / 'library_counts.tsv'}")
    print(f"wrote {OUT_DIR / 'library_composition.png'}")


if __name__ == "__main__":
    main()
