#!/usr/bin/env python3
"""
shared_regulators.py

For each gene set category from investigate_effects_on_hepatocyte_markers.py,
find perturbations (TFs) that significantly regulate multiple marker genes
within that category in a consistent direction (all up or all down).

For each gene set:
  - Load the per-TF result CSV for each marker gene
  - Apply FDR and log2FC thresholds to call a TF significant for each gene
  - Classify direction: "down" (log2fc < -fc_thresh), "up" (log2fc > fc_thresh)
  - For each TF, count how many markers it regulates and in which direction
  - Output a ranked table: one row per TF, columns = n_markers_regulated,
    direction_consistent (all same direction), list of regulated markers,
    median log2FC, min FDR

Output files (one per gene set):
  shared_regulators_{gene_set_stem}.csv

Usage:
    python shared_regulators.py \\
        --input_dir /path/to/output_dir \\
        --output_dir /path/to/output_dir \\
        --suffix perturbo_trans      # or wilcoxon
        --fdr_thresh 0.05 \\
        --fc_thresh 0.5 \\
        --min_markers 2              # min markers regulated to include a TF
"""

import argparse
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

# --- Gene sets (must match investigate_effects_on_hepatocyte_markers.py) ------
GENE_SETS = {
    "hepatocyte_markers": {
        "label": "Hepatocyte maturity markers",
        "markers": ["ALB", "SERPINA1", "APOE", "AFP", "TTR", "FGB", "ARG1"],
    },
    "cyp450": {
        "label": "CYP450 drug-metabolizing enzymes",
        "markers": ["CYP3A4", "CYP3A5", "CYP2D6", "CYP2C9", "CYP2C19", "CYP1A2", "CYP2E1"],
    },
    "ldl_cholesterol": {
        "label": "LDL and cholesterol metabolism",
        "markers": ["LDLR", "APOB", "PCSK9", "LDLRAP1", "SORT1", "CYP7A1", "ABCA1"],
    },
    "zonation": {
        "label": "Liver zonation markers",
        "markers": ["GLUL", "CYP2E1", "OAT", "SLC1A2", "CYP2A6",
                    "SDS", "CYP2F1", "HAL", "HSD17B13", "ALDH1B1",
                    "GSN", "COL1A2", "VIM"],
    },
    "cholangiocyte": {
        "label": "Cholangiocyte markers",
        "markers": ["KRT19", "KRT7"],
    },
    "mash": {
        "label": "MASH-associated genes (PMC11736312)",
        "markers": ["PNPLA3", "TM6SF2", "GCKR", "MBOAT7", "HSD17B13"],
    },
    "nash": {
        "label": "NASH-associated genes (Frontiers Genetics 2023)",
        "markers": ["FABP5", "SCD", "CCL20", "GPAT3", "PLIN1", "IL1RN"],
    },
}


# --- CLI ----------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--input_dir",   required=True,
                   help="Directory containing per-marker CSV files from "
                        "investigate_effects_on_hepatocyte_markers.py")
    p.add_argument("--output_dir",  default=None,
                   help="Output directory (default: same as input_dir)")
    p.add_argument("--suffix",      default="perturbo_trans",
                   choices=["perturbo_trans", "wilcoxon"],
                   help="File suffix used when generating CSVs (default: perturbo_trans)")
    p.add_argument("--fdr_thresh",  type=float, default=0.05)
    p.add_argument("--fc_thresh",   type=float, default=0.5)
    p.add_argument("--min_markers", type=int,   default=2,
                   help="Minimum number of markers a TF must regulate to be reported "
                        "(default: 2)")
    return p.parse_args()


# --- Helpers ------------------------------------------------------------------
def csv_path(input_dir, suffix, marker):
    """Construct the expected CSV path for a given marker."""
    if suffix == "perturbo_trans":
        fname = f"perturbo_trans_{marker.lower()}_per_element.csv"
    else:
        fname = f"wilcoxon_{marker.lower()}_per_tf.csv"
    return os.path.join(input_dir, fname)


def load_marker_results(input_dir, suffix, markers, fdr_thresh, fc_thresh):
    """
    Load per-TF results for each marker in the list.
    Returns a dict: {marker_name -> DataFrame with target_symbol, log2fc, padj, direction}
    Only includes rows passing fdr_thresh and abs(log2fc) > fc_thresh.
    """
    loaded = {}
    for marker in markers:
        path = csv_path(input_dir, suffix, marker)
        if not os.path.exists(path):
            print(f"  WARNING: {path} not found -- skipping {marker}")
            continue
        df = pd.read_csv(path)
        if "target_symbol" not in df.columns or "padj" not in df.columns:
            print(f"  WARNING: unexpected columns in {path}: {list(df.columns)}")
            continue
        sig = df[
            (df["padj"] < fdr_thresh) &
            (df["log2fc"].abs() > fc_thresh)
        ].copy()
        sig["direction"] = np.where(sig["log2fc"] > 0, "up", "down")
        loaded[marker] = sig
        print(f"  {marker}: {len(sig):,} significant TFs "
              f"(up={( sig['direction']=='up').sum()}, "
              f"down={( sig['direction']=='down').sum()})")
    return loaded


def find_shared_regulators(marker_results, min_markers):
    """
    For each TF that appears in at least min_markers marker result sets,
    compute:
      - n_markers_total: number of markers it significantly regulates
      - n_markers_down:  number regulated downward
      - n_markers_up:    number regulated upward
      - direction_consistent: True if all regulated markers go same direction
      - dominant_direction:   "down", "up", or "mixed"
      - markers_down:    comma-separated list of downregulated markers
      - markers_up:      comma-separated list of upregulated markers
      - median_log2fc:   median log2FC across regulated markers
      - min_padj:        smallest adjusted p-value across regulated markers
    """
    # Collect all TF symbols across all markers
    all_tfs = set()
    for df in marker_results.values():
        all_tfs.update(df["target_symbol"].unique())

    rows = []
    for tf in all_tfs:
        down_markers = []
        up_markers   = []
        lfcs         = []
        padjs        = []

        for marker, df in marker_results.items():
            row = df[df["target_symbol"] == tf]
            if len(row) == 0:
                continue
            # Take the most significant row if duplicates exist
            row = row.nsmallest(1, "padj").iloc[0]
            if row["direction"] == "down":
                down_markers.append(marker)
            else:
                up_markers.append(marker)
            lfcs.append(row["log2fc"])
            padjs.append(row["padj"])

        n_total = len(down_markers) + len(up_markers)
        if n_total < min_markers:
            continue

        consistent = len(down_markers) == 0 or len(up_markers) == 0
        if len(down_markers) >= len(up_markers):
            dominant = "down"
        else:
            dominant = "up"
        if not consistent:
            dominant = "mixed"

        rows.append({
            "target_symbol":         tf,
            "n_markers_regulated":   n_total,
            "n_markers_down":        len(down_markers),
            "n_markers_up":          len(up_markers),
            "direction_consistent":  consistent,
            "dominant_direction":    dominant,
            "markers_down":          ", ".join(sorted(down_markers)),
            "markers_up":            ", ".join(sorted(up_markers)),
            "median_log2fc":         float(np.median(lfcs)),
            "min_padj":              float(np.min(padjs)),
        })

    if not rows:
        return pd.DataFrame()

    df_out = (
        pd.DataFrame(rows)
          .sort_values(["n_markers_regulated", "min_padj"],
                       ascending=[False, True])
          .reset_index(drop=True)
    )
    return df_out


def plot_shared_regulators(df, gene_set_label, out_path, top_n=30):
    """
    Horizontal bar chart: top TFs by number of markers regulated.
    Bars colored by dominant direction (blue=down, red=up, grey=mixed).
    Separate bars for down and up counts, stacked.
    """
    if df.empty:
        return

    top = df.head(top_n).sort_values("n_markers_regulated", ascending=True)

    dir_colors = {"down": "#1f77b4", "up": "#d62728", "mixed": "#888888"}
    bar_colors = [dir_colors[d] for d in top["dominant_direction"]]

    fig, ax = plt.subplots(figsize=(9, max(4, len(top) * 0.4)))

    # Stacked bar: down (blue) + up (red)
    ax.barh(range(len(top)), top["n_markers_down"],
            color="#1f77b4", label="down", height=0.7)
    ax.barh(range(len(top)), top["n_markers_up"],
            left=top["n_markers_down"],
            color="#d62728", label="up", height=0.7)

    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top["target_symbol"], fontsize=9)
    ax.set_xlabel("Number of marker genes regulated", fontsize=11)
    ax.set_title(
        f"Shared regulators: {gene_set_label}\n"
        f"(blue = markers downregulated, red = upregulated)",
        fontsize=11, fontweight="bold",
    )
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(fontsize=9, frameon=False)

    # Annotate with total count
    for i, (_, row) in enumerate(top.iterrows()):
        ax.text(row["n_markers_regulated"] + 0.05, i,
                str(int(row["n_markers_regulated"])),
                va="center", fontsize=8)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Plot saved: {out_path}")


# --- Main ---------------------------------------------------------------------
def main():
    args    = parse_args()
    out_dir = args.output_dir or args.input_dir
    os.makedirs(out_dir, exist_ok=True)

    print(f"Input directory  : {args.input_dir}")
    print(f"Output directory : {out_dir}")
    print(f"Suffix           : {args.suffix}")
    print(f"FDR threshold    : {args.fdr_thresh}")
    print(f"log2FC threshold : {args.fc_thresh}")
    print(f"Min markers      : {args.min_markers}")

    for stem, info in GENE_SETS.items():
        label   = info["label"]
        markers = info["markers"]
        print(f"\n{'='*60}")
        print(f"{label}")
        print(f"{'='*60}")

        marker_results = load_marker_results(
            args.input_dir, args.suffix, markers,
            args.fdr_thresh, args.fc_thresh,
        )

        if len(marker_results) < 2:
            print(f"  Fewer than 2 markers loaded -- skipping shared regulator analysis")
            continue

        shared = find_shared_regulators(marker_results, args.min_markers)

        if shared.empty:
            print(f"  No TFs regulate >= {args.min_markers} markers at these thresholds")
            continue

        print(f"\n  {len(shared):,} TFs regulate >= {args.min_markers} markers")
        print(f"  Top 10:")
        print(shared.head(10)[
            ["target_symbol", "n_markers_regulated", "dominant_direction",
             "markers_down", "markers_up", "median_log2fc", "min_padj"]
        ].to_string(index=False))

        # Save CSV
        csv_out = os.path.join(out_dir, f"shared_regulators_{stem}.csv")
        shared.to_csv(csv_out, index=False)
        print(f"\n  Saved: {csv_out}")

        # Save plot
        png_out = os.path.join(out_dir, f"shared_regulators_{stem}.png")
        plot_shared_regulators(shared, label, png_out)

    print("\n=== Complete ===")


if __name__ == "__main__":
    main()