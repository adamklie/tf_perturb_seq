# -*- coding: utf-8 -*-
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
from statsmodels.stats.multitest import multipletests
import scanpy as sc
import anndata
from adjustText import adjust_text
from matplotlib.colors import to_hex


parser = argparse.ArgumentParser(description="Energy distance QC and UMAP plots.")
parser.add_argument("energy_distance_dir", nargs="?", default=".",
                    help="Path to energy_distance/ folder (default: current directory)")
parser.add_argument("--outdir", default=None,
                    help="Directory to write output files (default: same as energy_distance_dir)")
parser.add_argument("--fdr_thresh", type=float, default=1e-5,
                    help="FDR threshold for significant TSS in UMAP (default: 1e-5)")
parser.add_argument("--leiden_resolution", type=float, default=0.5,
                    help="Leiden clustering resolution (default: 0.5)")
args = parser.parse_args()

eddir  = Path(args.energy_distance_dir)
outdir = Path(args.outdir) if args.outdir else eddir
outdir.mkdir(parents=True, exist_ok=True)

# --- Gene symbol mapping ---
guides = pd.read_csv(eddir / "all_guides.csv")
guides["symbol"] = guides["guide_id"].str.split("#").str[0]
tss_to_symbol = guides.groupby("target_tss")["symbol"].first().to_dict()

def tss_to_gene(tss):
    return tss_to_symbol.get(tss, tss.split("|")[0])

# -- Figure 1 & 2: Outlier guide summary -------------------------------------
df = pd.read_csv(eddir / "all_guides.csv")
counts = (
    df.groupby(["target_tss", "is_outlier"])
    .size()
    .unstack(fill_value=0)
    .rename(columns={False: "inlier", True: "outlier"})
)
counts = counts.sort_values("outlier", ascending=False)

outlier_counts = counts["outlier"]
max_outliers   = int(outlier_counts.max())
bins = np.arange(-0.5, max_outliers + 1.5, 1)

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

ax = axes[0]
ax.hist(outlier_counts, bins=bins, color="#4C72B0", edgecolor="white", linewidth=0.6)
ax.set_xlabel("Number of outlier guides per TSS")
ax.set_ylabel("Number of TSS targets")
ax.set_title("Distribution of outlier guide counts per TSS")
ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
ax.spines[["top", "right"]].set_visible(False)
n_with_outlier = (outlier_counts >= 1).sum()
pct = 100 * n_with_outlier / len(outlier_counts)
ax.text(
    0.97, 0.97,
    f"{n_with_outlier} / {len(outlier_counts)} TSS\nhave >=1 outlier ({pct:.1f}%)",
    transform=ax.transAxes, ha="right", va="top",
    fontsize=9, color="#333333",
    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#cccccc", lw=0.8),
)

counts["total"]        = counts["inlier"] + counts["outlier"]
counts["outlier_frac"] = counts["outlier"] / counts["total"]
multi = counts[counts["total"] >= 2].copy()

ax2 = axes[1]
ax2.hist(multi["outlier_frac"], bins=20, color="#DD8452", edgecolor="white", linewidth=0.6)
ax2.set_xlabel("Fraction of guides that are outliers (TSS with >=2 guides)")
ax2.set_ylabel("Number of TSS targets")
ax2.set_title("Outlier fraction per TSS")
ax2.xaxis.set_major_formatter(ticker.PercentFormatter(xmax=1))
ax2.spines[["top", "right"]].set_visible(False)

plt.tight_layout()
plt.savefig(outdir / "outlier_guides_summary.pdf", bbox_inches="tight")
plt.savefig(outdir / "outlier_guides_summary.png", dpi=150, bbox_inches="tight")
plt.close()
print("Saved: outlier_guides_summary.pdf / .png")

# Intermediate TSV for R
outlier_r = counts.reset_index().rename(columns={"target_tss": "target_tss"})
outlier_r["gene_symbol"] = outlier_r["target_tss"].map(tss_to_symbol).fillna(
    outlier_r["target_tss"].str.split("|").str[0]
)
outlier_r.to_csv(outdir / "outlier_guides_summary.tsv", sep="\t", index=False)
print("Saved R-ready: outlier_guides_summary.tsv")

print("\nTop 20 TSS by outlier guide count:")
print(counts[["outlier", "total", "outlier_frac"]].head(20).to_string())

# -- UMAP + Leiden at TF level ------------------------------------------------
print(f"\nRunning UMAP with FDR < {args.fdr_thresh} threshold...")

# Identify significant TSS
summary = pd.read_csv(eddir / "edist_per_tss_summary.csv", index_col=0)
_, padj, _, _ = multipletests(summary["pval_mean"], method="fdr_bh")
summary["padj"] = padj
sig_tss = summary[summary["padj"] < args.fdr_thresh].index.tolist()
n_tfs   = len(set(tss_to_gene(t) for t in sig_tss))
print(f"Significant TSS: {len(sig_tss)}, unique TFs: {n_tfs}")

# Load and subset pairwise matrix
pw = pd.read_csv(eddir / "pt001_tss_estat_pairwise.csv", index_col=0)
missing = [t for t in sig_tss if t not in pw.index]
if missing:
    print(f"Warning: {len(missing)} significant TSS not found in pairwise matrix, dropping.")
    sig_tss = [t for t in sig_tss if t in pw.index]

pw_sub = pw.loc[sig_tss, sig_tss]
pw_sym = (pw_sub.values + pw_sub.values.T) / 2
np.fill_diagonal(pw_sym, 0)

# Collapse to TF level by averaging pairwise distances
tss_genes = [tss_to_gene(t) for t in sig_tss]
pw_df     = pd.DataFrame(pw_sym, index=tss_genes, columns=tss_genes)

unique_genes = sorted(set(tss_genes))
n_genes      = len(unique_genes)
tf_dist = pd.DataFrame(
    np.zeros((n_genes, n_genes)),
    index=unique_genes, columns=unique_genes,
)
for i, g1 in enumerate(unique_genes):
    for j, g2 in enumerate(unique_genes):
        if i == j:
            tf_dist.loc[g1, g2] = 0.0
        elif i < j:
            val = (
                pw_df.loc[g1, g2].values.mean()
                if hasattr(pw_df.loc[g1, g2], "values")
                else pw_df.loc[g1, g2]
            )
            tf_dist.loc[g1, g2] = val
            tf_dist.loc[g2, g1] = val

tf_mat = tf_dist.values.astype(float)
print(f"TF-level distance matrix: {tf_mat.shape[0]} x {tf_mat.shape[1]}")

# UMAP + Leiden
adata = anndata.AnnData(X=tf_mat)
adata.obs_names = unique_genes
sc.pp.neighbors(adata, use_rep="X", metric="precomputed",
                n_neighbors=min(15, n_genes - 1))
sc.tl.leiden(adata, resolution=args.leiden_resolution)
sc.tl.umap(adata)

leiden_labels = adata.obs["leiden"].values
unique_leiden = sorted(set(leiden_labels), key=int)
cmap          = plt.get_cmap("tab20", len(unique_leiden))
lcolors       = {c: cmap(i) for i, c in enumerate(unique_leiden)}
lcolors_hex   = {k: to_hex(v) for k, v in lcolors.items()}
colors        = [lcolors[c] for c in leiden_labels]

umap_coords = adata.obsm["X_umap"]

# Intermediate TSV for R (UMAP coordinates + Leiden + colors)
umap_r = pd.DataFrame({
    "gene_symbol": unique_genes,
    "UMAP1":       umap_coords[:, 0],
    "UMAP2":       umap_coords[:, 1],
    "leiden":      leiden_labels,
    "color_hex":   [lcolors_hex[c] for c in leiden_labels],
})
umap_r.to_csv(outdir / "tf_umap_leiden_ggplot.tsv", sep="\t", index=False)
print("Saved R-ready UMAP: tf_umap_leiden_ggplot.tsv")

# Taller than wide (16 wide x 20 tall)
fig, ax = plt.subplots(figsize=(16, 20))
ax.scatter(umap_coords[:, 0], umap_coords[:, 1],
           c=colors, s=40, alpha=0.8, linewidths=0)

texts = []
for i, gene in enumerate(unique_genes):
    texts.append(ax.text(umap_coords[i, 0], umap_coords[i, 1],
                         gene, fontsize=9, alpha=0.9))
adjust_text(
    texts,
    x=umap_coords[:, 0],
    y=umap_coords[:, 1],
    ax=ax,
    arrowprops=dict(arrowstyle="-", color="gray", lw=0.5),
    expand=(1.2, 1.4),
)

ax.set_xlabel("UMAP 1", fontsize=20)
ax.set_ylabel("UMAP 2", fontsize=20)
ax.tick_params(axis="both", labelsize=16)
ax.set_title(
    f"TF energy distance UMAP ({n_genes} TFs, FDR<{args.fdr_thresh}, "
    f"Leiden r={args.leiden_resolution})",
    fontsize=20,
)
ax.spines[["top", "right"]].set_visible(False)
plt.tight_layout()

fdr_str = f"{args.fdr_thresh:.0e}".replace("-0", "-")
plt.savefig(outdir / f"umap_tf_fdr{fdr_str}_leiden.pdf", bbox_inches="tight")
plt.savefig(outdir / f"umap_tf_fdr{fdr_str}_leiden.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved UMAP: umap_tf_fdr{fdr_str}_leiden.pdf / .png")

# Save annotated AnnData and Leiden classification
adata.write_h5ad(outdir / "tf_umap_leiden.h5ad")
print("Saved anndata: tf_umap_leiden.h5ad")

leiden_classification_sorted = (
    adata.obs["leiden"].astype(int).sort_values().to_frame(name="leiden")
)
leiden_classification_sorted["color"] = (
    leiden_classification_sorted["leiden"].astype(str).map(lcolors_hex)
)
leiden_classification_sorted.to_csv(
    outdir / "tf_umap_leiden_classification_wcolors.csv"
)
print("Saved: tf_umap_leiden_classification_wcolors.csv")
