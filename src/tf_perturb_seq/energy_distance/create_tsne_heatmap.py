import argparse
import matplotlib
matplotlib.use("Agg")
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
import umap
import scanpy as sc
import anndata
from pathlib import Path

parser = argparse.ArgumentParser(description="Remake pt001 t-SNE/UMAP and heatmap with gene symbols.")
parser.add_argument("energy_distance_dir", nargs="?", default=".",
                    help="Path to energy_distance/ folder (default: current directory)")
parser.add_argument("--outdir", default=None,
                    help="Directory to write output files (default: same as energy_distance_dir)")
parser.add_argument("--leiden_resolution", type=float, default=0.5,
                    help="Leiden clustering resolution for UMAP plot (default: 0.5)")
args = parser.parse_args()

eddir = Path(args.energy_distance_dir)
outdir = Path(args.outdir) if args.outdir else eddir
outdir.mkdir(parents=True, exist_ok=True)

# --- Build TSS -> gene symbol mapping from all_guides.csv ---
guides = pd.read_csv(eddir / "all_guides.csv")
guides["gene_symbol"] = guides["guide_id"].str.split("#").str[0]
tss_to_symbol = guides.groupby("target_tss")["gene_symbol"].first().to_dict()

def tss_label(tss):
    sym = tss_to_symbol.get(tss)
    coord = tss.split("|")[1] if "|" in tss else tss
    return f"{sym}|{coord}" if sym else tss

# -- Figure 1: t-SNE using Anish's precomputed coordinates -------------------
print("Loading precomputed t-SNE...")
tsne_df = pd.read_csv(eddir / "pt001_tss_estat_emb_tsne_data.csv")
tsne_df["label"] = tsne_df["tss"].apply(tss_label)
n = len(tsne_df)

# Map cluster IDs to colors
clusters = tsne_df["cluster"].values
unique_clusters = sorted(set(clusters))
cmap = plt.get_cmap("tab20", len(unique_clusters))
cluster_colors = {c: cmap(i) for i, c in enumerate(unique_clusters)}
colors = [cluster_colors[c] for c in clusters]

fig, ax = plt.subplots(figsize=(14, 12))
ax.scatter(tsne_df["x"], tsne_df["y"], c=colors, s=18, alpha=0.8, linewidths=0)
for _, row in tsne_df.iterrows():
    ax.annotate(row["label"], (row["x"], row["y"]),
                fontsize=4.5, alpha=0.75,
                xytext=(3, 3), textcoords="offset points")
ax.set_xlabel("t-SNE 1")
ax.set_ylabel("t-SNE 2")
ax.set_title(f"Energy distance embedding ({n} TSS targets, p<0.001)", fontsize=11)
ax.spines[["top", "right"]].set_visible(False)
plt.tight_layout()
plt.savefig(outdir / "pt001_tsne_symbols.pdf", bbox_inches="tight")
plt.savefig(outdir / "pt001_tsne_symbols.png", dpi=150, bbox_inches="tight")
plt.close()
print("Saved t-SNE.")

# -- Figure 2: UMAP with Leiden clustering -----------------------------------
print("Running UMAP + Leiden clustering...")
pw = pd.read_csv(eddir / "pt001_tss_estat_pairwise.csv", index_col=0)
pw_sym = (pw.values + pw.values.T) / 2
np.fill_diagonal(pw_sym, 0)
tss_ids = pw.index.tolist()
labels = [tss_label(t) for t in tss_ids]

# Build AnnData from distance matrix for scanpy Leiden
adata = anndata.AnnData(X=pw_sym)
adata.obs_names = labels
sc.pp.neighbors(adata, use_rep="X", metric="precomputed", n_neighbors=15)
sc.tl.leiden(adata, resolution=args.leiden_resolution)
sc.tl.umap(adata)

leiden_labels = adata.obs["leiden"].values
unique_leiden = sorted(set(leiden_labels), key=int)
cmap2 = plt.get_cmap("tab20", len(unique_leiden))
lcluster_colors = {c: cmap2(i) for i, c in enumerate(unique_leiden)}
lcolors = [lcluster_colors[c] for c in leiden_labels]

umap_coords = adata.obsm["X_umap"]
fig, ax = plt.subplots(figsize=(14, 12))
ax.scatter(umap_coords[:, 0], umap_coords[:, 1],
           c=lcolors, s=18, alpha=0.8, linewidths=0)
for i, lab in enumerate(labels):
    ax.annotate(lab, (umap_coords[i, 0], umap_coords[i, 1]),
                fontsize=4.5, alpha=0.75,
                xytext=(3, 3), textcoords="offset points")
ax.set_xlabel("UMAP 1")
ax.set_ylabel("UMAP 2")
ax.set_title(f"Energy distance UMAP ({n} TSS, p<0.001, Leiden r={args.leiden_resolution})", fontsize=11)
ax.spines[["top", "right"]].set_visible(False)
plt.tight_layout()
plt.savefig(outdir / "pt001_umap_leiden_symbols.pdf", bbox_inches="tight")
plt.savefig(outdir / "pt001_umap_leiden_symbols.png", dpi=150, bbox_inches="tight")
plt.close()
print("Saved UMAP.")

# -- Figure 3: Clustered heatmap (fixed colorscale) --------------------------
print("Clustering heatmap...")
condensed = squareform(pw_sym, checks=False)
Z = linkage(condensed, method="average")
order = dendrogram(Z, no_plot=True)["leaves"]
pw_ordered = pw_sym[np.ix_(order, order)]
ordered_labels = [labels[i] for i in order]

fig, ax = plt.subplots(figsize=(max(14, n * 0.12), max(12, n * 0.12)))
im = ax.imshow(pw_ordered, aspect="auto",
               norm=mcolors.PowerNorm(gamma=0.4, vmin=0, vmax=pw_sym.max()),
               cmap="magma")
plt.colorbar(im, ax=ax, fraction=0.02, pad=0.02, label="Energy distance")
tick_fs = max(3, min(7, 200 // n))
ax.set_xticks(range(n))
ax.set_xticklabels(ordered_labels, rotation=90, fontsize=tick_fs)
ax.set_yticks(range(n))
ax.set_yticklabels(ordered_labels, fontsize=tick_fs)
ax.set_title(f"Pairwise energy distance heatmap ({n} TSS, p<0.001)", fontsize=11)
plt.tight_layout()
plt.savefig(outdir / "pt001_heatmap_symbols.pdf", bbox_inches="tight")
plt.savefig(outdir / "pt001_heatmap_symbols.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"All outputs written to: {outdir}")