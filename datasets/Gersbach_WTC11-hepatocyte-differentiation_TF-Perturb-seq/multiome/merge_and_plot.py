"""
Loads expression files for multiple runs, subsets each to its mapped cells,
concatenates into a single AnnData, then plots coverage diagnostics and
hepatocyte marker violin plots.

Usage:
    python merge_and_plot.py \
        --runs   IGVFFI2923WJHO.h5ad:IGVFSM2174QALS \
                 IGVFFI9930KPDX.h5:IGVFSM0825LHRL \
                 IGVFFI4195HOBE.h5:IGVFSM1624SZQN \
        --mapping       cell_barcode_mapping_all.tsv \
        --raw-runs      IGVFFI2923WJHO.h5ad:IGVFSM2174QALS \
                        IGVFFI9930KPDX.h5:IGVFSM0825LHRL \
                        IGVFFI4195HOBE.h5:IGVFSM1624SZQN \
        --output        merged.h5ad \
        --coverage-plot coverage_plots.png \
        --marker-plot   hepatocyte_markers.png \
        --min-genes     300

Arguments:
    --runs          One or more <expression_file>:<igvfsm_id> pairs.
                    Accepts .h5ad or Cell Ranger .h5 files.
    --mapping       Combined cell_barcode_mapping TSV from assign_multiseq_barcodes.py.
    --raw-runs      Same format as --runs but pointing to raw-count files.
                    Used for coverage QC metrics. If omitted, uses --runs files
                    (only accurate if those files contain raw counts).
    --output        Output merged h5ad path (default: merged.h5ad).
    --coverage-plot Output path for coverage PNG (default: coverage_plots.png).
    --marker-plot   Output path for marker violin PNG (default: hepatocyte_markers.png).
    --skip-norm     Skip normalization (if input files are already normalized).
    --min-genes     Minimum genes detected per cell for marker plots (default: None).
    --genes         Ensembl gene IDs to plot (space-separated).
                    Defaults to AFP/ALB/SERPINA1/APOE/TTR/FGB.
"""

import argparse
import sys

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DEFAULT_GENES = ["AFP", "ALB", "SERPINA1", "APOE", "TTR", "FGB"]
ENSG_TO_SYMBOL = {
    "ENSG00000081051": "AFP",
    "ENSG00000163631": "ALB",
    "ENSG00000197249": "SERPINA1",
    "ENSG00000130203": "APOE",
    "ENSG00000118271": "TTR",
    "ENSG00000171564": "FGB",
}


# ---------------------------------------------------------------------------
# Loading helpers
# ---------------------------------------------------------------------------

def parse_run_specs(specs):
    runs = []
    for spec in specs:
        parts = spec.rsplit(":", 1)
        if len(parts) != 2:
            print(f"ERROR: entries must be <file>:<igvfsm_id>, got: {spec}", file=sys.stderr)
            sys.exit(1)
        runs.append((parts[0], parts[1]))
    return runs


def load_expression(path: str, igvfsm_id: str) -> sc.AnnData:
    """Load a .h5ad or Cell Ranger .h5 file and standardise barcode format.

    For .h5 files, var_names are set to Ensembl IDs (feature/id field) rather
    than gene symbols (feature/name), so gene names are consistent across runs
    regardless of how the h5ad was originally built.
    """
    if path.endswith(".h5ad"):
        adata = sc.read_h5ad(path)
        # Strip version suffixes if present (e.g. ENSG00000081051.2 -> ENSG00000081051)
        adata.var_names = [v.split(".")[0] for v in adata.var_names]
    elif path.endswith(".h5"):
        import h5py
        # Read with gex_only=False to avoid losing features; we filter below
        adata = sc.read_10x_h5(path, gex_only=True)
        # Override var_names with Ensembl IDs from the id field
        with h5py.File(path) as f:
            feature_ids   = f["matrix/features/id"][:].astype(str)
            feature_types = f["matrix/features/feature_type"][:].astype(str)
        # Keep only Gene Expression features (same filter as gex_only=True)
        gex_mask = feature_types == "Gene Expression"
        adata.var_names = feature_ids[gex_mask]
        # Standardise barcodes: strip -1 suffix, append IGVFSM ID
        adata.obs_names = [
            f"{bc.rsplit('-', 1)[0]}_{igvfsm_id}" for bc in adata.obs_names
        ]
    else:
        raise ValueError(f"Unsupported file format: {path}. Expected .h5ad or .h5")
    adata.var_names_make_unique()
    return adata


# ---------------------------------------------------------------------------
# Coverage plot
# ---------------------------------------------------------------------------

def plot_coverage(obs_df: pd.DataFrame, output: str) -> None:
    umi    = obs_df["total_counts"].values
    ngenes = obs_df["n_genes_by_counts"].values

    groupby = "sample_description"
    if groupby in obs_df.columns and obs_df[groupby].notna().any():
        groups  = obs_df[groupby].astype(str)
        labels  = sorted(groups.unique())
        cmap    = plt.get_cmap("tab20", len(labels))
        color_map = {l: cmap(i) for i, l in enumerate(labels)}
        colors  = groups.map(color_map).values
    else:
        groups, labels, color_map, colors = None, None, None, ["steelblue"] * len(umi)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("Sequencing coverage - all runs merged", fontsize=15, fontweight="bold")

    # Knee plot
    ax = axes[0, 0]
    ranked = np.sort(umi)[::-1]
    ax.plot(np.arange(1, len(ranked) + 1), ranked, color="steelblue", linewidth=1.2)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Barcode rank"); ax.set_ylabel("UMI counts")
    ax.set_title("Barcode rank plot (knee plot)", fontweight="bold")
    ax.axhline(np.median(umi), color="tomato", linestyle="--", linewidth=1,
               label=f"Median: {np.median(umi):,.0f}")
    ax.legend(fontsize=9)

    # UMI histogram
    ax = axes[0, 1]
    ax.hist(umi, bins=80, color="steelblue", edgecolor="white", linewidth=0.3)
    ax.axvline(np.median(umi), color="tomato", linestyle="--", linewidth=1.2,
               label=f"Median: {np.median(umi):,.0f}")
    ax.set_xlabel("UMI counts per cell"); ax.set_ylabel("Number of cells")
    ax.set_title("UMI counts per cell", fontweight="bold"); ax.legend(fontsize=9)

    # Genes histogram
    ax = axes[1, 0]
    ax.hist(ngenes, bins=80, color="mediumseagreen", edgecolor="white", linewidth=0.3)
    ax.axvline(np.median(ngenes), color="tomato", linestyle="--", linewidth=1.2,
               label=f"Median: {np.median(ngenes):,.0f}")
    ax.set_xlabel("Genes detected per cell"); ax.set_ylabel("Number of cells")
    ax.set_title("Genes detected per cell", fontweight="bold"); ax.legend(fontsize=9)

    # Scatter
    ax = axes[1, 1]
    if groups is not None:
        for label in labels:
            mask = (groups == label).values
            ax.scatter(umi[mask], ngenes[mask], s=2, alpha=0.4,
                       color=color_map[label], label=label, rasterized=True)
        ax.legend(title="Sample", fontsize=7, title_fontsize=8,
                  markerscale=4, bbox_to_anchor=(1.01, 1), loc="upper left")
    else:
        ax.scatter(umi, ngenes, s=2, alpha=0.3, color="steelblue", rasterized=True)
    ax.set_xlabel("UMI counts per cell"); ax.set_ylabel("Genes detected per cell")
    ax.set_title("UMI counts vs genes detected", fontweight="bold")

    plt.tight_layout()
    plt.savefig(output, dpi=150, bbox_inches="tight")
    print(f"  Saved coverage plot: {output}")


# ---------------------------------------------------------------------------
# Marker violin plot
# ---------------------------------------------------------------------------

def plot_markers(adata: sc.AnnData, genes: list, skip_norm: bool,
                 min_genes: int, output: str) -> None:
    adata_plot = adata
    if min_genes is not None:
        if "n_genes_by_counts" not in adata_plot.obs.columns:
            sc.pp.calculate_qc_metrics(adata_plot, inplace=True)
        before = adata_plot.n_obs
        adata_plot = adata_plot[adata_plot.obs["n_genes_by_counts"] >= min_genes].copy()
        print(f"  --min-genes {min_genes}: {before:,} -> {adata_plot.n_obs:,} cells "
              f"({before - adata_plot.n_obs:,} removed, "
              f"{adata_plot.n_obs/before*100:.1f}% retained for plotting)")

    unversioned_to_full = {v.split(".")[0]: v for v in adata_plot.var_names}
    resolved = {g: unversioned_to_full.get(g, g) for g in genes}
    present  = [resolved[g] for g in genes if resolved[g] in adata_plot.var_names]
    missing  = [g for g in genes if resolved[g] not in adata_plot.var_names]

    if missing:
        print(f"  WARNING: genes not found, skipping: {missing}", file=sys.stderr)
    if not present:
        print("  ERROR: none of the marker genes found.", file=sys.stderr)
        sys.exit(1)

    groupby = None
    if "sample_description" in adata_plot.obs.columns and \
            adata_plot.obs["sample_description"].notna().any():
        groupby = "sample_description"
        print(f"  Grouping by: {groupby}")

    y_label = "log1p(CPM)" if not skip_norm else "log-normalized expression"
    n_genes = len(present)
    fig, axes = plt.subplots(n_genes, 1, figsize=(12, 3 * n_genes))
    if n_genes == 1:
        axes = [axes]

    for ax, gene in zip(axes, present):
        sc.pl.violin(adata_plot, keys=gene, groupby=groupby,
                     ax=ax, show=False, stripplot=False, inner="box")
        legend = ax.get_legend()
        if legend is not None:
            legend.remove()
        unversioned = gene.split(".")[0]
        symbol = ENSG_TO_SYMBOL.get(unversioned)
        ax.set_title(f"{symbol} ({unversioned})" if symbol else gene,
                     fontsize=13, fontweight="bold")
        ax.set_xlabel("")
        ax.set_ylabel(y_label, fontsize=10)
        if groupby:
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)

    fig.suptitle("Hepatocyte marker gene expression", fontsize=14, y=1.01)
    plt.tight_layout()
    plt.savefig(output, dpi=150, bbox_inches="tight")
    print(f"  Saved marker plot: {output}")



# ---------------------------------------------------------------------------
# Marker expression fraction bar plot
# ---------------------------------------------------------------------------

def plot_marker_fractions(adata: sc.AnnData, genes: list, min_genes: int,
                          output: str) -> None:
    """
    For each marker gene, plot the fraction of cells with non-zero expression
    as a grouped bar plot, broken down by timepoint and coloured by run.
    """
    adata_plot = adata
    if min_genes is not None:
        if "n_genes_by_counts" not in adata_plot.obs.columns:
            sc.pp.calculate_qc_metrics(adata_plot, inplace=True)
        adata_plot = adata_plot[
            adata_plot.obs["n_genes_by_counts"] >= min_genes
        ].copy()
        print(f"  Fraction plot: {adata_plot.n_obs:,} cells after --min-genes {min_genes} filter")

    # Resolve versioned var_names
    unversioned_to_full = {v.split(".")[0]: v for v in adata_plot.var_names}
    resolved = {g: unversioned_to_full.get(g, g) for g in genes}
    present  = [resolved[g] for g in genes if resolved[g] in adata_plot.var_names]
    if not present:
        print("  WARNING: no marker genes found for fraction plot.", file=sys.stderr)
        return

    # Extract timepoint from sample_description (e.g. wtc11_hepato_t3_1 -> t3)
    obs = adata_plot.obs.copy()
    obs["timepoint"] = obs["sample_description"].str.extract(r"(t\d+)", expand=False)
    obs["run"]       = obs["igvfsm_id"] if "igvfsm_id" in obs.columns else obs["run"]

    timepoints = sorted(obs["timepoint"].dropna().unique(),
                        key=lambda x: int(x[1:]))  # sort t3 < t8 < t13 < t18
    runs       = sorted(obs["run"].dropna().unique())
    run_colors = {r: plt.get_cmap("tab10")(i) for i, r in enumerate(runs)}

    n_genes = len(present)
    fig, axes = plt.subplots(1, n_genes, figsize=(4 * n_genes, 5), sharey=True)
    if n_genes == 1:
        axes = [axes]

    bar_width = 0.8 / len(runs)

    for ax, gene in zip(axes, present):
        # Get expression vector (already normalized/log; non-zero = expressed)
        expr = pd.Series(
            np.asarray(adata_plot[:, gene].X.todense()).ravel()
            if hasattr(adata_plot[:, gene].X, "todense")
            else adata_plot[:, gene].X.ravel(),
            index=adata_plot.obs_names
        )
        obs["expressed"] = (expr > 0).values

        for r_idx, run in enumerate(runs):
            fracs = []
            for tp in timepoints:
                mask = (obs["timepoint"] == tp) & (obs["run"] == run)
                n_total = mask.sum()
                frac    = obs.loc[mask, "expressed"].mean() if n_total > 0 else np.nan
                fracs.append(frac)

            x = np.arange(len(timepoints)) + r_idx * bar_width
            ax.bar(x, fracs, width=bar_width, color=run_colors[run],
                   label=run, edgecolor="white", linewidth=0.5)

        ax.set_xticks(np.arange(len(timepoints)) + bar_width * (len(runs) - 1) / 2)
        ax.set_xticklabels(timepoints, fontsize=10)
        ax.set_xlabel("Timepoint", fontsize=10)
        ax.set_ylim(0, 1)
        ax.set_yticks(np.arange(0, 1.1, 0.2))

        unversioned = gene.split(".")[0]
        symbol = ENSG_TO_SYMBOL.get(unversioned)
        ax.set_title(f"{symbol}\n({unversioned})" if symbol else gene,
                     fontsize=11, fontweight="bold")

    axes[0].set_ylabel("Fraction of cells with non-zero expression", fontsize=10)

    # Single shared legend
    handles = [plt.Rectangle((0, 0), 1, 1, color=run_colors[r]) for r in runs]
    fig.legend(handles, runs, title="Run (IGVFSM ID)", fontsize=8,
               title_fontsize=9, loc="lower center",
               ncol=len(runs), bbox_to_anchor=(0.5, -0.08))

    fig.suptitle("Fraction of cells expressing hepatocyte markers",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig(output, dpi=150, bbox_inches="tight")
    print(f"  Saved fraction plot: {output}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--runs",     required=True, nargs="+",
                        metavar="FILE:IGVFSM_ID",
                        help="Expression files: <file.h5ad|h5>:<igvfsm_id>")
    parser.add_argument("--mapping",  required=True,
                        help="Combined cell_barcode_mapping_all.tsv")
    parser.add_argument("--raw-runs", nargs="+", default=None,
                        metavar="FILE:IGVFSM_ID",
                        help="Raw-count files for coverage QC (same format as --runs). "
                             "If omitted, --runs files are used for QC metrics.")
    parser.add_argument("--output",        default="merged.h5ad")
    parser.add_argument("--coverage-plot", default="coverage_plots.png")
    parser.add_argument("--marker-plot",   default="hepatocyte_markers.png")
    parser.add_argument("--skip-norm",     action="store_true")
    parser.add_argument("--fraction-plot", default="marker_fractions.png",
                        help="Output path for marker expression fraction bar plot")
    parser.add_argument("--min-genes",     type=int, default=None)
    parser.add_argument("--genes",         nargs="+", default=None)
    return parser.parse_args()


def main():
    args = parse_args()
    runs     = parse_run_specs(args.runs)
    raw_runs = parse_run_specs(args.raw_runs) if args.raw_runs else runs
    genes    = args.genes if args.genes else DEFAULT_GENES

    # Load combined mapping
    print(f"[1/5] Loading mapping: {args.mapping}")
    mapping = pd.read_csv(args.mapping, sep="\t")
    mapping = mapping.dropna(subset=["h5ad_barcode"])
    mapping = mapping[~mapping["multiseq_barcode"].isin(["Doublet", "Unassigned"])]
    print(f"      {len(mapping):,} valid mapped cells across all runs")

    # Load, subset, and concatenate expression files
    print(f"\n[2/5] Loading and subsetting {len(runs)} expression file(s)")
    adatas = []
    for path, igvfsm_id in runs:
        print(f"  {path} ({igvfsm_id})")
        adata = load_expression(path, igvfsm_id)
        print(f"    Loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

        run_mapping = mapping[mapping["igvfsm_id"] == igvfsm_id]
        valid       = set(run_mapping["h5ad_barcode"])
        adata       = adata[adata.obs_names.isin(valid)].copy()
        print(f"    Subset: {adata.n_obs:,} cells")

        # Attach metadata
        meta = run_mapping.set_index("h5ad_barcode")[
            ["multiseq_barcode", "sample_accession", "sample_description", "igvfsm_id"]
        ]
        adata.obs = adata.obs.drop(
            columns=[c for c in meta.columns if c in adata.obs.columns]
        )
        adata.obs = adata.obs.join(meta, how="left")
        adata.obs["run"] = igvfsm_id
        adatas.append(adata)

    print(f"\n[3/5] Concatenating runs")
    for i, adata in enumerate(adatas):
        print(f"  Run {i+1} var_names sample: {list(adata.var_names[:3])}")

    merged = sc.concat(adatas, label="run", keys=[r[1] for r in runs],
                       join="inner", merge="same")
    merged.var_names_make_unique()
    print(f"      Merged: {merged.n_obs:,} cells x {merged.n_vars:,} genes")

    # Normalize (on merged object so scale is consistent)
    if args.skip_norm:
        print(f"\n[4/5] Skipping normalization (--skip-norm set)")
    else:
        print(f"\n[4/5] Normalizing (library-size + log1p)")
        sc.pp.normalize_total(merged, target_sum=1e4)
        sc.pp.log1p(merged)

    merged.write_h5ad(args.output)
    print(f"      Saved merged h5ad: {args.output}")

    # Coverage QC - use raw counts if separate raw files provided
    print(f"\n[5/5] Plotting coverage and markers")
    if args.raw_runs:
        print(f"  Loading raw files for coverage QC metrics")
        raw_adatas = []
        for path, igvfsm_id in raw_runs:
            raw = load_expression(path, igvfsm_id)
            run_mapping = mapping[mapping["igvfsm_id"] == igvfsm_id]
            valid = set(run_mapping["h5ad_barcode"])
            raw   = raw[raw.obs_names.isin(valid)].copy()
            meta  = run_mapping.set_index("h5ad_barcode")[
                ["sample_description", "igvfsm_id"]
            ]
            raw.obs = raw.obs.drop(
                columns=[c for c in meta.columns if c in raw.obs.columns]
            )
            raw.obs = raw.obs.join(meta, how="left")
            raw_adatas.append(raw)
        raw_merged = sc.concat(raw_adatas, label="run", keys=[r[1] for r in raw_runs],
                               join="inner", merge="same")
        sc.pp.calculate_qc_metrics(raw_merged, inplace=True)
        qc_obs = raw_merged.obs
    else:
        print("  NOTE: No --raw-runs provided; coverage metrics reflect normalized values.")
        sc.pp.calculate_qc_metrics(merged, inplace=True)
        qc_obs = merged.obs

    # Per-replicate summary
    if "sample_description" in qc_obs.columns:
        print(f"\n  {'Sample':<25} {'N cells':>8} {'Total UMIs':>12} "
              f"{'Median UMI/cell':>16} {'Median genes/cell':>18}")
        print("  " + "-" * 82)
        for sample, grp in qc_obs.groupby("sample_description", sort=True):
            print(f"  {sample:<25} {len(grp):>8,} "
                  f"{grp['total_counts'].sum():>12,.0f} "
                  f"{grp['total_counts'].median():>16,.0f} "
                  f"{grp['n_genes_by_counts'].median():>18,.0f}")
        tsv_path = args.coverage_plot.replace(".png", "_per_replicate.tsv")
        (qc_obs.groupby("sample_description")
               .agg(n_cells=("total_counts","count"),
                    total_umis=("total_counts","sum"),
                    mean_umis_per_cell=("total_counts","mean"),
                    median_umis_per_cell=("total_counts","median"),
                    median_genes_per_cell=("n_genes_by_counts","median"))
               .reset_index()
               .to_csv(tsv_path, sep="\t", index=False, float_format="%.1f"))
        print(f"\n  Per-replicate table: {tsv_path}")

    plot_coverage(qc_obs, args.coverage_plot)

    # Marker violins use the normalized merged object
    if "n_genes_by_counts" not in merged.obs.columns:
        sc.pp.calculate_qc_metrics(merged, inplace=True)
    plot_markers(merged, genes, args.skip_norm, args.min_genes, args.marker_plot)

    plot_marker_fractions(merged, genes, args.min_genes, args.fraction_plot)

    print("\nDone.")


if __name__ == "__main__":
    main()