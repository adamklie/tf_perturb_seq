"""
cNMF Program Comparison Across Datasets
========================================
For each available dataset:
  1. Load gene spectra TPM (program gene weights)
  2. Load cell usages (cell x program)
  3. Load cell -> guide assignment from h5mu
  4. Compute mean program usage per TF perturbation

Then:
  - Match programs across datasets by cosine similarity (Hungarian algorithm)
  - For a chosen TF of interest, plot mean usage across matched programs,
    one column per dataset
  - Plot cross-dataset program similarity heatmap (how well programs match)
  - Plot cross-dataset TF-program correlation (do same TFs regulate same programs?)

Usage:
  python cnmf_comparison.py
  # Edit the CONFIG section below to set paths, k, dt, and TF of interest
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import mudata
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
from scipy.stats import zscore

# -
# CONFIG - edit this section
# -

BASE_DIR = "/hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets"
OUT_DIR  = "/hpc/group/gersbachlab/seg95/tf_perturb_seq/src/tf_perturb_seq/gene_program_discovery/benchmark_figures_cnmf"
os.makedirs(OUT_DIR, exist_ok=True)

# One entry per dataset that has cNMF results.
# run_name:    the prefix before .gene_spectra_tpm.k_*.dt_*.txt (also the
#              PerturbNMF/Result subfolder these files live in)
# k:           chosen number of programs
# dt:          chosen distance threshold (use underscores: "2_0" not "2.0")
# eval_subdir: the "{k}_{eval_subdir}" folder under PerturbNMF/Result/<run_name>/Eval
#              that the Evaluation step writes to (confirmed from the
#              run_torch_cnmf_apptainer_*.sh scripts' --perturb_path_base,
#              e.g. "Eval/30_2_0/30_CRT" -> eval_subdir="2_0")
DATASETS = [
    {
        "name":        "CC-Perturb (Engreitz)",
        "dir":         "Engreitz_WTC11-benchmark_TF-Perturb-seq",
        "run_name":    "070226_torchcNMF_engreitz",
        "k":           30,
        "dt":          "2_0",
        "eval_subdir": "2_0",
    },
    {
        "name":        "GEM-X v3 10xv3 5' (Gersbach)",
        "dir":         "Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3",
        "run_name":    "070226_torchcNMF_gersbach_gemx",
        "k":           30,
        "dt":          "2_0",
        "eval_subdir": "2_0",
    },
    {
        "name":        "HTv2 10xv3 5' (Gersbach)",
        "dir":         "Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2",
        "run_name":    "070226_torchcNMF_gersbach_htv2",
        "k":           30,
        "dt":          "2_0",
        "eval_subdir": "2_0",
    },
    {
        "name":        "HTO 10xv3 5' (Hon)",
        "dir":         "Hon_WTC11-benchmark_TF-Perturb-seq",
        "run_name":    "070226_torchcNMF_hon",
        "k":           30,
        "dt":          "2_0",
        "eval_subdir": "2_0",
    },
    {
        "name":        "10xv3 3' (Huangfu)",
        "dir":         "Huangfu_WTC11-benchmark_TF-Perturb-seq",
        "run_name":    "070226_torchcNMF_huangfu",
        "k":           30,
        "dt":          "2_0",
        "eval_subdir": "2_0",
    },
]

# Color palette matching the energy distance script
PALETTE = {
    "CC-Perturb":        "#7030a0",
    "GEM-x v3 10xv3 5'": "#00b050",
    "HTv2 10xv3 5'":     "#ff4b3d",
    "HTO 10xv3 5'":      "#00b0f0",
    "10xv3 3'":          "#ebce9a",
}

# TF of interest for the per-TF usage plot.
# Must match gene_name values in the guide modality.
# Change this to explore different TFs.
TF_OF_INTEREST = "EZH2"

# Non-targeting label (used to compute baseline usage for normalization)
NON_TARGETING_LABEL = "non-targeting"

# Number of top genes to use for program matching (cosine similarity)
N_TOP_GENES_MATCH = 200

# Minimum number of cells required per TF perturbation to include it
MIN_CELLS_PER_TF = 10

# -
# Helper functions
# -

def load_cnmf_files(ds):
    """Load gene spectra TPM and usages for one dataset."""
    result_dir = os.path.join(BASE_DIR, ds["dir"], "PerturbNMF", "Result", ds["run_name"])
    k, dt, run = ds["k"], ds["dt"], ds["run_name"]

    tpm_path    = os.path.join(result_dir,
                               f"{run}.gene_spectra_tpm.k_{k}.dt_{dt}.txt")
    usage_path  = os.path.join(result_dir,
                               f"{run}.usages.k_{k}.dt_{dt}.consensus.txt")

    if not os.path.exists(tpm_path):
        raise FileNotFoundError(f"TPM file not found:\n  {tpm_path}")
    if not os.path.exists(usage_path):
        raise FileNotFoundError(f"Usage file not found:\n  {usage_path}")

    # gene_spectra_tpm: rows = programs, columns = genes (ENSG IDs)
    spectra = pd.read_csv(tpm_path, sep="\t", index_col=0)
    # usages: rows = cells, columns = programs
    usages  = pd.read_csv(usage_path, sep="\t", index_col=0)

    print(f"  [{ds['name']}] spectra: {spectra.shape}  usages: {usages.shape}")
    return spectra, usages


def load_cell_guide_assignments(ds):
    """
    Load cell -> gene_name assignment from h5mu.
    Returns a Series indexed by cell barcode with gene_name values.
    Cells with no single assigned guide are labeled NA.
    """
    h5mu_path = os.path.join(BASE_DIR, ds["dir"], "sceptre",
                             "inference_mudata_cleaned.h5mu")
    if not os.path.exists(h5mu_path):
        h5mu_path = os.path.join(BASE_DIR, ds["dir"], "sceptre",
                             "inference_mudata.h5mu")

    print(f"  [{ds['name']}] loading h5mu (may take a moment)...")
    mdata = mudata.read_h5mu(h5mu_path, backed="r")

    guide_mod  = mdata["guide"]
    assignment = guide_mod.layers["guide_assignment"]   # cell x guide sparse
    gene_names = guide_mod.var["gene_name"].values       # per guide
    cell_bcs   = guide_mod.obs_names.tolist()

    # For each cell, find guides with assignment == 1
    # guide_assignment is typically a binary sparse matrix
    import scipy.sparse as sp
    if sp.issparse(assignment):
        assignment = assignment.toarray()

    # Each cell gets the gene_name of its assigned guide(s).
    # If multiple guides assigned, take the first non-NT one;
    # if all NT or none, label as non-targeting / NA.
    assigned_tf = []
    for i in range(assignment.shape[0]):
        row      = assignment[i]
        hit_idx  = np.where(row > 0)[0]
        if len(hit_idx) == 0:
            assigned_tf.append(np.nan)
        else:
            names = [gene_names[j] for j in hit_idx]
            # Prefer non-targeting label if that's all there is
            non_nt = [n for n in names if n != NON_TARGETING_LABEL]
            if non_nt:
                assigned_tf.append(non_nt[0])
            else:
                assigned_tf.append(NON_TARGETING_LABEL)

    cell_tf = pd.Series(assigned_tf, index=cell_bcs, name="gene_name")

    mdata.file.close()
    print(f"  [{ds['name']}] {cell_tf.notna().sum()} / {len(cell_tf)} cells assigned")
    return cell_tf


def load_ensembl_symbol_map(ds):
    """Load ENSG -> gene symbol mapping from the Data directory."""
    csv_path = os.path.join(BASE_DIR, ds["dir"], "sceptre",
                            "ensembl_to_symbol.csv")
    if os.path.exists(csv_path):
        m = pd.read_csv(csv_path, index_col=0)
        # Expect columns: ensembl_id, gene_name (or similar)
        # Use the first column as key, second as value
        cols = m.columns.tolist()
        return dict(zip(m.index, m.iloc[:, 0]))
    return {}

def load_perturbation_associations(ds):

    path = os.path.join(
        BASE_DIR,
        ds["dir"],
        "PerturbNMF",
        "Result",
        ds["run_name"],
        "Eval",
        f"{ds['k']}_{ds['eval_subdir']}",
        f"{ds['k']}_perturbation_association_results_WTC.txt"
    )

    df = pd.read_csv(path, sep="\t")

    return df

def match_programs(spectra_a, spectra_b, n_top=N_TOP_GENES_MATCH):
    """
    Match programs from dataset B to dataset A using cosine similarity
    on the top-N genes by mean TPM weight across both datasets.
    Returns:
      sim_matrix  : (n_programs_a x n_programs_b) cosine similarity matrix
      b_to_a      : array of length n_programs_b giving matched program in A
      match_scores: cosine similarity for each matched pair
    """
    # Common genes
    common_genes = spectra_a.columns.intersection(spectra_b.columns)
    if len(common_genes) == 0:
        raise ValueError("No common genes between datasets for program matching.")

    a = spectra_a[common_genes].values   # programs_a x genes
    b = spectra_b[common_genes].values   # programs_b x genes

    # Restrict to top genes by mean weight across both datasets
    mean_weight = (a.mean(axis=0) + b.mean(axis=0)) / 2
    top_idx     = np.argsort(mean_weight)[::-1][:n_top]
    a_top = a[:, top_idx]
    b_top = b[:, top_idx]

    # Cosine similarity: 1 - cosine distance
    sim = 1 - cdist(a_top, b_top, metric="cosine")   # shape: (n_a, n_b)

    # Hungarian: maximise similarity -> minimise negative similarity
    row_ind, col_ind = linear_sum_assignment(-sim)

    # b_to_a[j] = which program in A is matched to program j in B
    b_to_a      = np.full(b.shape[0], -1, dtype=int)
    match_scores = np.full(b.shape[0], np.nan)
    for r, c in zip(row_ind, col_ind):
        b_to_a[c]      = r
        match_scores[c] = sim[r, c]

    return sim, b_to_a, match_scores


def compute_tf_program_usage(usages, cell_tf, min_cells=MIN_CELLS_PER_TF):
    """
    Compute mean program usage per TF perturbation.
    Returns a DataFrame: rows = TFs, columns = programs.
    Only includes TFs with >= min_cells cells.
    """
    # Align cells: usages index vs cell_tf index
    common_cells = usages.index.intersection(cell_tf.index)
    if len(common_cells) == 0:
        # Try stripping dataset suffix from usage index (e.g. "_0")
        usage_stripped = usages.copy()
        usage_stripped.index = usage_stripped.index.str.rsplit("_", n=1).str[0]
        common_cells = usage_stripped.index.intersection(cell_tf.index)
        if len(common_cells) == 0:
            raise ValueError(
                "No overlapping cell barcodes between usages and h5mu.\n"
                f"  Usage index example:  {usages.index[:3].tolist()}\n"
                f"  Cell TF index example: {cell_tf.index[:3].tolist()}"
            )
        usages = usage_stripped

    u = usages.loc[common_cells]
    t = cell_tf.loc[common_cells]

    tf_usage = {}
    for tf, grp in t.groupby(t):
        if grp.shape[0] < min_cells:
            continue
        tf_usage[tf] = u.loc[grp.index].mean(axis=0).values

    result = pd.DataFrame(tf_usage,
                          index=usages.columns).T   # TFs x programs
    return result


# -
# Load all datasets
# -

all_spectra   = {}   # dataset_name -> spectra DataFrame (programs x genes)
all_usages    = {}   # dataset_name -> usages DataFrame (cells x programs)
all_cell_tf   = {}   # dataset_name -> Series (cell -> gene_name)
all_tf_usage  = {}   # dataset_name -> tf_usage DataFrame (TFs x programs)
all_ensg_map  = {}   # dataset_name -> {ensg: symbol}
all_assoc     = {}   # dataset_name -> TF-to-program mapping

for ds in DATASETS:
    print(f"\nLoading {ds['name']}...")
    spectra, usages          = load_cnmf_files(ds)
    cell_tf                  = load_cell_guide_assignments(ds)
    ensg_map                 = load_ensembl_symbol_map(ds)
    assoc		     = load_perturbation_associations(ds)

    all_spectra[ds["name"]]  = spectra
    all_usages[ds["name"]]   = usages
    all_cell_tf[ds["name"]]  = cell_tf
    all_ensg_map[ds["name"]] = ensg_map
    all_assoc[ds["name"]]    = assoc

    tf_usage = compute_tf_program_usage(usages, cell_tf)
    all_tf_usage[ds["name"]] = tf_usage

    print(f"  [{ds['name']}] TF-program matrix: {tf_usage.shape}")
    print(f"  TFs included: {sorted(tf_usage.index.tolist())}")

dataset_names = list(all_spectra.keys())
n_datasets    = len(dataset_names)

# -
# Figure C1: Program-program similarity matrix between all dataset pairs
# Shows how well each program in dataset A matches to a program in dataset B
# -

if n_datasets >= 2:
    n_pairs = n_datasets * (n_datasets - 1) // 2
    # With 5 datasets there are 10 pairs - a single row would be absurdly
    # wide, so wrap into a grid instead.
    n_cols  = min(5, n_pairs)
    n_rows  = int(np.ceil(n_pairs / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(5.5 * n_cols, 5 * n_rows),
                             squeeze=False)
    ax_idx = 0

    for i in range(n_datasets):
        for j in range(i + 1, n_datasets):
            ds_a = dataset_names[i]
            ds_b = dataset_names[j]

            sim, b_to_a, scores = match_programs(
                all_spectra[ds_a], all_spectra[ds_b]
            )

            ax = axes[ax_idx // n_cols, ax_idx % n_cols]
            sns.heatmap(
                sim,
                ax         = ax,
                cmap       = "Blues",
                vmin       = 0, vmax = 1,
                xticklabels= False,
                yticklabels= False,
                cbar_kws   = {"label": "Cosine similarity"},
            )
            ax.set_title(f"{ds_a}\nvs\n{ds_b}", fontsize = 10, fontweight = "bold")
            ax.set_xlabel(f"Programs ({ds_b})", fontsize = 9)
            ax.set_ylabel(f"Programs ({ds_a})", fontsize = 9)

            # Annotate median match score
            ax.text(0.97, 0.03,
                    f"Median match: {np.nanmedian(scores):.2f}",
                    transform   = ax.transAxes,
                    ha          = "right", va = "bottom",
                    fontsize    = 8, color = "white",
                    fontweight  = "bold")
            ax_idx += 1

    # Hide any unused trailing axes in the grid
    for empty_idx in range(ax_idx, n_rows * n_cols):
        axes[empty_idx // n_cols, empty_idx % n_cols].axis("off")

    fig.suptitle("C1: Program-Program Cosine Similarity Between Datasets\n"
                 f"(top {N_TOP_GENES_MATCH} genes; Hungarian-matched diagonal = best pairing)",
                 fontsize = 12, fontweight = "bold", y = 1.02)
    plt.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "C1_program_similarity_heatmap.png"),
                dpi = 300, bbox_inches = "tight")
    plt.close(fig)
    print("Figure C1 saved.")

else:
    print("Only one dataset loaded - skipping C1 (needs >= 2 datasets).")

# -
# Figure C2: Matched program match-score distribution
# Violin of cosine similarity scores for all matched pairs across dataset pairs
# Higher = programs are more reproducible
# -

if n_datasets >= 2:
    score_records = []
    for i in range(n_datasets):
        for j in range(i + 1, n_datasets):
            ds_a, ds_b = dataset_names[i], dataset_names[j]
            _, _, scores = match_programs(all_spectra[ds_a], all_spectra[ds_b])
            for s in scores:
                if not np.isnan(s):
                    score_records.append({
                        "pair":  f"{ds_a}\nvs\n{ds_b}",
                        "score": s
                    })

    score_df = pd.DataFrame(score_records)

    fig, ax = plt.subplots(figsize=(max(5, 3 * score_df["pair"].nunique()), 5))
    sns.violinplot(data=score_df, x="pair", y="score",
                   palette="Blues", inner="box", ax=ax)
    ax.axhline(0.5, color="grey", linestyle="--", linewidth=0.8, label="0.5 threshold")
    ax.set_xlabel(None)
    ax.set_ylabel("Cosine Similarity (matched program pairs)")
    ax.set_title("C2: Distribution of Best-Match Scores Across Dataset Pairs\n"
                 "Higher = more reproducible programs",
                 fontweight="bold")
    ax.legend(frameon=False)
    plt.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "C2_match_score_distribution.png"),
                dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Figure C2 saved.")

# ------------------------------------------------------------
# Figure C3: TF-program perturbation effects
# Rows = TFs shared across datasets
# Columns = matched programs
# Values = log2FC from perturbation association analysis
# ------------------------------------------------------------

shared_tfs = sorted(
    set.intersection(
        *[set(df["target_name"].unique())
          for df in all_assoc.values()]
    )
)

ref_name = dataset_names[0]
n_programs = all_spectra[ref_name].shape[0]

# program mappings
ds_mappings = {}
for ds_name in dataset_names[1:]:
    _, b_to_a, _ = match_programs(
        all_spectra[ref_name],
        all_spectra[ds_name]
    )
    ds_mappings[ds_name] = b_to_a


def build_tf_program_matrix(ds_name):

    assoc = all_assoc[ds_name]

    mat = np.full(
        (len(shared_tfs), n_programs),
        np.nan
    )

    for tf_idx, tf in enumerate(shared_tfs):

        tf_df = assoc.loc[
            assoc["target_name"] == tf
        ]

        effect = np.full(n_programs, np.nan)

        for _, row in tf_df.iterrows():
            prog = int(row["program_name"])
            effect[prog] = row["log2FC"]

        if ds_name == ref_name:
            mat[tf_idx] = effect

        else:

            b_to_a = ds_mappings[ds_name]

            a_to_b = np.full(
                n_programs,
                -1,
                dtype=int
            )

            for b_idx, a_idx in enumerate(b_to_a):
                if 0 <= a_idx < n_programs:
                    a_to_b[a_idx] = b_idx

            aligned = np.array([
                effect[a_to_b[i]]
                if a_to_b[i] >= 0
                else np.nan
                for i in range(n_programs)
            ])

            mat[tf_idx] = aligned

    return mat


fig, axes = plt.subplots(
    1,
    n_datasets,
    figsize=(4*n_datasets + 1,
             max(6, 0.25*len(shared_tfs))),
    sharey=True,
    squeeze=False
)

all_vals = []

aligned_mats = {}

for ds_name in dataset_names:

    mat = build_tf_program_matrix(ds_name)

    aligned_mats[ds_name] = mat

    all_vals.extend(
        mat[np.isfinite(mat)]
    )

vmax = np.nanpercentile(
    np.abs(all_vals),
    99
)

for idx, ds_name in enumerate(dataset_names):

    ax = axes[0, idx]

    im = ax.imshow(
        aligned_mats[ds_name],
        aspect="auto",
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
        interpolation="nearest"
    )

    ax.set_title(
        ds_name,
        fontsize=10,
        fontweight="bold"
    )

    ax.set_xticks(range(n_programs))
    ax.set_xticklabels(
        [f"P{i+1}" for i in range(n_programs)],
        rotation=90,
        fontsize=6
    )

    if idx == 0:
        ax.set_yticks(range(len(shared_tfs)))
        ax.set_yticklabels(
            shared_tfs,
            fontsize=7
        )

fig.subplots_adjust(right=0.88)

cbar_ax = fig.add_axes(
    [0.90, 0.15, 0.015, 0.7]
)

fig.colorbar(
    im,
    cax=cbar_ax,
    label="Perturbation log2FC"
)

fig.suptitle(
    "C3: TF ? Program Perturbation Effects\n"
    f"Programs aligned to {ref_name}",
    fontsize=11,
    fontweight="bold"
)

fig.savefig(
    os.path.join(
        OUT_DIR,
        "C3_tf_program_log2fc_heatmap.png"
    ),
    dpi=300,
    bbox_inches="tight"
)

plt.close(fig)


# -
# Figure C4: TF of interest - mean usage across matched programs, per dataset
# One column per dataset; rows = matched programs; color = mean usage
# Also shows the non-targeting baseline for reference
# -

def make_tf_effect_plot(
    tf_name,
    dataset_names,
    all_assoc,
    all_spectra,
    ref_name,
    out_path
):

    n_programs = all_spectra[ref_name].shape[0]

    ds_mappings = {}

    for ds_name in dataset_names[1:]:

        _, b_to_a, _ = match_programs(
            all_spectra[ref_name],
            all_spectra[ds_name]
        )

        ds_mappings[ds_name] = b_to_a

    effect_vectors = {}

    for ds_name in dataset_names:

        assoc = all_assoc[ds_name]

        tf_df = assoc.loc[
            assoc["target_name"] == tf_name
        ]

        if len(tf_df) == 0:
            continue

        effect = np.full(
            n_programs,
            np.nan
        )

        for _, row in tf_df.iterrows():

            prog = int(row["program_name"])

            effect[prog] = row["log2FC"]

        if ds_name != ref_name:

            b_to_a = ds_mappings[ds_name]

            a_to_b = np.full(
                n_programs,
                -1,
                dtype=int
            )

            for b_idx, a_idx in enumerate(b_to_a):
                if 0 <= a_idx < n_programs:
                    a_to_b[a_idx] = b_idx

            effect = np.array([
                effect[a_to_b[i]]
                if a_to_b[i] >= 0
                else np.nan
                for i in range(n_programs)
            ])

        effect_vectors[ds_name] = effect

    if len(effect_vectors) == 0:

        print(
            f"{tf_name} not found."
        )

        return

    effect_df = pd.DataFrame(
        effect_vectors,
        index=[
            f"P{i+1}"
            for i in range(n_programs)
        ]
    )

    vmax = np.nanpercentile(
        np.abs(effect_df.values),
        99
    )

    fig, ax = plt.subplots(
        figsize=(
            max(4, n_datasets*2),
            max(6, n_programs*0.2)
        )
    )

    sns.heatmap(
        effect_df,
        cmap="RdBu_r",
        center=0,
        vmin=-vmax,
        vmax=vmax,
        linewidths=0.25,
        ax=ax
    )

    ax.set_title(
        f"{tf_name} perturbation effects",
        fontweight="bold"
    )

    ax.set_xlabel("Dataset")
    ax.set_ylabel("Matched Program")

    fig.suptitle(
        f"C4: {tf_name} ? Program Perturbation Profile",
        fontsize=11,
        fontweight="bold"
    )

    fig.savefig(
        out_path,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)

    print(
        f"C4 saved for {tf_name}"
    )


make_tf_effect_plot(
    tf_name=TF_OF_INTEREST,
    dataset_names=dataset_names,
    all_assoc=all_assoc,
    all_spectra=all_spectra,
    ref_name=dataset_names[0],
    out_path=os.path.join(
        OUT_DIR,
        f"C4_{TF_OF_INTEREST}_program_log2fc.png"
    )
)

# -
# Figure C5: Cross-dataset TF-program correlation
# For each pair of datasets: correlate the TF-program usage profiles
# (one value per TF: Spearman correlation of its program usage vector
# between the two datasets, after program matching)
# Shows whether the same TFs activate the same programs across datasets
# -

if n_datasets >= 2:
    from scipy.stats import spearmanr

    corr_records = []
    ref_name    = dataset_names[0]
    ref_spectra = all_spectra[ref_name]

    # Recompute mappings once
    ds_mappings = {}
    for ds_name in dataset_names[1:]:
        _, b_to_a, _ = match_programs(ref_spectra, all_spectra[ds_name])
        ds_mappings[ds_name] = b_to_a

    n_programs = ref_spectra.shape[0]

    def get_aligned_tf_vec(ds_name, tf_name):
        if tf_name not in all_tf_usage[ds_name].index:
            return None
        raw = all_tf_usage[ds_name].loc[tf_name].values
        if ds_name == ref_name:
            return raw
        b_to_a = ds_mappings[ds_name]
        a_to_b = np.full(n_programs, -1, dtype=int)
        for j, a_i in enumerate(b_to_a):
            if 0 <= a_i < n_programs:
                a_to_b[a_i] = j
        return np.array([raw[a_to_b[i]] if a_to_b[i] >= 0 else np.nan
                         for i in range(n_programs)])

    for i in range(n_datasets):
        for j in range(i + 1, n_datasets):
            ds_a, ds_b = dataset_names[i], dataset_names[j]
            pair_label = f"{ds_a}\nvs\n{ds_b}"

            for tf in shared_tfs:
                va = get_aligned_tf_vec(ds_a, tf)
                vb = get_aligned_tf_vec(ds_b, tf)
                if va is None or vb is None:
                    continue
                # Keep only positions where both are not nan
                mask = ~(np.isnan(va) | np.isnan(vb))
                if mask.sum() < 3:
                    continue
                rho, _ = spearmanr(va[mask], vb[mask])
                corr_records.append({
                    "pair": pair_label,
                    "tf":   tf,
                    "rho":  rho
                })

    corr_df = pd.DataFrame(corr_records)

    fig, ax = plt.subplots(
        figsize=(max(6, 3 * corr_df["pair"].nunique()), 5)
    )
    sns.violinplot(data=corr_df, x="pair", y="rho",
                   palette="Set2", inner="box", ax=ax)
    ax.axhline(0, color="grey", linestyle="--", linewidth=0.8)
    ax.axhline(0.5, color="#00b050", linestyle=":", linewidth=0.8,
               label="- = 0.5")

    # Annotate each violin with median
    for idx, pair in enumerate(corr_df["pair"].unique()):
        med = corr_df.loc[corr_df["pair"] == pair, "rho"].median()
        ax.text(idx, med + 0.03, f"{med:.2f}",
                ha="center", va="bottom", fontsize=8, fontweight="bold")

    ax.set_xlabel(None)
    ax.set_ylabel("Spearman - (TF program-usage profile, A vs B)")
    ax.set_title(
        "C5: Cross-Dataset TF-Program Profile Correlation\n"
        "Each point = one TF; high - = same TF activates same programs in both datasets",
        fontweight="bold"
    )
    ax.legend(frameon=False, fontsize=8)
    plt.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "C5_tf_program_correlation.png"),
                dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Figure C5 saved.")

    # Also save the correlation table
    corr_pivot = corr_df.pivot_table(index="tf", columns="pair", values="rho")
    corr_pivot.to_csv(os.path.join(OUT_DIR, "C5_tf_program_correlation_table.csv"))
    print("C5 correlation table saved.")

print(f"\nAll cNMF figures written to: {OUT_DIR}")