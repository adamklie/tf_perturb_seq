"""Build a per-meta-cluster summary of top genes and top GO terms.

For each consensus (meta) program at k=20, emit:
  - top 10 genes of the meta program (from median meta-spectra)
  - for every composite (member) cNMF program contributing to that meta program:
      * top 10 genes (from its dataset's gene_spectra_score)
      * top 10 GO BP terms (from its dataset's precomputed GO enrichment)

Outputs (TSV) written next to this script:
  - Shared_GO_summary.tsv : one row per meta_cluster_id (wide)
  - Shared_GO_long.tsv    : one row per (meta_cluster_id, member program)
"""

#%%
from pathlib import Path
import pandas as pd

META_K = 60
TOP_N = 10

META_DIR = Path(f"/oak/stanford/groups/engreitz/Users/ymo/Testing_Scripts/Meta_program/Result/meta_consensus/k_{META_K}")
META_SPECTRA = META_DIR / f"meta_spectra.k_{META_K}.median.txt"
ASSIGN = META_DIR / f"shared_vs_specific.k_{META_K}.tsv"

# NOTE: K values here match what was actually used to produce
# shared_vs_specific.k_20.tsv (program_id ranges in that file):
#   Hon_CM 1..80, Huangfu_embryonic 1..50, Huangfu_definitive 1..50.
# This differs from the current scripts/meta_consensus/config.yaml, which has
# Hon_CM k=50 / Huangfu_definitive k=80 swapped — keep this file in sync with
# whatever the meta_consensus run actually used.
DATASETS = {
    "Hon_CM": {
        "k": 80,
        "spectra": Path("/oak/stanford/groups/engreitz/Users/ymo/Project/IGVF_Hon_CM/Result/050626_1M_CM_torch_halsvar_dataloader/Inference/Inference.gene_spectra_score.k_80.dt_2_0.txt"),
        "go": Path("/oak/stanford/groups/engreitz/Users/ymo/Project/IGVF_Hon_CM/Result/050626_1M_CM_torch_halsvar_dataloader/Evaluation/80_2_0/80_GO_term_enrichment.txt"),
        "regulators": Path("/oak/stanford/groups/engreitz/Users/ymo/Project/IGVF_Hon_CM/Result/050626_1M_CM_torch_halsvar_dataloader/Evaluation/Calibration/CRT_pctMito_pctRibo_logGenes_top50_logCounts_logGuides/80_2_0/80_CRT_batch_percent_mito_pct_counts_ribo_log_n_counts_log_total_gene_umis_log_guides_per_cell_CM.txt"),
    },
    "Huangfu_embryonic": {
        "k": 50,
        "spectra": Path("/oak/stanford/groups/engreitz/Users/ymo/Project/IGVF_Huangfu_embryonic-stemcell/Result/Adam_run/Inference/Inference.gene_spectra_score.k_50.dt_2_0.txt"),
        "go": Path("/oak/stanford/groups/engreitz/Users/ymo/Project/IGVF_Huangfu_embryonic-stemcell/Result/Adam_run/Evaluation/50_2_0/50_GO_term_enrichment.txt"),
        "regulators": Path("/oak/stanford/groups/engreitz/Users/ymo/Project/IGVF_Huangfu_embryonic-stemcell/Result/Adam_run/Evaluation/50_2_0/50_perturbation_association_results_all.txt"),
    },
    "Huangfu_definitive": {
        "k": 50,
        "spectra": Path("/oak/stanford/groups/engreitz/Users/ymo/Project/IGVF_Huangfu_definitive-endoderm/Result/Adam_run/Inference/Inference.gene_spectra_score.k_50.dt_2_0.txt"),
        "go": Path("/oak/stanford/groups/engreitz/Users/ymo/Project/IGVF_Huangfu_definitive-endoderm/Result/Adam_run/Evaluation/50_2_0/50_GO_term_enrichment.txt"),
        "regulators": Path("/oak/stanford/groups/engreitz/Users/ymo/Project/IGVF_Huangfu_definitive-endoderm/Result/Adam_run/Evaluation/50_2_0/50_perturbation_association_results_all.txt"),
    },
}

DATASET_ORDER = ["Hon_CM", "Huangfu_embryonic", "Huangfu_definitive"]

OUT_DIR = Path("/oak/stanford/groups/engreitz/Users/ymo/Testing_Scripts/Meta_program/Result/Shared_GO")
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_WIDE = OUT_DIR / f"Shared_GO_summary.k_{META_K}.tsv"
OUT_LONG = OUT_DIR / f"Shared_GO_long.k_{META_K}.tsv"
OUT_XLSX = OUT_DIR / f"Shared_GO_summary.k_{META_K}.xlsx"

SIG_THRESH = 0.05

#%%
def top_n_genes_from_spectra(spectra_df: pd.DataFrame, n: int = TOP_N) -> dict:
    """Return {row_index: [g1, ..., gn]} sorted by descending row value."""
    out = {}
    for idx, row in spectra_df.iterrows():
        out[idx] = row.sort_values(ascending=False).head(n).index.tolist()
    return out


def top_n_go_from_table(go_df: pd.DataFrame, n: int = TOP_N) -> dict:
    """Return {program_name: [t1, ..., tn]} of GO terms with Adjusted P-value < SIG_THRESH,
    sorted by ascending Adjusted P-value. Returns fewer than n if too few significant hits."""
    sig = go_df[go_df["Adjusted P-value"] < SIG_THRESH]
    out = {}
    for prog, sub in sig.groupby("program_name"):
        terms = sub.sort_values("Adjusted P-value", ascending=True).head(n)["Term"].tolist()
        out[prog] = terms
    return out


def top_n_regulators(reg_df: pd.DataFrame, n: int = TOP_N) -> tuple[dict, dict]:
    """Return ({prog: top n positive}, {prog: top n negative}) regulators.

    Filters to adj_pval < SIG_THRESH, then ranks by adj_pval ascending
    within sign(log2FC). Returns fewer than n if too few significant hits.
    """
    sig = reg_df[reg_df["adj_pval"] < SIG_THRESH]
    pos_out, neg_out = {}, {}
    for prog, sub in sig.groupby("program_name"):
        pos = sub[sub["log2FC"] > 0].sort_values("adj_pval", ascending=True).head(n)["target_name"].tolist()
        neg = sub[sub["log2FC"] < 0].sort_values("adj_pval", ascending=True).head(n)["target_name"].tolist()
        pos_out[prog] = pos
        neg_out[prog] = neg
    return pos_out, neg_out


def pad_to_n(items: list, n: int = TOP_N) -> list:
    return items + [""] * (n - len(items)) if len(items) < n else items[:n]


def join_semi(items: list) -> str:
    return ";".join(str(x) for x in items)


#%%
print(f"Loading meta-spectra: {META_SPECTRA}")
meta_spectra = pd.read_csv(META_SPECTRA, sep="\t", index_col=0)
print(f"  shape={meta_spectra.shape} (rows=meta_cluster_id, cols=genes)")
meta_top10 = top_n_genes_from_spectra(meta_spectra, TOP_N)

print(f"Loading cluster assignments: {ASSIGN}")
assign = pd.read_csv(ASSIGN, sep="\t")
print(f"  shape={assign.shape}")
print(f"  category counts:\n{assign['category'].value_counts()}")


#%%
ds_top10_genes: dict = {}      # {(dataset, program_id): [genes]}
ds_top10_go: dict = {}         # {(dataset, program_id): [GO terms]}
ds_top10_pos_reg: dict = {}    # {(dataset, program_id): [positive regulators]}
ds_top10_neg_reg: dict = {}    # {(dataset, program_id): [negative regulators]}

for ds_name in DATASET_ORDER:
    cfg = DATASETS[ds_name]

    print(f"\n[{ds_name}] gene_spectra_score: {cfg['spectra']}")
    spectra = pd.read_csv(cfg["spectra"], sep="\t", index_col=0)
    print(f"  shape={spectra.shape}")
    top_genes = top_n_genes_from_spectra(spectra, TOP_N)
    for prog_id, gene_list in top_genes.items():
        ds_top10_genes[(ds_name, int(prog_id))] = gene_list

    print(f"[{ds_name}] GO enrichment: {cfg['go']}")
    go = pd.read_csv(cfg["go"], sep="\t")
    print(f"  rows={len(go)}, unique programs={go['program_name'].nunique()}")
    top_go = top_n_go_from_table(go, TOP_N)
    # Ensure every program 1..K has an entry, even if filter dropped all rows
    for pid in range(1, cfg["k"] + 1):
        ds_top10_go[(ds_name, pid)] = top_go.get(pid, [])

    print(f"[{ds_name}] regulators: {cfg['regulators']}")
    reg = pd.read_csv(cfg["regulators"], sep="\t")
    print(f"  rows={len(reg)}, unique programs={reg['program_name'].nunique()}")
    pos, neg = top_n_regulators(reg, TOP_N)
    for pid in range(1, cfg["k"] + 1):
        ds_top10_pos_reg[(ds_name, pid)] = pos.get(pid, [])
        ds_top10_neg_reg[(ds_name, pid)] = neg.get(pid, [])


#%%
def member_sort_key(row):
    return (DATASET_ORDER.index(row["dataset"]), row["program_id"])


assign["_sort"] = assign.apply(member_sort_key, axis=1)
assign_sorted = assign.sort_values(["meta_cluster_id", "_sort"]).drop(columns="_sort")

n_members_per_cluster = assign_sorted.groupby("meta_cluster_id").size()
K_MAX = int(n_members_per_cluster.max())
print(f"\nMax members in any cluster: K_MAX={K_MAX}")
print(f"n_members per cluster:\n{n_members_per_cluster.sort_index()}")


#%%
# Build the wide (one row per meta-cluster) and long (one row per member program)
# output tables in a single pass over the cluster assignments.
wide_rows = []
long_rows = []

for cid, group in assign_sorted.groupby("meta_cluster_id"):
    # `members` is the set of composite cNMF programs assigned to this meta-cluster,
    # already ordered (dataset, program_id) thanks to the upstream sort.
    members = group.reset_index(drop=True)
    category = members["category"].iloc[0]  # all members share the same category
    datasets_present = members["dataset"].drop_duplicates().tolist()
    # Re-order to canonical dataset order (Hon_CM, Huangfu_embryonic, Huangfu_definitive)
    # so the `datasets` column is comparable across rows.
    datasets_present_sorted = [d for d in DATASET_ORDER if d in datasets_present]

    # Sanity check: meta_cluster_id from the assignment file must exist as a row
    # in meta_spectra. If not, the two inputs are mismatched.
    if cid not in meta_top10:
        raise KeyError(
            f"meta_cluster_id={cid} not found in meta_spectra rows "
            f"(available: {list(meta_top10.keys())[:5]}...)"
        )

    # -- Build per-member top-10 sets for overlap statistics --
    # member_gene_sets / member_go_sets : ordered list, one set per member
    #   (parallel to member_datasets), used for pairwise Jaccard.
    # per_dataset_*_union : dataset -> union of all that dataset's members'
    #   top-10 sets. Used for the cross-dataset intersection (more lenient than
    #   intersecting every individual member — a term only has to appear in
    #   *some* program of each dataset, not in *all* of them).
    member_gene_sets = []
    member_go_sets = []
    per_dataset_gene_union: dict = {d: set() for d in datasets_present_sorted}
    per_dataset_go_union: dict = {d: set() for d in datasets_present_sorted}
    for _, m in members.iterrows():
        ds, pid = m["dataset"], int(m["program_id"])
        g = set(ds_top10_genes[(ds, pid)])
        t = set(ds_top10_go[(ds, pid)])
        member_gene_sets.append(g)
        member_go_sets.append(t)
        per_dataset_gene_union[ds] |= g
        per_dataset_go_union[ds] |= t

    def mean_pairwise_jaccard(sets, restrict_cross_dataset_idx=None):
        """Average Jaccard over all unordered pairs of `sets`.

        If `restrict_cross_dataset_idx` is given (list of dataset labels
        parallel to `sets`), only pairs from *different* datasets are
        averaged — that's the right summary for "do members from
        different datasets converge on the same biology?".

        Returns NaN if there are no pairs (e.g. single-member cluster, or
        a specific_* cluster when restricting to cross-dataset pairs).
        """
        from itertools import combinations
        pairs = list(combinations(range(len(sets)), 2))
        if restrict_cross_dataset_idx is not None:
            pairs = [(i, j) for i, j in pairs
                     if restrict_cross_dataset_idx[i] != restrict_cross_dataset_idx[j]]
        if not pairs:
            return float("nan")
        vals = []
        for i, j in pairs:
            u = sets[i] | sets[j]
            vals.append(len(sets[i] & sets[j]) / len(u) if u else 0.0)
        return sum(vals) / len(vals)

    # Compute four flavors of mean Jaccard:
    #   GO/genes × (all pairs / cross-dataset pairs only)
    member_datasets = members["dataset"].tolist()
    mean_jacc_go_all = mean_pairwise_jaccard(member_go_sets)
    mean_jacc_go_xds = mean_pairwise_jaccard(member_go_sets, member_datasets)
    mean_jacc_g_all = mean_pairwise_jaccard(member_gene_sets)
    mean_jacc_g_xds = mean_pairwise_jaccard(member_gene_sets, member_datasets)

    # Shared-set logic:
    #   - n_datasets > 1 : intersect the per-dataset unions (cross-dataset shared)
    #   - n_datasets == 1, n_members > 1 : intersect member top-10 sets within
    #     the single dataset (no per-dataset bucketing to consolidate first)
    #   - n_members == 1 : just use the lone member's top-10 as the "shared" set
    if len(per_dataset_go_union) > 1:
        go_inter_ds = set.intersection(*per_dataset_go_union.values())
    elif len(member_go_sets) > 1:
        go_inter_ds = set.intersection(*member_go_sets)
    else:
        go_inter_ds = member_go_sets[0] if member_go_sets else set()

    if len(per_dataset_gene_union) > 1:
        g_inter_ds = set.intersection(*per_dataset_gene_union.values())
    elif len(member_gene_sets) > 1:
        g_inter_ds = set.intersection(*member_gene_sets)
    else:
        g_inter_ds = member_gene_sets[0] if member_gene_sets else set()

    # Assemble the wide-table row. Note the `x == x` trick on Jaccard values:
    # NaN != NaN in Python, so this emits an empty string for NaN (single-dataset
    # cluster) instead of writing the literal "nan" into the TSV.
    row = {
        "meta_cluster_id": cid,
        "category": category,
        "n_members": len(members),
        "n_datasets": len(datasets_present),
        "datasets": ",".join(datasets_present_sorted),
        "meta_top10_genes": join_semi(pad_to_n(meta_top10[cid], TOP_N)),
        "mean_jaccard_GO_all_pairs": round(mean_jacc_go_all, 4),
        "mean_jaccard_GO_cross_dataset_pairs": (
            round(mean_jacc_go_xds, 4) if mean_jacc_go_xds == mean_jacc_go_xds else ""
        ),
        "n_GO_intersect_across_datasets": len(go_inter_ds),
        "GO_intersect_across_datasets": join_semi(sorted(go_inter_ds)),
        "mean_jaccard_genes_all_pairs": round(mean_jacc_g_all, 4),
        "mean_jaccard_genes_cross_dataset_pairs": (
            round(mean_jacc_g_xds, 4) if mean_jacc_g_xds == mean_jacc_g_xds else ""
        ),
        "n_genes_intersect_across_datasets": len(g_inter_ds),
        "genes_intersect_across_datasets": join_semi(sorted(g_inter_ds)),
    }

    # -- Fill the member_<slot>_* columns of the wide row, and emit one
    #    corresponding row into long_rows per real member. --
    # K_MAX is the max member count across ALL clusters, so every wide row
    # has the same number of columns; clusters with fewer members get blank
    # strings in the trailing slots.
    for slot in range(1, K_MAX + 1):
        if slot <= len(members):
            m = members.iloc[slot - 1]
            ds, pid = m["dataset"], int(m["program_id"])
            label = f"{ds}|{pid}"
            # Fail loudly if the program is missing from the precomputed caches.
            # This usually means DATASETS[*]["k"] doesn't match what the
            # meta_consensus run actually used (see top-of-file comment).
            if (ds, pid) not in ds_top10_genes:
                raise KeyError(f"No gene_spectra_score row for {ds}|{pid} — check K in DATASETS")
            if (ds, pid) not in ds_top10_go:
                raise KeyError(f"No GO enrichment for {ds}|{pid} — check K in DATASETS")
            top_genes = ds_top10_genes[(ds, pid)]
            top_go = ds_top10_go[(ds, pid)]
            row[f"member_{slot}_label"] = label
            row[f"member_{slot}_top10_genes"] = join_semi(pad_to_n(top_genes, TOP_N))
            row[f"member_{slot}_top10_GO"] = join_semi(pad_to_n(top_go, TOP_N))

            # Regulators use .get() with [] default because top_n_regulators
            # silently drops programs with no significant hits — but the
            # populate loop earlier pre-fills ds_top10_*_reg for every pid,
            # so this is belt-and-suspenders.
            top_pos = ds_top10_pos_reg.get((ds, pid), [])
            top_neg = ds_top10_neg_reg.get((ds, pid), [])
            long_rows.append({
                "meta_cluster_id": cid,
                "category": category,
                "member_index": slot,
                "dataset": ds,
                "program_id": pid,
                "top10_genes": row[f"member_{slot}_top10_genes"],
                "top10_GO": row[f"member_{slot}_top10_GO"],
                "top10_positive_regulators": join_semi(pad_to_n(top_pos, TOP_N)),
                "top10_negative_regulators": join_semi(pad_to_n(top_neg, TOP_N)),
            })
        else:
            # Pad unused slots with empty strings so all wide rows have equal length.
            row[f"member_{slot}_label"] = ""
            row[f"member_{slot}_top10_genes"] = ""
            row[f"member_{slot}_top10_GO"] = ""

    wide_rows.append(row)


#%%
wide_df = pd.DataFrame(wide_rows)
long_df = pd.DataFrame(long_rows)

print(f"\nWriting {OUT_WIDE}  shape={wide_df.shape}")
wide_df.to_csv(OUT_WIDE, sep="\t", index=False)

print(f"Writing {OUT_LONG}  shape={long_df.shape}")
long_df.to_csv(OUT_LONG, sep="\t", index=False)


#%%
# Excel workbook:
#   - "Summary": wide table, one row per meta-cluster
#   - "Long":    long table, one row per (meta-cluster, member)
#   - "cluster_<id>": per-cluster sheet with the meta-program's top genes at top
#                    and members as rows (top10_genes, top10_GO side by side)
print(f"\nWriting {OUT_XLSX}")
with pd.ExcelWriter(OUT_XLSX, engine="openpyxl") as writer:
    wide_df.to_excel(writer, sheet_name="Summary", index=False)
    long_df.to_excel(writer, sheet_name="Long", index=False)

    for cid in sorted(wide_df["meta_cluster_id"].unique()):
        row = wide_df[wide_df["meta_cluster_id"] == cid].iloc[0]
        members = long_df[long_df["meta_cluster_id"] == cid].sort_values("member_index")

        header = pd.DataFrame({
            "field": [
                "meta_cluster_id", "category",
                "n_members", "n_datasets", "datasets",
                "meta_top10_genes",
            ],
            "value": [
                row["meta_cluster_id"], row["category"],
                row["n_members"], row["n_datasets"], row["datasets"],
                row["meta_top10_genes"],
            ],
        })

        sheet = f"cluster_{cid}"
        header.to_excel(writer, sheet_name=sheet, index=False, startrow=0)
        members[["member_index", "dataset", "program_id",
                 "top10_genes", "top10_GO",
                 "top10_positive_regulators", "top10_negative_regulators"]].to_excel(
            writer, sheet_name=sheet, index=False, startrow=len(header) + 2
        )

print("\nDone.")

# %%
