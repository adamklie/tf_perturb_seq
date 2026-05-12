"""WG1 examples — data QC, transcriptome-wide significance, trans-target overlap.

Run from the jamboree root:

    uv run python working_groups/wg1_data_qc/examples.py

Or open in an editor and execute sections individually. Each section uses only
artifacts already in this folder (no Synapse pull needed).
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd  # noqa: E402

from _lib import (  # noqa: E402
    DATASETS,
    landed_datasets,
    load_edistance_summary_wg1,
    load_qc_summary,
    load_significant_tfs,
    load_tf_cross_lineage,
    load_tf_metadata,
)


# §1 — Per-dataset QC at-a-glance ---------------------------------------------
print("\n§1 QC summary")
qc = load_qc_summary()
print(qc[["dataset_name", "n_cells", "gene_umi_median", "mito_pct_median",
          "guides_per_cell_mean", "intended_auroc"]].to_string(index=False))


# §2 — Transcriptome-wide significance per dataset ----------------------------
# Use `n_sig_distance_gt_NC_max` (calibration-robust) as the headline number;
# `n_sig_pval_lt_0p05` is shown alongside but flagged for the Huangfu runs.
print("\n§2 Energy-distance significance per dataset")
ed = load_edistance_summary_wg1()
cols = ["dataset_id", "n_targeting", "n_sig_distance_gt_NC_max",
        "frac_sig_distance_gt_NC_max", "n_sig_pval_lt_0p05", "calibration_state"]
print(ed[cols].to_string(index=False))


# §3 — Cross-lineage TF classification ---------------------------------------
print("\n§3 Cross-lineage classification counts")
cl = load_tf_cross_lineage()
print(cl["classification"].value_counts().to_string())

print("\n   Convergent-significant TFs (sig in every lineage with data):")
convergent = cl[cl["classification"] == "convergent_significant"]
print(convergent[["gene_symbol", "jaspar_tf_family", "lambert_2018_dbd"]].to_string(index=False))


# §4 — UpSet-ready boolean matrix --------------------------------------------
# Build a TF × dataset membership table that drops directly into upsetplot.
print("\n§4 UpSet boolean matrix (head)")
sig_cols = [c for c in cl.columns if c.startswith("sig_dist_gt_NC_max_")]
upset = cl.set_index("gene_symbol")[sig_cols].fillna(False).astype(bool)
upset.columns = [c.replace("sig_dist_gt_NC_max_", "") for c in upset.columns]
upset = upset[upset.any(axis=1)]
print(f"   {len(upset)} TFs significant in ≥1 lineage")
print(upset.head().to_string())

try:
    from upsetplot import UpSet, from_indicators  # noqa: F401
    print("   upsetplot installed — render with:")
    print("       from upsetplot import UpSet, from_indicators")
    print("       UpSet(from_indicators(upset.columns.tolist(), upset)).plot()")
except ImportError:
    print("   upsetplot not installed; `uv add upsetplot` to render the plot.")


# §5 — Per-dataset top-N by distance, joined with TF metadata ----------------
landed = landed_datasets("wg1_significant_tfs")
print(f"\n§5 Per-dataset top-10 targets by distance_mean (landed: {landed})")
for short in landed:
    df = load_significant_tfs(short)
    targeting = df[df["type"] == "targeting"].sort_values("distance_mean", ascending=False)
    print(f"\n   {short}: {len(df)} rows, {df['sig_distance_gt_NC_max'].sum()} sig")
    print(targeting.head(10)[
        ["gene_symbol", "distance_mean", "distance_rank_targeting",
         "sig_distance_gt_NC_max", "jaspar_tf_family"]
    ].to_string(index=False))


# §6 — Optional plot: bar of significance counts -----------------------------
# Uncomment to render:
# import matplotlib.pyplot as plt
# ax = ed.set_index("dataset_id")[["n_sig_distance_gt_NC_max", "n_sig_pval_lt_0p05"]].plot.bar()
# ax.set_ylabel("n TFs called significant")
# plt.tight_layout()
# plt.savefig("wg1_sig_counts.png", dpi=150)
print("\n(Plot snippets in §6 commented out; uncomment to render.)")
