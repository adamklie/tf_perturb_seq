"""Core library for meta-consensus across cNMF runs on separate datasets.

Mirrors the consensus logic in torch_cnmf/cnmf.py:1485-1763 but applied
across separately-trained cNMF runs instead of replicates within one run.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score


@dataclass
class DatasetSpec:
    name: str
    inference_dir: str
    k: int
    dt: str  # e.g. "2_0"


def spectra_path(spec: DatasetSpec, input_file_type: str) -> str:
    """Path to Inference.<file_type>.k_<k>.dt_<dt>.txt for a dataset."""
    fname = f"Inference.{input_file_type}.k_{spec.k}.dt_{spec.dt}.txt"
    return os.path.join(spec.inference_dir, fname)


def load_spectra(spec: DatasetSpec, input_file_type: str) -> pd.DataFrame:
    """Load one dataset's gene-spectra file.

    Returns DataFrame indexed by f"{dataset_name}|{program_id}", columns=genes.
    """
    path = spectra_path(spec, input_file_type)
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing spectra file for {spec.name}: {path}")
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = [f"{spec.name}|{i}" for i in df.index]
    return df


def intersect_genes(dfs: List[pd.DataFrame]) -> pd.Index:
    """Common gene columns across all dataset DataFrames (string intersection)."""
    common = dfs[0].columns
    for d in dfs[1:]:
        common = common.intersection(d.columns)
    return common


def stack_programs(dfs: List[pd.DataFrame], common_genes: pd.Index) -> pd.DataFrame:
    """Subset each df to common_genes and concatenate row-wise."""
    aligned = [d.loc[:, common_genes] for d in dfs]
    return pd.concat(aligned, axis=0)


def normalize_rows(S: pd.DataFrame, method: str) -> pd.DataFrame:
    """Row-wise normalization, configurable.

    l2:        unit-norm rows (matches cnmf.py:1557)
    row_zscore: (row - row.mean()) / row.std()
    row_sum:   row / row.sum()  (L1)
    none:      pass through
    """
    X = S.values.astype(float)
    if method == "l2":
        norms = np.linalg.norm(X, axis=1, keepdims=True)
        norms[norms == 0] = 1.0
        N = X / norms
    elif method == "row_zscore":
        mu = X.mean(axis=1, keepdims=True)
        sd = X.std(axis=1, keepdims=True)
        sd[sd == 0] = 1.0
        N = (X - mu) / sd
    elif method == "row_sum":
        s = X.sum(axis=1, keepdims=True)
        s[s == 0] = 1.0
        N = X / s
    elif method == "none":
        N = X
    else:
        raise ValueError(f"Unknown normalization method: {method}")
    return pd.DataFrame(N, index=S.index, columns=S.columns)


def run_kmeans_consensus(
    N: pd.DataFrame,
    k: int,
    n_init: int = 10,
    random_state: int = 1,
) -> tuple[pd.Series, pd.DataFrame, float]:
    """KMeans on normalized program matrix, then per-cluster median.

    Returns: (cluster_labels, meta_spectra, silhouette).
    cluster_labels are 1-indexed to match cnmf.py:1590.
    meta_spectra is the per-cluster median of the (normalized) input rows.
    """
    km = KMeans(n_clusters=k, n_init=n_init, random_state=random_state)
    km.fit(N.values)
    labels = pd.Series(km.labels_ + 1, index=N.index, name="meta_cluster_id")

    meta_spectra = N.groupby(labels).median()

    sil = float(silhouette_score(N.values, labels.values, metric="euclidean"))
    return labels, meta_spectra, sil


def categorize_clusters(
    labels: pd.Series,
    dataset_names: List[str],
) -> pd.DataFrame:
    """Per-cluster summary: counts per dataset and category label.

    category:
      - "shared_all"     if programs from every dataset in cluster
      - "shared_partial" if 2 datasets contribute
      - "specific_<ds>"  if only 1 dataset contributes
    """
    df = pd.DataFrame({
        "program": labels.index,
        "dataset": [p.split("|", 1)[0] for p in labels.index],
        "program_id": [p.split("|", 1)[1] for p in labels.index],
        "meta_cluster_id": labels.values,
    })

    rows = []
    n_ds_total = len(dataset_names)
    for cid, sub in df.groupby("meta_cluster_id"):
        counts = sub["dataset"].value_counts()
        per_ds = {f"n_{ds}": int(counts.get(ds, 0)) for ds in dataset_names}
        n_ds_present = int((counts > 0).sum())
        if n_ds_present == n_ds_total:
            category = "shared_all"
        elif n_ds_present >= 2:
            category = "shared_partial"
        else:
            only_ds = counts.index[0]
            category = f"specific_{only_ds}"
        rows.append({
            "meta_cluster_id": int(cid),
            "n_programs": int(len(sub)),
            "n_datasets": n_ds_present,
            **per_ds,
            "category": category,
            "member_programs": ",".join(sub["program"].tolist()),
        })
    cluster_summary = pd.DataFrame(rows).sort_values("meta_cluster_id").reset_index(drop=True)
    return cluster_summary


def build_shared_vs_specific(
    labels: pd.Series,
    cluster_summary: pd.DataFrame,
) -> pd.DataFrame:
    """Per-input-program table inheriting cluster category."""
    cat_map = dict(zip(cluster_summary["meta_cluster_id"], cluster_summary["category"]))
    rows = []
    for prog, cid in labels.items():
        ds, pid = prog.split("|", 1)
        rows.append({
            "dataset": ds,
            "program_id": pid,
            "meta_cluster_id": int(cid),
            "category": cat_map[int(cid)],
        })
    return pd.DataFrame(rows)


def write_category_counts(
    cluster_summary: pd.DataFrame,
    dataset_names: List[str],
    path: str,
) -> None:
    """Top-level tally to a small text file."""
    cat_counts = cluster_summary["category"].value_counts()
    with open(path, "w") as f:
        f.write(f"# meta-cluster category tally (n clusters per category)\n")
        f.write(f"shared_all\t{int(cat_counts.get('shared_all', 0))}\n")
        f.write(f"shared_partial\t{int(cat_counts.get('shared_partial', 0))}\n")
        for ds in dataset_names:
            key = f"specific_{ds}"
            f.write(f"{key}\t{int(cat_counts.get(key, 0))}\n")
        f.write(f"\n# total meta-clusters: {len(cluster_summary)}\n")
        f.write(f"# total input programs: {int(cluster_summary['n_programs'].sum())}\n")
