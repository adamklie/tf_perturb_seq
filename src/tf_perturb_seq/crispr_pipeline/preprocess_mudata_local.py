"""Local-input variant of `preprocess_mudata.py` from Chikara-Takeuchi/energy_dist_TFperturb.

Reads a local inference_mudata.h5mu (no Synapse fetch), runs scanpy preprocessing
+ PCA, builds the gRNA -> cells dict, and writes the four artifacts the pipeline
expects in --output-dir:
  - preprocessed.h5ad      (anndata with obsm["X_pca"])
  - gRNA_dict.pickle       (dict of gRNA_name -> list of cell barcodes)
  - pca_dataframe.pickle   (DataFrame of PCA coords; what util_functions.load_files reads)
  - annotation_table.csv   (columns: guide_id, intended_target_name, type, spacer, [chr/start/end])

Mirrors the upstream `preprocess_mudata.py` behavior except for the input
source (local file path instead of Synapse syn ID) and the config-file
generation (handled by the calling shell wrapper, not here).

Used by scripts/run_energy_distance_pipeline.sh. The upstream Synapse variant
lives at external/energy_dist_TFperturb/preprocess_mudata.py.
"""
import argparse
import os
import pickle

import muon
import numpy as np
import pandas as pd
import scanpy as sc


def get_promoter_name(row):
    if row["type"] == "non-targeting":
        return "non-targeting"
    return f"{row['intended_target_name']}|{row['intended_target_chr']}:{int(row['intended_target_start'])}-{int(row['intended_target_end'])}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mudata-path", required=True)
    ap.add_argument("--output-dir", required=True)
    args = ap.parse_args()

    out = args.output_dir
    os.makedirs(out, exist_ok=True)

    print(f"[preprocess] reading {args.mudata_path}")
    mdata = muon.read_h5mu(args.mudata_path)
    print(mdata)

    rna_key = "gene" if "gene" in mdata.mod else "rna"
    gRNA_key = "gRNA" if "gRNA" in mdata.mod else "guide"
    print(f"[preprocess] modality keys: rna={rna_key}, guide={gRNA_key}")

    adata_exp = mdata[rna_key].copy()
    sc.pp.filter_genes(adata_exp, min_counts=1)
    sc.pp.normalize_total(adata_exp)
    sc.pp.log1p(adata_exp)
    sc.pp.scale(adata_exp)
    print("[preprocess] running PCA (n_comps=50)")
    sc.tl.pca(adata_exp, random_state=0, n_comps=50)

    print("[preprocess] building gRNA -> cells dict")
    adata_g = mdata[gRNA_key]
    g_x = adata_g.X
    try:
        g_x = g_x.tocsc()
    except Exception:
        pass
    gRNA_dict = {}
    for j, gname in enumerate(adata_g.var_names):
        col = g_x[:, j]
        if hasattr(col, "toarray"):
            col = col.toarray().ravel()
        else:
            col = np.asarray(col).ravel()
        cells = list(adata_g.obs_names[col > 0])
        gRNA_dict[gname] = cells

    print("[preprocess] building annotation table")
    # IGVF MuData's guide modality already has `guide_id` as a column (var_names is also guide_id).
    # Don't reset_index — that would try to insert a `guide_id` column that already exists.
    g_var = adata_g.var.copy()
    if "guide_id" not in g_var.columns:
        g_var["guide_id"] = adata_g.var_names
    keep_cols = ["guide_id"]
    for col in ("intended_target_name", "type", "spacer", "intended_target_chr",
                "intended_target_start", "intended_target_end"):
        if col in g_var.columns:
            keep_cols.append(col)
    annotation = g_var[keep_cols].copy()
    if {"intended_target_chr", "intended_target_start", "intended_target_end"} <= set(annotation.columns):
        annotation["intended_target_name"] = annotation.apply(
            lambda r: get_promoter_name(r) if r.get("type") not in ("non-targeting",) else "non-targeting",
            axis=1,
        )
    annotation.to_csv(os.path.join(out, "annotation_table.csv"), index=False)

    print("[preprocess] writing preprocessed.h5ad + gRNA_dict.pickle + pca_dataframe.pickle")
    adata_exp.write_h5ad(os.path.join(out, "preprocessed.h5ad"))
    with open(os.path.join(out, "gRNA_dict.pickle"), "wb") as fh:
        pickle.dump(gRNA_dict, fh)
    pca_df = pd.DataFrame(adata_exp.obsm["X_pca"], index=adata_exp.obs_names)
    pca_df.to_pickle(os.path.join(out, "pca_dataframe.pickle"))

    print("[preprocess] done")


if __name__ == "__main__":
    main()
