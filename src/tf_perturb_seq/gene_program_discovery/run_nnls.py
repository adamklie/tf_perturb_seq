#!/usr/bin/env python3
"""
run_nnls.py  -  Fixed-H NMF (NNLS) projection of a reference h5ad dataset
onto cNMF spectra.

Usage:
    python run_nnls.py \\
        --inpath   /path/to/data \\
        --outpath  /path/to/output \\
        --spectra  spectra_file.txt \\
        --refdata  reference_data.h5ad \\
        --outfile  usage_matrix.csv
"""

import argparse
import os
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from sklearn.decomposition import non_negative_factorization as nnf


def parse_args():
    parser = argparse.ArgumentParser(
        description="Project reference data onto cNMF spectra via fixed-H NMF (NNLS)."
    )
    parser.add_argument("--inpath",  required=True, help="Directory containing input files.")
    parser.add_argument("--outpath", required=True, help="Directory for output files.")
    parser.add_argument("--spectra", required=True, help="Filename of cNMF spectra (TSV, K x G).")
    parser.add_argument("--refdata", required=True, help="Filename of reference h5ad dataset.")
    parser.add_argument("--outfile", required=True, help="Output filename for the usage matrix (CSV).")
    parser.add_argument("--max-iter", type=int, default=1000, help="Max NMF iterations (default: 1000).")
    parser.add_argument("--tol",      type=float, default=1e-4, help="NMF convergence tolerance (default: 1e-4).")
    return parser.parse_args()


def main():
    args = parse_args()

    spectra_path = os.path.join(args.inpath, args.spectra)
    refdata_path = os.path.join(args.inpath, args.refdata)
    os.makedirs(args.outpath, exist_ok=True)
    out_path = os.path.join(args.outpath, args.outfile)

    # -- Load cNMF spectra (K x G) --------------------------------------------
    print(f"[1/6] Loading spectra: {spectra_path}")
    spectra = pd.read_csv(spectra_path, sep="\t", index_col=0)

    # -- Load reference h5ad --------------------------------------------------
    print(f"[2/6] Loading reference data: {refdata_path}")
    adata = sc.read_h5ad(refdata_path)

    # -- Gene name harmonisation -----------------------------------------------
    print("[3/6] Harmonising gene names ...")
    if "feature_name" in adata.var.columns:
        adata.var_names = adata.var["feature_name"].values
    adata.var_names_make_unique()

    shared = sorted(set(spectra.columns) & set(adata.var_names))
    print(
        f"  cNMF genes   : {len(spectra.columns)}\n"
        f"  Reference genes: {adata.n_vars}\n"
        f"  Shared genes : {len(shared)}"
    )
    if len(shared) == 0:
        raise ValueError("No shared genes between spectra and reference data. Check gene name formats.")

    # -- Subset & align -------------------------------------------------------
    print("[4/6] Subsetting to shared genes ...")
    H_fixed   = spectra[shared].values.astype(np.float32)   # (K, G_shared)
    adata_sub = adata[:, shared].copy()
    print(f"  H shape: {H_fixed.shape}  |  X shape: {adata_sub.shape}")

    # -- Normalize (TPM normalization on raw counts) ------------------
    print("[5/6] Normalising (TPM) ...")
    if "counts" in adata_sub.layers:
        X = (adata_sub.layers["counts"].toarray().astype(np.float32)
             if sp.issparse(adata_sub.layers["counts"])
             else np.array(adata_sub.layers["counts"], dtype=np.float32))
        print("  Using 'counts' layer.")
    else:
        X = (adata_sub.X.toarray().astype(np.float32)
             if sp.issparse(adata_sub.X)
             else np.array(adata_sub.X, dtype=np.float32))
        print("  WARNING: 'counts' layer not found, falling back to .X - confirm this is raw counts.")

    # TPM normalization: scale each cell to 1e6 total counts
    cell_totals = X.sum(axis=1, keepdims=True)
    cell_totals[cell_totals == 0] = 1.0          # avoid divide-by-zero for empty cells
    Xn = (X / cell_totals) * 1e6
    print(f"  Xn range: [{Xn.min():.3f}, {Xn.max():.3f}]")

    # -- Fixed-H NMF ----------------------------------------------------------
    print("[6/6] Running fixed-H NMF ...")
    U, _, n_iter = nnf(
        Xn,
        n_components=H_fixed.shape[0],
        H=H_fixed,
        update_H=False,
        init="custom",
        solver="cd",
        beta_loss="frobenius",
        max_iter=args.max_iter,
        tol=args.tol,
    )
    print(f"  U shape: {U.shape}  |  converged in {n_iter} iterations")

    # -- Reconstruction quality -----------------------------------------------
    R = Xn - U @ H_fixed
    rel_sse  = np.sum(R ** 2) / np.sum(Xn ** 2)
    rel_frob = np.linalg.norm(R, "fro") / np.linalg.norm(Xn, "fro")
    print(f"  Relative SSE       : {rel_sse:.4f}")
    print(f"  Relative Frobenius : {rel_frob:.4f}")

    # -- Save usage matrix ----------------------------------------------------
    usage_df = pd.DataFrame(
        U,
        index=adata.obs_names,
        columns=[str(x) for x in spectra.index],
    )
    usage_df.to_csv(out_path)
    print(f"\nUsage matrix written to: {out_path}")


if __name__ == "__main__":
    main()