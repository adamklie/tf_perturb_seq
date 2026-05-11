"""Adapter for upstream Stage 3e Excel summary (cNMF_compile_excel_table.ipynb).

Wraps the notebook's compile flow as a CLI so it can run via SLURM. Customizes for
the Huangfu HUES8 datasets where:
  - Single-value categorical (`sample == 'all'`) so `Sample = ['all']`
  - Perturbation files are `<K>_perturbation_association_results_<sample>.txt`
  - No categorical association (`--Perform_categorical` was skipped by design)
  - Guide info already embedded in the h5mu (no separate guide_h5ad needed)
"""

import argparse
import os
import sys
import warnings

import muon as mu
import pandas as pd

# Silence noisy mudata FutureWarnings during read
warnings.filterwarnings("ignore", category=FutureWarning, module="mudata")
warnings.filterwarnings("ignore", category=FutureWarning, module="muon")

sys.path.insert(0, "/cellar/users/aklie/projects/tf_perturb_seq/external/PerturbNMF/src")

from Stage3_Interpretation.B_Summarization.src import (
    compile_Program_loading_score_sheet_long,
    compile_Program_loading_score_sheet_flat,
    Compile_GO_sheet,
    Compile_Geneset_sheet,
    Compile_Trait_sheet,
    Compile_Perturbation_sheet,
    Compile_Association_sheet,
    Compile_Explained_variance,
    Compile_Target_Summary_sheet,
    Compile_Summary_sheet,
    load_simple_sheets,
    add_specificity_scores_file,
    check_program_name_match,
)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out_dir", required=True, help="Parent of <run_name>/")
    p.add_argument("--run_name", required=True)
    p.add_argument("--save_path", required=True, help="Intermediate working dir for specificity scores")
    p.add_argument("--K", type=int, default=200)
    p.add_argument("--sel_thresh", type=float, default=2.0)
    p.add_argument("--samples", nargs="+", default=["all"])
    p.add_argument("--perturbation_file_name", default="perturbation_association_results")
    p.add_argument("--non_targeting_key", nargs="+", default=["non-targeting"])
    p.add_argument("--categorical_key", default="sample")
    p.add_argument("--prog_key", default="cNMF")
    p.add_argument("--data_key", default="rna")
    p.add_argument("--guide_targets_key", default="guide_targets")
    p.add_argument("--effect_size", default="log2FC")
    p.add_argument("--adjusted_pval_key", default="Adjusted P-value")
    p.add_argument("--num_gene", type=int, default=300)
    args = p.parse_args()

    k = args.K
    sel_thresh = args.sel_thresh
    sel_str = str(sel_thresh).replace(".", "_")

    output_folder = f"{args.out_dir}/{args.run_name}/Interpretation/Summary_table/{k}_{sel_str}"
    os.makedirs(output_folder, exist_ok=True)
    os.makedirs(args.save_path, exist_ok=True)

    h5mu = f"{args.out_dir}/{args.run_name}/Inference/adata/cNMF_{k}_{sel_str}.h5mu"
    print(f"Reading {h5mu}...")
    mdata = mu.read(h5mu)

    print("Loading per-K eval results...")
    (df_loading_long, df_loading_flat, df_GO, df_Geneset, df_Trait,
     df_Perturbation, df_Association, df_ExplainedVariance,
     df_Perturbation_significant_only) = load_simple_sheets(
        mdata, args.out_dir, args.run_name, k, sel_thresh,
        num_gene=args.num_gene,
        Sample=args.samples,
        perturbation_file_name=args.perturbation_file_name,
        GO_Term_key="Term", GO_Genes_key="Genes",
        Geneset_Term_key="Term", Geneset_Genes_key="Genes",
        Trait_Term_key="Term", Trait_Genes_key="Genes",
        Perturbation_Sample_key="Sample",
    )

    check_program_name_match(
        mdata, prog_key=args.prog_key,
        dataframes=[df_GO, df_Geneset, df_Trait, df_Perturbation,
                    df_Association, df_ExplainedVariance, df_Perturbation_significant_only],
    )

    print("Building Target Summary...")
    perturb_path_base = f"{args.out_dir}/{args.run_name}/Evaluation/{k}_{sel_str}/{k}_{args.perturbation_file_name}"
    df_TargetSummary = Compile_Target_Summary_sheet(
        mdata, perturb_path_base,
        Sample=args.samples,
        categorical_key=args.categorical_key,
        prog_key=args.prog_key,
        data_key=args.data_key,
        guide_targets_key=args.guide_targets_key,
        save_path=args.save_path,
        effect_size=args.effect_size,
    )

    # Workaround for upstream bug: Compile_Summary_sheet at line 773 of
    # Compile_excel_sheet.py does `df['variance_explained'] = df_Explained_Variance`
    # which fails when df_Explained_Variance has multiple columns (ours has
    # VarianceExplained + ProgramID + program_name). Reduce to a Series indexed
    # by program_name aligned to the row order Summary expects.
    # `Compile_Explained_variance` (loaded by load_simple_sheets) does
    # df.set_index('program_name'), so program_name is the INDEX (not a column)
    # and df_ExplainedVariance has remaining columns ['VarianceExplained', 'ProgramID'].
    # Reduce to a 1D array of variance values aligned to mdata['cNMF'].var_names order
    # so `df['variance_explained'] = df_ev_for_summary` doesn't choke on the
    # multi-column DataFrame.
    df_ev_for_summary = df_ExplainedVariance
    if df_ExplainedVariance is not None and "VarianceExplained" in getattr(df_ExplainedVariance, "columns", []):
        programs = [str(p) for p in mdata["cNMF"].var_names]
        ev_index_str = df_ExplainedVariance.index.astype(str)
        df_ev_for_summary = (
            df_ExplainedVariance.assign(_idx=ev_index_str)
            .set_index("_idx")["VarianceExplained"]
            .reindex(programs)
            .values
        )
        print(f"  Reshaped df_ExplainedVariance ({df_ExplainedVariance.shape}) -> 1D array len={len(df_ev_for_summary)} for Summary sheet")

    print("Building Summary sheet...")
    df_Summary = Compile_Summary_sheet(
        mdata, df_GO, df_Geneset, df_Perturbation, df_loading_flat, df_ev_for_summary,
        Sample=args.samples,
        specicicity_path=args.save_path,  # original arg name has typo, kept as-is
        categorical_key=args.categorical_key,
        non_tagerting_key=args.non_targeting_key,
        effect_size=args.effect_size,
        adjusted_pval_key=args.adjusted_pval_key,
    )

    out_xlsx = f"{output_folder}/cNMF_{k}_{sel_str}.xlsx"
    print(f"Writing {out_xlsx}...")
    MAX_ROWS = 1048575
    with pd.ExcelWriter(out_xlsx) as writer:
        df_Summary.to_excel(writer, sheet_name="Summary", index=True)
        df_loading_long.to_excel(writer, sheet_name="Program Loadings", index=True)
        df_TargetSummary.to_excel(writer, sheet_name="Targets Summary", index=True)
        if df_Association is not None:
            df_Association.to_excel(writer, sheet_name="Sample Association", index=True)
        if df_Perturbation is not None:
            combined = []
            for samp in args.samples:
                df_p = add_specificity_scores_file(args.save_path, perturb_path_base, samp)
                df_p["Sample"] = samp
                combined.append(df_p)
            df_perturb_full = pd.concat(combined, ignore_index=True)
            for i in range(0, len(df_perturb_full), MAX_ROWS):
                n = i // MAX_ROWS + 1
                df_perturb_full.iloc[i:i + MAX_ROWS].to_excel(
                    writer, sheet_name=f"Perturbation Association {n}", index=True)
            for i in range(0, len(df_Perturbation_significant_only), MAX_ROWS):
                n = i // MAX_ROWS + 1
                df_Perturbation_significant_only.iloc[i:i + MAX_ROWS].to_excel(
                    writer, sheet_name=f"significant regulators only {n}", index=True)
        if df_Trait is not None:
            df_Trait.to_excel(writer, sheet_name="Trait Enrichment", index=True)
        if df_GO is not None:
            df_GO.to_excel(writer, sheet_name="GO Term Enrichment", index=True)
        if df_Geneset is not None:
            df_Geneset.to_excel(writer, sheet_name="Geneset Enrichment", index=True)

    print(f"Done. {out_xlsx}")


if __name__ == "__main__":
    main()
