"""Upload PerturbNMF outputs for one Huangfu HUES8 dataset to Synapse.

Usage (from project root):

    SYNAPSE_AUTH_TOKEN=... python upload_to_synapse.py --dataset DE
    SYNAPSE_AUTH_TOKEN=... python upload_to_synapse.py --dataset ESC

Synapse layout (under project syn64423137 → "2026_UTSW" or top-level — see DATASETS):

    PerturbNMF/<dataset>/<run_name>/
        Inference/
            adata/
                cNMF_<K>_2_0.h5mu                    # ×8 K values
            Inference.k_selection_stats.df.npz
        Evaluation/
            <K>_2_0/                                 # ×8 K values
                <K>_perturbation_association_results_all.txt
                <K>_GO_term_enrichment.txt
                <K>_geneset_enrichment.txt
                <K>_trait_enrichment.txt
                <K>_Explained_Variance.txt
                <K>_Explained_Variance_Summary.txt
                <K>_fake_perturbation_association_results.txt
        Plot/
            k_selection/                             # full panel + per-metric panels
            Perturb_gene_200_2_0/
                merged_perturbed_genes_<dataset>_K200.pdf  # 583 MB / 888 MB
        Interpretation/Summary_table/200_2_0/
            cNMF_200_2_0.xlsx                        # Stage 3e

The h5mu files are the largest artifacts (a few GB each × 8 K). Skip uploading
all 8 if storage is a concern; the K=200 h5mu is sufficient for downstream
plotting reproduction. Set --upload-only-k200-h5mu to upload just that one.
"""

import argparse
import os
import sys
from pathlib import Path

DATASETS = {
    "DE": {
        "name": "Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq",
        "run":  "042926_huangfu_de_torchcnmf_KskillA",
    },
    "ESC": {
        "name": "Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq",
        "run":  "042926_huangfu_esc_torchcnmf_KskillA",
    },
}

PROJECT_ROOT = "/cellar/users/aklie/projects/tf_perturb_seq"
SYNAPSE_PROJECT_ID = "syn64423137"
PARENT_FOLDER_NAME = "PerturbNMF"

K_VALUES = [30, 50, 60, 80, 100, 200, 250, 300]


def get_children_map(syn, parent_id):
    return {c["name"]: c["id"] for c in syn.getChildren(parent_id)}


def get_or_create_folder_id(syn, parent_id, name):
    import synapseclient
    children = get_children_map(syn, parent_id)
    if name in children:
        return children[name]
    folder = syn.store(synapseclient.Folder(name=name, parent=parent_id))
    return folder.id


def upload_file(syn, local_path, parent_id, *, dry_run=False):
    """Upload a single file. Skips if a child of the same name already exists."""
    import synapseclient
    name = os.path.basename(local_path)
    children = get_children_map(syn, parent_id)
    if name in children:
        print(f"  [skip] {name} (already on Synapse: {children[name]})")
        return children[name]
    if dry_run:
        print(f"  [dry-run] would upload {name} ({os.path.getsize(local_path):,} B) -> {parent_id}")
        return None
    print(f"  [upload] {name} ({os.path.getsize(local_path):,} B) -> {parent_id}")
    f = syn.store(synapseclient.File(path=local_path, parent=parent_id))
    return f.id


def upload_dir_recursive(syn, local_dir, parent_id, *, dry_run=False, skip_patterns=("logs",)):
    """Upload directory tree recursively. Skip dirs matching skip_patterns; skip 0-byte files."""
    for entry in sorted(Path(local_dir).iterdir()):
        if entry.name in skip_patterns or entry.name.startswith("."):
            continue
        if entry.is_dir():
            sub_id = get_or_create_folder_id(syn, parent_id, entry.name) if not dry_run else f"<folder:{entry.name}>"
            print(f"-> folder {entry.name} ({sub_id})")
            upload_dir_recursive(syn, entry, sub_id, dry_run=dry_run, skip_patterns=skip_patterns)
        else:
            if entry.stat().st_size == 0:
                print(f"  [skip] {entry.name} (0 bytes)")
                continue
            upload_file(syn, str(entry), parent_id, dry_run=dry_run)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", choices=["DE", "ESC"], required=True)
    parser.add_argument("--dry-run", action="store_true",
                        help="Print what would be uploaded; don't actually upload.")
    parser.add_argument("--upload-only-k200-h5mu", action="store_true",
                        help="Skip h5mu files except K=200 (saves ~30 GB upload).")
    parser.add_argument("--skip-merged-pdf", action="store_true",
                        help="Skip the ~600-900 MB merged Perturb_gene PDF.")
    parser.add_argument("--token-env", default="SYNAPSE_AUTH_TOKEN",
                        help="Env var holding the Synapse auth token.")
    args = parser.parse_args()

    token = os.environ.get(args.token_env)
    if not token and not args.dry_run:
        print(f"ERROR: ${args.token_env} not set. Run with --dry-run to preview, "
              f"or export the token first.", file=sys.stderr)
        sys.exit(1)

    ds = DATASETS[args.dataset]
    base = Path(PROJECT_ROOT) / "datasets" / ds["name"] / "PerturbNMF" / "Result" / ds["run"]
    if not base.is_dir():
        print(f"ERROR: run dir not found: {base}", file=sys.stderr)
        sys.exit(1)

    if args.dry_run:
        print(f"[dry-run] no Synapse login")
        syn = None
        perturbnmf_id = "<folder:PerturbNMF>"
        ds_id = f"<folder:{ds['name']}>"
        run_id = f"<folder:{ds['run']}>"
    else:
        import synapseclient
        syn = synapseclient.Synapse()
        syn.login(authToken=token)
        print(f"Logged in as {syn.getUserProfile().userName}")

        perturbnmf_id = get_or_create_folder_id(syn, SYNAPSE_PROJECT_ID, PARENT_FOLDER_NAME)
        ds_id = get_or_create_folder_id(syn, perturbnmf_id, ds["name"])
        run_id = get_or_create_folder_id(syn, ds_id, ds["run"])

    print(f"\n=== Uploading {ds['name']}/{ds['run']} -> Synapse {run_id} ===\n")

    # 1) Inference/adata h5mu files
    inference_id = get_or_create_folder_id(syn, run_id, "Inference") if syn else "<folder:Inference>"
    adata_id = get_or_create_folder_id(syn, inference_id, "adata") if syn else "<folder:adata>"
    adata_dir = base / "Inference" / "adata"
    print(f"-> Inference/adata ({adata_id})")
    for k in K_VALUES:
        h5mu = adata_dir / f"cNMF_{k}_2_0.h5mu"
        if not h5mu.exists():
            print(f"  [missing] {h5mu.name}")
            continue
        if args.upload_only_k200_h5mu and k != 200:
            print(f"  [skip] {h5mu.name} (--upload-only-k200-h5mu)")
            continue
        upload_file(syn, str(h5mu), adata_id, dry_run=args.dry_run)

    # k_selection_stats artifact
    ks_npz = base / "Inference" / "Inference.k_selection_stats.df.npz"
    if ks_npz.exists():
        upload_file(syn, str(ks_npz), inference_id, dry_run=args.dry_run)

    # 2) Evaluation/<K>_2_0/* (perturbation, geneset, GO, trait, explained_variance, fake)
    eval_id = get_or_create_folder_id(syn, run_id, "Evaluation") if syn else "<folder:Evaluation>"
    for k in K_VALUES:
        per_k_local = base / "Evaluation" / f"{k}_2_0"
        if not per_k_local.is_dir():
            print(f"-> [missing] Evaluation/{k}_2_0/")
            continue
        per_k_id = get_or_create_folder_id(syn, eval_id, f"{k}_2_0") if syn else f"<folder:{k}_2_0>"
        print(f"-> Evaluation/{k}_2_0/ ({per_k_id})")
        for f in sorted(per_k_local.glob("*.txt")):
            if f.stat().st_size == 0:
                continue
            upload_file(syn, str(f), per_k_id, dry_run=args.dry_run)

    # 3) Plot/k_selection/* (panels)
    plot_id = get_or_create_folder_id(syn, run_id, "Plot") if syn else "<folder:Plot>"
    ksel_local = base / "Plot" / "k_selection"
    if ksel_local.is_dir():
        ksel_id = get_or_create_folder_id(syn, plot_id, "k_selection") if syn else "<folder:k_selection>"
        print(f"-> Plot/k_selection/ ({ksel_id})")
        for f in sorted(ksel_local.glob("*")):
            if f.is_file() and f.stat().st_size > 0:
                upload_file(syn, str(f), ksel_id, dry_run=args.dry_run)

    # 4) Plot/Perturb_gene_200_2_0/merged_perturbed_genes_*.pdf (and per-target PDFs)
    perturb_local = base / "Plot" / "Perturb_gene_200_2_0"
    if perturb_local.is_dir():
        perturb_id = get_or_create_folder_id(syn, plot_id, "Perturb_gene_200_2_0") if syn else "<folder:Perturb_gene_200_2_0>"
        print(f"-> Plot/Perturb_gene_200_2_0/ ({perturb_id})")
        # Always upload the merged combined PDF (high-value summary)
        for merged in perturb_local.glob("merged_perturbed_genes_*.pdf"):
            if not args.skip_merged_pdf and merged.stat().st_size > 0:
                upload_file(syn, str(merged), perturb_id, dry_run=args.dry_run)
        # Per-target PDFs are optional bulk; uploading 1500-2000 PDFs may be slow.
        # Skipped by default. Uncomment to include.
        # for f in sorted(perturb_local.glob("*.pdf")):
        #     if f.name.startswith("merged_") or f.stat().st_size == 0:
        #         continue
        #     upload_file(syn, str(f), perturb_id, dry_run=args.dry_run)

    # 5) Interpretation/Summary_table/200_2_0/cNMF_200_2_0.xlsx
    summary_local = base / "Interpretation" / "Summary_table" / "200_2_0" / "cNMF_200_2_0.xlsx"
    if summary_local.exists() and summary_local.stat().st_size > 0:
        interp_id = get_or_create_folder_id(syn, run_id, "Interpretation") if syn else "<folder:Interpretation>"
        st_id = get_or_create_folder_id(syn, interp_id, "Summary_table") if syn else "<folder:Summary_table>"
        sub_id = get_or_create_folder_id(syn, st_id, "200_2_0") if syn else "<folder:200_2_0>"
        print(f"-> Interpretation/Summary_table/200_2_0/ ({sub_id})")
        upload_file(syn, str(summary_local), sub_id, dry_run=args.dry_run)

    print("\nDone.")


if __name__ == "__main__":
    main()
