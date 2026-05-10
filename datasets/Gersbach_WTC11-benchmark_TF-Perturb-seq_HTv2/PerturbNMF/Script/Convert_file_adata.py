"""Convert HTv2 inference_mudata.h5mu (gene + guide MuData) to AnnData for cNMF input.

Mirrors the pattern in Hon's PerturbNMF/Script/Convert_file_adata.py — copies the
gene mod as adata, attaches guide_assignment from guide.layers, and stages
guide_names/guide_targets in adata.uns.

Input:  HTv2 cleanser_800 inference MuData (HPC)
Output: AnnData ready to feed torch-cNMF_inference_pipeline.py
"""

import muon as mu

INPUT = "/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/runs/cleanser_800_mito_15pc/pipeline_dashboard/inference_mudata.h5mu"
OUTPUT = "/cellar/users/aklie/projects/tf_perturb_seq/datasets/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/PerturbNMF/Data/inference_mudata_cleaned.h5ad"

mudata = mu.read(INPUT)
adata = mudata["gene"].copy()

adata.obsm["guide_assignment"] = mudata["guide"].layers["guide_assignment"].copy()

adata.uns["guide_names"] = list(mudata["guide"].var["guide_id"])
adata.uns["guide_targets"] = list(mudata["guide"].var["gene_name"])

adata.write(OUTPUT)
print(f"wrote {OUTPUT}")
print(f"adata shape: {adata.shape}")
print(f"guide_assignment shape: {adata.obsm['guide_assignment'].shape}")
print(f"n guide_names: {len(adata.uns['guide_names'])}, n guide_targets: {len(adata.uns['guide_targets'])}")
