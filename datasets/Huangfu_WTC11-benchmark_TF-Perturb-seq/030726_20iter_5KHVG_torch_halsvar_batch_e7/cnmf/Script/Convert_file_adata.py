#%%
import muon as mu 
import scanpy as sc
# %%

mudata= mu.read('/oak/stanford/groups/engreitz/Users/ymo/Project/IGVF_Huangfu_WTC11/Data/inference_mudata.h5mu')
# %%
adata=mudata['gene'].copy()

# %%
adata.obsm['guide_assignment'] = mudata['guide'].layers['guide_assignment'].copy()

# %%
adata.uns['guide_names'] =  list(mudata['guide'].var['guide_id'])
adata.uns['guide_targets'] =  list( mudata['guide'].var['gene_name'])

# %%
adata.write('/oak/stanford/groups/engreitz/Users/ymo/Project/IGVF_Huangfu_WTC11/Data/inference_mudata.h5ad')
# %%
