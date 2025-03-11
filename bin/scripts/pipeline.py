#!/usr/bin/env python
# coding: utf-8

# # Set-up

# In[1]:


import time
import pickle as pkl
import pandas as pd
import scanpy as sc
import cupy as cp
import rapids_singlecell as rsc
import warnings
warnings.filterwarnings("ignore")

import matplotlib.pyplot as plt
import seaborn as sns


# In[ ]:


import rmm
from rmm.allocators.cupy import rmm_cupy_allocator

rmm.reinitialize(
    managed_memory=False,  # Allows oversubscription
    pool_allocator=False,  # default is False
    devices=0,  # GPU device IDs to register. By default registers only GPU 0.
)
cp.cuda.set_allocator(rmm_cupy_allocator)


# In[9]:


subset = 1000


# # Load and prepare data

# In[3]:


data_load_start = time.time()


# In[ ]:


get_ipython().run_cell_magic('time', '', 'adata = sc.read_10x_h5("/cellar/users/aklie/data/datasets/tf_perturb_seq/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/Cell_Ranger_Output/filtered_feature_bc_matrix.h5")\nadata.var_names_make_unique()\n')


# In[ ]:


if subset is not None:
    assert isinstance(subset, int)
    sc.pp.subsample(adata, n_obs=1000, random_state=1234)
adata_sub


# In[12]:


get_ipython().run_cell_magic('time', '', 'rsc.get.anndata_to_GPU(adata)\n')


# In[13]:


data_load_time = time.time()
print("Total data load and format time: %s" % (data_load_time - data_load_start))


# # Preprocessing

# In[15]:


preprocess_start = time.time()


# ## Quality control 

# In[16]:


get_ipython().run_cell_magic('time', '', 'rsc.pp.flag_gene_family(adata, gene_family_name="MT", gene_family_prefix="MT")\n')


# In[17]:


get_ipython().run_cell_magic('time', '', 'rsc.pp.flag_gene_family(adata, gene_family_name="RIBO", gene_family_prefix="RPS")\n')


# In[18]:


get_ipython().run_cell_magic('time', '', 'rsc.pp.calculate_qc_metrics(adata, qc_vars=["MT", "RIBO"])\n')


# In[19]:


sc.pl.scatter(adata, x="total_counts", y="pct_counts_MT")
sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts")


# In[27]:


sc.pl.violin(adata, "n_genes_by_counts", jitter=0.4)
sc.pl.violin(adata, "total_counts", jitter=0.4)
sc.pl.violin(adata, "pct_counts_MT", jitter=0.4)


# ## Filter

# In[28]:


get_ipython().run_cell_magic('time', '', 'adata = adata[adata.obs["n_genes_by_counts"] < 5000]\nadata = adata[adata.obs["pct_counts_MT"] < 20]\n')


# In[32]:


get_ipython().run_cell_magic('time', '', 'rsc.pp.filter_genes(adata, min_count=3)\n')


# In[34]:


adata.shape


# In[ ]:


adata.write()


# ## Normalize

# In[35]:


get_ipython().run_cell_magic('time', '', 'rsc.pp.normalize_total(adata, target_sum=1e4)\n')


# In[36]:


get_ipython().run_cell_magic('time', '', 'rsc.pp.log1p(adata)\n')


# ## Variable gene selection

# In[37]:


get_ipython().run_cell_magic('time', '', 'rsc.pp.highly_variable_genes(adata, n_top_genes=5000, flavor="cell_ranger")\n')


# In[38]:


get_ipython().run_cell_magic('time', '', 'adata.raw = adata\n')


# In[39]:


get_ipython().run_cell_magic('time', '', 'rsc.pp.filter_highly_variable(adata)\n')


# In[40]:


get_ipython().run_cell_magic('time', '', 'rsc.pp.regress_out(adata, keys=["total_counts", "pct_counts_MT"])\n')


# ## Scale

# In[44]:


get_ipython().run_cell_magic('time', '', 'rsc.pp.scale(adata, max_value=10)\n')


# In[45]:


get_ipython().run_cell_magic('time', '', 'rsc.get.anndata_to_CPU(adata)\n')


# In[ ]:


adata.write()


# # Clustering and Visualization

# ## Principal component analysis

# In[53]:


get_ipython().run_cell_magic('time', '', 'rsc.tl.pca(adata, n_comps=100)\n')


# In[54]:


sc.pl.pca_variance_ratio(adata, log=True, n_pcs=100)


# ## Computing the neighborhood graph and UMAP

# In[55]:


get_ipython().run_cell_magic('time', '', 'rsc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)\n')


# In[56]:


get_ipython().run_cell_magic('time', '', 'rsc.tl.umap(adata)\n')


# # Clustering

# In[57]:


get_ipython().run_cell_magic('time', '', 'rsc.tl.leiden(adata, resolution=0.6)\n')


# In[59]:


get_ipython().run_cell_magic('time', '', 'sc.pl.umap(adata, color=["leiden"], legend_loc="on data")\n')


# In[ ]:


adata.write()


# # DONE

# ---
