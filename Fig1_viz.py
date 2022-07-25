import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt

sc.set_figure_params(figsize=(4, 4), fontsize=20)
gse = sc.read('/Users/katebridges/Downloads/GSE142471.h5ad')

# lumping samples by expr condition
conditionID = np.zeros(gse.shape[0])
for g in np.arange(gse.shape[0]):
    if gse.obs['sample'][g].str.contains('NW_'):
        conditionID[g] = 1
    else:
        conditionID[g] = 2

gse.obs['conditionID'] = conditionID

# PLOTTING data in UMAP space by cell type, separated by expr condition
# conditions together - using default colormap for celltypes
lim_celltype = ['Basal', 'Dendritic cell', 'Endothelial', 'Fibroblast', 'Langerhans', 'Skeletal muscle', 'Macrophage',
                'Myofibroblast', 'Spinous', 'HF/HFSC', 'T cell']
sc.pl.umap(gse, color='celltype', s=24, groups=lim_celltype)
plt.xlim([np.min(gse.obsm['X_umap'], axis=0)[0]-1, np.max(gse.obsm['X_umap'], axis=0)[0]+1])
plt.ylim([np.min(gse.obsm['X_umap'], axis=0)[1]-1, np.max(gse.obsm['X_umap'], axis=0)[1]+1])

# highlighting naive samples
sc.pl.umap(gse[np.where(gse.obs['conditionID'] == 1)], color='celltype', s=24, groups=lim_celltype)
plt.xlim([np.min(gse.obsm['X_umap'], axis=0)[0]-1, np.max(gse.obsm['X_umap'], axis=0)[0]+1])
plt.ylim([np.min(gse.obsm['X_umap'], axis=0)[1]-1, np.max(gse.obsm['X_umap'], axis=0)[1]+1])

# highlighting wounded samples
sc.pl.umap(gse[np.where(gse.obs['conditionID'] == 2)], color='celltype', s=24, groups=lim_celltype)
plt.xlim([np.min(gse.obsm['X_umap'], axis=0)[0]-1, np.max(gse.obsm['X_umap'], axis=0)[0]+1])
plt.ylim([np.min(gse.obsm['X_umap'], axis=0)[1]-1, np.max(gse.obsm['X_umap'], axis=0)[1]+1])

# VISUALIZING expression of genes of interest by cell type + condition
cond_dict = {1.0: 'NW', 2.0: 'W'}

# pulling out keratinocytes, fibroblasts, immune cell populations for paired/stacked violin plots (W vs. NW)
idx1 = np.where(gse.obs['ident'] != -1)[0]
idx2 = np.where(gse.obs['ident'] != 9)[0]
idx3 = np.where(gse.obs['ident'] != 10)[0]
idx = np.intersect1d(idx1, np.intersect1d(idx2, idx3))

gse_lim = gse[idx, :]

lim_dict = {0.0: 'Keratinocyte',
            1.0: 'Keratinocyte',
            2.0: 'Keratinocyte',
            3.0: 'Fibroblast',
            4.0: 'Fibroblast',
            5.0: 'Macrophage',
            6.0: 'Dendritic cell',
            7.0: 'Langerhans',
            8.0: 'T cell'}

gse_lim.obs['ident_type'] = gse_lim.obs['ident'].map(lim_dict)

sample_celltype = []
for k in range(gse_lim.shape[0]):
    sample_celltype.append(gse_lim.obs['ident_type'][k] + ' ' + gse_lim.obs['conditionID'].map(cond_dict)[k])

gse_lim.obs['sample_celltype'] = sample_celltype

# genes of interest
fig1 = ['Vegfa', 'Fgf2', 'Fgf1', 'Tnf', 'Tgfb1', 'Angpt2', 'Pdgfb', 'Fn1', 'Dll4']

sc.pl.dotplot(gse_lim, fig1, groupby='sample_celltype', dendrogram=True, swap_axes=True, cmap='Reds', use_raw=False)
sc.pl.stacked_violin(gse_lim, fig1, groupby='sample_celltype', dendrogram=True, swap_axes=True, cmap='Reds', use_raw=False, vmin=-0.5)

# UMAP PLOTS - feature plots of genes of interest separated by expr condition
genes_feature = ['Vegfa', 'Fn1', 'Tnf', 'Ptgs1']

for b in genes_feature:
    for f in np.unique(gse.obs['conditionID']):
        sc.pl.umap(gse[np.where(gse.obs['conditionID'] == f)], color=b, s=24, cmap='Reds', vmax=gse[:, b].X.max(), use_raw=False)
        plt.xlim([np.min(gse.obsm['X_umap'], axis=0)[0]-1, np.max(gse.obsm['X_umap'], axis=0)[0]+1])
        plt.ylim([np.min(gse.obsm['X_umap'], axis=0)[1]-1, np.max(gse.obsm['X_umap'], axis=0)[1]+1])
        plt.title(b + ' ' + cond_dict[f])
