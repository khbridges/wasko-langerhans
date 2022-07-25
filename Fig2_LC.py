import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import diffxpy.api as de
from scipy.stats import ranksums

# reading in sc object
gse = sc.read('/Users/katebridges/Downloads/GSE142471.h5ad')
cond_dict = {1.0: 'NW', 2.0: 'W'}

# isolating Langerhans cells for separate analysis
gse_LC = gse[np.where(gse.obs['ident'] == 7)[0], :]
# gse_LC.write('/Users/katebridges/Downloads/gse142471_LC.h5ad')

# DIFFERENTIAL EXPR analysis with diffxpy
test_tt = de.test.t_test(
    data=gse_LC.raw,
    grouping=gse_LC.obs['conditionID']
)

# need to omit NaN rows, order by pvalue, write results to csv
res = test_tt.summary().dropna()
res = res.sort_values(by='pval')
res.to_csv('/Users/katebridges/Downloads/LC_DEres_20200519.csv')

# isolating significant genes (FDR < 0.05, abs(log2FC) > 1)
sig_res = res.iloc[np.where(res['qval'] < 0.05)[0], :]
sig_res = sig_res.iloc[np.where(abs(sig_res['log2fc']) > 1)[0], :]

# norm_sig = gse_LC[:, sig_res.index.values]

# reformat data as pandas dataframe?
df = pd.DataFrame(gse_LC[:, sig_res.index.values].X.T.todense())
df.columns = [gse_LC.obs['conditionID'], gse_LC.obs['sample']]
df.index = sig_res['gene']

# need to format to visualize by sample, condition
network_pal = sns.color_palette("husl", len(np.unique(df.columns.get_level_values('conditionID'))))
network_lut = dict(zip(np.unique(df.columns.get_level_values('conditionID')), network_pal))
network_colors = pd.Series(df.columns.get_level_values('conditionID'), index=df.columns).map(network_lut)

network_pal_ = sns.color_palette("husl", len(np.unique(df.columns.get_level_values('sample'))))
network_lut_ = dict(zip(np.unique(df.columns.get_level_values('sample')), network_pal_))
network_colors_ = pd.Series(df.columns.get_level_values('sample'), index=df.columns).map(network_lut_)

all_colors = pd.DataFrame(network_colors).join(pd.DataFrame(network_colors_))

g = sns.clustermap(df, z_score=0, cmap='seismic', vmin=-10, vmax=10, xticklabels='', yticklabels='',
                   col_colors=network_colors_.values, col_cluster=False)
plt.xlabel('Row z-score')
g.ax_heatmap.set_xlabel('Single Langerhans cells')
g.ax_heatmap.set_ylabel('Significant DE genes')

# create custom legend
samples = ['unwounded_1', 'unwounded_2', 'wounded_1', 'wounded_2', 'wounded_3']
for label in np.unique(df.columns.get_level_values('sample')):
    g.ax_col_dendrogram.bar(0, 0, color=network_lut_[label.astype('int')], label=samples[(label-1).astype('int')], linewidth=0)

l1 = g.ax_col_dendrogram.legend(title='Sample ID', loc="center", ncol=3, bbox_to_anchor=(0.6, 0.8),
                                bbox_transform=plt.gcf().transFigure)
plt.show()

# VISUALIZING expression of genes of interest in LCs before/after wounding
# reading in genes of interest
genes_oi = pd.read_excel('/Users/katebridges/Downloads/LC_Genes-for-plots_ssRNAseq_Haensel.xlsx', sheet_name=1, header=None)

# next need to plot violin plots (LCs, nonwounded vs wounded, with statistics, for top 50 genes)
# can do bootstrapping analysis for p vals and attach in csv for Renee

# figure params
sns.set(rc={'figure.figsize': (5.25, 5)})
sns.set_style("whitegrid")
plt.tight_layout()

cmap = matplotlib.cm.get_cmap('viridis')
my_pal0 = {'Nonwounded': cmap(0.01),
           'Wounded': cmap(0.99)}

violinfig_dir = '/Users/katebridges/Documents/figures/forRenee/violin plots/'
for g in range(genes_oi.shape[0]):
    expr_df = pd.DataFrame({'Condition': gse_LC.obs['conditionID'].map(cond_dict).values,
                            'Normalized expression': np.array(gse_LC[:, genes_oi.iloc[g].values[0]].X.todense()).flatten()})
    ax = sns.violinplot(x="Condition", y="Normalized expression", data=expr_df, palette=my_pal0, fig_size=(5, 5)).set(title=genes_oi.iloc[g].values[0])
    plt.tight_layout()
    plt.savefig(violinfig_dir + genes_oi.iloc[g].values[0] + '.png'), plt.close()

# last step: 2 sided wilcoxon rank sum tests to get stats for LC expr of top 50 priority genes, wounded vs. unwounded
lc_stat = np.zeros((np.unique(genes_oi).shape[0], 2))
for g in range(np.unique(genes_oi).shape[0]):
    expr = np.array(gse_LC[:, np.unique(genes_oi)[g]].X.todense()).flatten()
    stat, pval = ranksums(expr[np.where(gse_LC.obs['conditionID'] == 1.0)[0]], expr[np.where(gse_LC.obs['conditionID'] == 2.0)[0]], alternative='two-sided')
    lc_stat[g, 0] = stat
    lc_stat[g, 1] = pval

# write results to xlsx
excelpath = '/Users/katebridges/Documents/LC_wilcoxon_top50.xlsx'
writer = pd.ExcelWriter(excelpath, engine='xlsxwriter')

d = {'gene': np.unique(genes_oi),
     'statistic': lc_stat[:, 0],
     'pvalue': lc_stat[:, 1]
     }
df = pd.DataFrame(data=d)
df.to_excel(writer)
writer.save()
