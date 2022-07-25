import os
import warnings
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import pandas as pd
import scipy
import seaborn as sb

# for scran normalization (future version)
# from rpy2.rinterface import RRuntimeWarning
# from rpy2.robjects import pandas2ri
# from rpy2.robjects.packages import importr

""" Loads sample, performs QC and normalization and concatenates into one anndata object
"""

# Set settings
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(scanpy=True, dpi=80, dpi_save=600, color_map='viridis', vector_friendly=False)
sc.settings.autoshow = False  #set to True if running interactively

fig_dir = 'figures/'
sc.settings.figdir = fig_dir

if not os.path.isdir(fig_dir):
    os.makedirs(fig_dir)

# Print versions
sc.logging.print_version_and_date()
sc.logging.print_versions()

results_path = '/Users/katebridges/Downloads/'
results_file = 'GSE142471.h5ad'

if not os.path.isdir(results_path):
    os.makedirs(results_path)

# Automatically convert rpy2 outputs to pandas dataframes
# pandas2ri.activate()
# warnings.filterwarnings("ignore", category=RRuntimeWarning)

########################################################################################################################
########################################################################################################################
# ###### Quality control and merging of samples ######
# Load files and merge them together
# ### Note: Here is where you specify the path to the folder containing the 10x matrices (output from cellranger).
# These are located in the GEO submission accompanying the paper: GSE142471
paths = ['/Users/katebridges/Downloads//GSE142471/NW_0/',
         '/Users/katebridges/Downloads//GSE142471/NW_1/',
         '/Users/katebridges/Downloads//GSE142471/W_0/',
         '/Users/katebridges/Downloads//GSE142471/W_1/',
         '/Users/katebridges/Downloads//GSE142471/W_2/']

samples = ['NW_0', 'NW_1', 'W_0', 'W_1', 'W_2']

adatas = []
###################
# Set filtering parameters based on previous plots

filt_param = {key: '' for key in samples}
filt_param['NW_0'] = {'min_counts': 200,
                      'max_counts': 45000,
                      'min_genes': 500,
                      'max_genes': 5000}

filt_param['NW_1'] = {'min_counts': 200,
                      'max_counts': 45000,
                      'min_genes': 500,
                      'max_genes': 5000}

filt_param['W_0'] = {'min_counts': 200,
                     'max_counts': 45000,
                     'min_genes': 500,
                     'max_genes': 5000}

filt_param['W_1'] = {'min_counts': 200,
                     'max_counts': 45000,
                     'min_genes': 500,
                     'max_genes': 5000}

filt_param['W_2'] = {'min_counts': 200,
                     'max_counts': 45000,
                     'min_genes': 500,
                     'max_genes': 5000}


for path1, sample in zip(paths, samples):
    adata = sc.read_10x_mtx(path1, var_names='gene_symbols', cache=True)
    adata.var_names_make_unique()  # to make variable names unique

    mito_genes = adata.var_names.str.startswith('mt-')
    # for each cell compute fraction of counts in mito genes vs. all genes
    # the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
    adata.obs['percent_mito'] = adata[:, mito_genes].X.sum(axis=1).A1 / adata.X.sum(axis=1).A1
    # add the total counts per cell as observations-annotation to adata
    adata.obs['n_counts'] = adata.X.sum(axis=1).A1
    adata.obs['n_genes'] = np.sum(adata.X > 0, axis=1).A1

    # QC plots
    # sc.pl.highest_expr_genes(adata, n_top=20)
    # sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'], multi_panel=True)
    # sc.pl.violin(adata, 'n_counts', log=True, cut=0)
    # sc.pl.scatter(adata, 'n_counts', 'n_genes', color='percent_mito')
    # sc.pl.violin(adata, 'percent_mito')
    # plt.axhline(0.10, color='orange')
    #
    # plt.figure()
    # sb.distplot(adata.obs['n_counts'], kde=False, bins=60)
    # plt.axvline(filt_param[sample]['min_counts'])
    # plt.axvline(filt_param[sample]['max_counts'])
    #
    # plt.figure()
    # sb.distplot(adata.obs['n_counts'][adata.obs['n_counts'] < 45000], kde=False, bins=60)
    # plt.axvline(filt_param[sample]['min_counts'])
    #
    # plt.figure()
    # p3 = sb.distplot(adata.obs['n_genes'], kde=False, bins=60)
    # plt.axvline(filt_param[sample]['min_genes'])
    # plt.axvline(filt_param[sample]['max_genes'])
    #
    # plt.figure()
    # sb.distplot(adata.obs['n_genes'][adata.obs['n_genes'] < 6000], kde=False, bins=60)
    # plt.axvline(filt_param[sample]['min_genes'])

    # Filter cells
    print('Filter sample: {}'.format(sample))
    print('Number of cells before filters: {:d}'.format(adata.n_obs))

    adata = adata[adata.obs['percent_mito'] < 0.15, :]
    sc.pp.filter_cells(adata, min_genes=filt_param[sample]['min_genes'])
    sc.pp.filter_cells(adata, max_genes=filt_param[sample]['max_genes'])
    print('Number of cells after filters: {:d}'.format(adata.n_obs))
    adata.var_names_make_unique()

    # ####################### Normalization ################
    # Normalization via built-in functions (option to do scran normalization in future versions)
    # hold on to raw, non-norm counts
    adata.raw = adata

    # total-count normalization
    sc.pp.normalize_total(adata, target_sum=1e4)

    # Replot QC parameters after normalization
    adata.obs['n_counts'] = adata.X.sum(axis=1).A1
    adata.obs['n_genes'] = np.sum(adata.X > 0, axis=1).A1
    # sc.pl.scatter(adata, 'n_counts', 'n_genes', color='percent_mito')
    # sc.pl.scatter(adata, 'size_factors', 'n_counts')
    # sc.pl.scatter(adata, 'size_factors', 'n_genes')

    adatas.append(adata.copy())

del adata
adata = adatas[0].concatenate(adatas[1:], batch_key='sample', batch_categories=samples)
del adatas

print('Final dataset:')
print(adata.obs['sample'].value_counts())

# Re-set figure directory to plot all samples together
fig_dir = 'figures/QC/'
sc.settings.figdir = fig_dir

# Plot post-normalization plots for all samples
# sc.pl.scatter(adata, 'n_counts', 'n_genes', color='sample')
# sc.pl.scatter(adata, 'size_factors', 'n_counts', color='sample')
# sc.pl.scatter(adata, 'size_factors', 'n_genes', color='sample')

# Log-normalize the data
sc.pp.log1p(adata)

########################################################################################################################
# General metrics of data quality
########################################################################################################################
X = adata.layers['counts']

# Get mean read count per cell
# X1 = np.expm1(X)
print(X.sum(axis=1).A1.mean())

# Get mean number of detected genes
X2 = X > 0
print(X2.sum(axis=1).A1.mean())

# FYI - no metadata on gene ids (e.g. 'protein coding') can cause save error

# Write results to file
adata.write(results_path + results_file)
