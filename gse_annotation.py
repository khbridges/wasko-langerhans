import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib

# import from pipeline
from pipeline_functions import training_data_select
from pipeline_functions import viz_training_data
from pipeline_functions import one_hot_encode
from pipeline_functions import cell_type_classifier
from pipeline_functions import process_label

# reading in preprocessed dataset
gse = sc.read('/Users/katebridges/Downloads/GSE142471.h5ad')

# generate 2D embedding
sc.tl.pca(gse, svd_solver='auto')
sc.pl.pca(gse, color='ident')
sc.pp.neighbors(gse)
sc.tl.umap(gse)

# NN BASED CELL TYPE CLASSIFICATION
# training data selection
cell_types = ['Basal', 'Spinous', 'HF/HFSC', 'Fibroblast', 'Myofibroblast', 'Macrophage', 'Dendritic cell',
              'Langerhans', 'T cell', 'Endothelial', 'Skeletal muscle']

# marker genes (master list)
markers = ['Krt14', 'Krt5', 'Mt2', 'Krt15', 'Dnajb1', 'Fosl1', 'Fcgbp',
           'Krt1', 'Krt10', 'Krtdap', 'Mt4', 'Krt77',
           'Krt17', 'Krt79', 'Cd34',
           'Col1a2', 'Col1a1', 'Col3a1', 'Mfap4', 'Mfap5',
           'Acta2', 'Rgs5', 'Tagln',
           'Itgam', 'Cd14', 'Cd68', 'Ccl9', 'Adgre1', 'Atp6v0d2',
           'Cd74', 'Cd80', 'Ccr7', 'Cd209a',
           'Cd207', 'H2-M2', 'Epcam', 'Ltc4s',
           'Cd3g', 'Nkg7', 'Ctla2a', 'Cd3e', 'Cd3d', 'Cd28', 'Cd7',
           'Icam1', 'Icam2', 'Aqp1', 'Pecam1',
           'Tnnc2', 'Acta1', 'Mylpf', 'Myl1', 'Tnnt3', 'Tnni2', 'Ckm', 'Pvalb',
           'Cpz', 'Gpx3', 'Ctsk']

# encoding which marker genes belong to which cell type
celltypes = np.zeros((len(cell_types), len(markers)))
celltypes[0, :7] = [1, 1, 1, 1, 1, 1, 1]  # basal
celltypes[1, 7:12] = [1, 1, 1, 1, 1]  # spinous
celltypes[2, 12:15] = [1, 1, 1]  # HF/HFSC
celltypes[3, 15:20] = [1, 1, 1, 1, 1]  # fibroblast
celltypes[3, 56:] = [1, 1, 1]
celltypes[4, 20:23] = [1, 1, 1]  # myofibroblast
celltypes[5, 23:29] = [1, 1, 1, 1, 1, 0]  # macrophage
celltypes[6, 29:33] = [1, 1, 1, 1]  # dendritic cell
celltypes[7, 33:37] = [1, 1, 1, 1]  # langerhans
celltypes[8, 37:44] = [1, 1, 1, 1, 1, 1, 1]  # T cell
celltypes[9, 44:48] = [1, 1, 1, 1]  # endothelial
celltypes[10, 48:56] = [1, 1, 1, 1, 1, 1, 1, 1]  # skeletal muscle

correct_ind = np.arange(11)

tot_lab, tot_ideal_ind, tot_traindata, tot_testdata = training_data_select(gse, markers, celltypes, cell_types, correct_ind)

cmap = matplotlib.cm.get_cmap('inferno')
viz_training_data(gse, tot_lab, tot_ideal_ind, cell_types, gse.obsm['X_umap'], cmap, '', (6, 5), 0.75)

# neural network training for cell type classification
learning_rate = 0.025  # altering learning rate to change how much neural net can adjust during each training epoch
training_epochs = 500
batch_size = 100
display_step = 5

# using aggregate data for training to bolster cell type abundances in training sets
tot_lab_onehot = one_hot_encode(tot_lab)
all_train_ind = np.array([])
ideal_ = np.argmax(tot_lab_onehot, axis=1)
train_split = 0.5
for k in np.unique(ideal_):
    # randomly select half of each cell type for training, other half goes to validation
    all_ind = np.where(ideal_ == k)[0]
    train_ind = np.random.choice(all_ind, round(train_split*len(all_ind)), replace=False)
    all_train_ind = np.concatenate((all_train_ind, train_ind))

total_predicted_lab, tot_prob, colorm, pred = cell_type_classifier(tot_lab_onehot, tot_traindata,
                                                                   tot_testdata,
                                                                   all_train_ind,
                                                                   learning_rate, training_epochs, batch_size,
                                                                   display_step)

# visualizing confusion matrix, gives context to accuracy of training
fig, ax = plt.subplots()
ax.matshow(colorm, cmap=plt.cm.Blues)
ticks = np.arange(0, 11)
ax.set_xticks(ticks)
ax.set_xticklabels(cell_types, rotation=60)
ax.set_yticks(ticks)
ax.set_yticklabels(cell_types)
ax.set_xlabel('Predicted')
ax.set_ylabel('Actual')
#
for i in range(11):
    for j in range(11):
        c = colorm[j, i]
        ax.text(i, j, str(c), va='center', ha='center')
plt.show()


# processing results & writing to original object
total_lab1, total_prob1 = process_label(tot_prob, tot_lab, total_predicted_lab, tot_ideal_ind, gse, 0.95)

# storing as metadata
gse.obs['ident'] = total_lab1
gse.obs['nn_prob'] = total_prob1

# map results
celltype_dict_ = {0.0: 'Basal', 1.0: 'Spinous', 2.0: 'HF/HFSC', 3.0: 'Fibroblast', 4.0: 'Myofibroblast',
                  5.0: 'Macrophage', 6.0: 'Dendritic cell', 7.0: 'Langerhans', 8.0: 'T cell', 9.0: 'Endothelial',
                  10.0: 'Skeletal muscle', -1.0: 'Poorly classified'}
gse.obs['celltype'] = gse.obs['ident'].map(celltype_dict_)

# write
gse.write('/Users/katebridges/Downloads/GSE142471.h5ad')
