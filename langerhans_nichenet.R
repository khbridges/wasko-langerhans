# loading necessary packages
library(dplyr)
library(Seurat) # v4 incoming?
library(ggplot2)
library(tidyr)
library(readr)
library(pheatmap)
library(tibble)
library(writexl)
library(anndata)
library(readxl)

# install.packages("devtools")
# devtools::install_github("saeyslab/nichenetr")

library(nichenetr)
library(tidyverse)

# converted h5ad to rds with sceasy. would recommend using h5seurat for future analyses
gse <- readRDS("~/Downloads/gse142471.rds")

# creating new Seurat object (version mismatch)
gse_2 <- CreateSeuratObject(counts = gse@assays$RNA)
gse_2@meta.data <- gse@meta.data
gse <- gse_2

# need to create mapping for conditionID & cell type
expr_ID <- c(1, 2)
expr <- c("Unwounded", "Wounded")
type_ID <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, -1)
types <- c('Basal', 'Spinous', 'HF/HFSC', 'Fibroblast', 'Myofibroblast', 'Macrophage', 'Dendritic cell',
'Langerhans', 'T cell', 'Endothelial', 'Skeletal muscle', 'Poorly classified')

expr_df <- data.frame(expr_ID, expr)
type_df <- data.frame(type_ID, types)

gse@meta.data$exprID <- expr_df$expr[match(gse$conditionID, expr_df$expr_ID)]
gse@meta.data$CellType <- type_df$types[match(gse$ident, type_df$type_ID)]

# NicheNet analyses follow prior tutorial: https://github.com/saeyslab/nichenetr/blob/master/vignettes/ligand_activity_geneset.md
# reading in nichenet prior models
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

# converting nichenet prior models from human to mouse
lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()

ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]

weighted_networks_lr = weighted_networks_lr %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()

# perform the nichenet analysis - need to define sender & receiver populations
receiver = "Endothelial"
Idents(object = gse) <- 'CellType'
expressed_genes_receiver = get_expressed_genes(receiver, gse, pct = 0.10)

background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

sender_celltypes = c('Basal', 'Spinous', 'HF/HFSC', 'Fibroblast', 'Myofibroblast', 'Macrophage', 'Dendritic cell',
                     'Langerhans', 'T cell', 'Skeletal muscle')

list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, gse, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

# use DE to pull out genes of interest (sig upregulated after treatment)
seurat_obj_receiver= subset(x = gse, idents = receiver)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["exprID"]])

condition_oi = "Wounded"
condition_reference = "Unwounded" 

DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10) %>% rownames_to_column("gene")

geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

# defining potential ligands based on geneset of interest
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

# perform nichenet analysis - ranking ligands by presence of target genes in geneset of interest
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))

# visualizing best ligand results (keeping all to start)
best_upstream_ligands = ligand_activities %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()

# only interested in LC-specific (& other cell type) signals moving forward...
# best upstream ligands were written to file and manually evaluated by cell type specificity
ligands_sorted <- read.csv(file = '/Users/katebridges/Documents/EC_ligands.csv')
keratin_ligand <- ligands_sorted$Ligands[which(ligands_sorted$Specificity == 'Keratinocyte')]
fib_ligand <- ligands_sorted$Ligands[which(ligands_sorted$Specificity == 'Fibroblast')]
mac_ligand <- ligands_sorted$Ligands[grepl('Macrophage', ligands_sorted$Specificity)]
lc_ligand <- ligands_sorted$Ligands[grepl('Langerhans', ligands_sorted$Specificity)]
# lc_ligand <- c("Mmp9", "Cd320", "Cxcl16", "Itgb2", "Efna1", "Pvr", "Pgf", "Wnt9a", "Arf1", "Pcdh7", "Sema7a", "Ocln", "Flrt3", "Cd24a", "Cxcl1", "Cxcl2", "Vegfa")

DotPlot(gse, features = lc_ligand, cols = "RdYlBu") + RotatedAxis()

# inference & visualization of active target genes
active_ligand_target_links_df = lc_ligand %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 100) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands", "Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
p_ligand_target_network

# receptor for top ligands
lr_network_top = lr_network %>% filter(from %in% lc_ligand & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% lc_ligand & to %in% best_upstream_receptors)
# write.csv(lr_network_top_df_large, '/Users/katebridges/Documents/Fib_EC_LR.csv')

lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()

p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network

# # DE analysis for each sender cell type
# # this uses a new nichenetr function - reinstall nichenetr if necessary!
# gse$celltype <- gse$CellType
# DE_table_all = Idents(gse) %>% levels() %>% intersect(sender_celltypes) %>% lapply(get_lfc_celltype, seurat_obj = gse, condition_colname = "exprID", condition_oi = condition_oi, condition_reference = condition_reference, expression_pct = 0.10) %>% reduce(full_join)
# DE_table_all[is.na(DE_table_all)] = 0
# 
# # Combine ligand activities with DE information
# all_lig <- unique(append(append(append(fib_ligand[1:25], keratin_ligand), mac_ligand), lc_ligand))
# ligand_activities_de = ligand_activities[match(all_lig, ligand_activities$test_ligand),] %>% select(test_ligand, pearson) %>% rename(ligand = test_ligand) %>% left_join(DE_table_all %>% rename(ligand = gene))
# ligand_activities_de[is.na(ligand_activities_de)] = 0
# 
# # make LFC heatmap
# lfc_matrix = ligand_activities_de  %>% select(-ligand, -pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand)
# rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()
# 
# order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
# vis_ligand_lfc = lfc_matrix[order_ligands,]
# 
# colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()
# 
# #vis_ligand_lfc <- vis_ligand_lfc[,8]
# 
# p_ligand_lfc = t(lfc_matrix) %>% make_threecolor_heatmap_ggplot("Prioritized ligands","LFC in Sender", low_color = "midnightblue",mid_color = "white", mid = median(lfc_matrix), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"))
# p_ligand_lfc

