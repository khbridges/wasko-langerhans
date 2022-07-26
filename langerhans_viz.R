# reference: https://github.com/saeyslab/nichenetr/blob/master/vignettes/circos.md

library(nichenetr)
library(tidyverse)
library(circlize)

# reading in LR links
circos_links <- read.csv(file ='/Users/katebridges/Documents/LC manuscript/TO_EC_LR.csv')
circos_links = circos_links %>% filter(weight > 0.5)

# taking out a specific cell type
circos_links = circos_links[circos_links$ligand_type != 'Langerhans',]

# giving LR links specific color and order
# "Langerhans" = "dark orange",
grid_col_ligand =c("Macrophage" = "gold",
                   "Keratinocyte" = "lawngreen",
                   "Fibroblast" = "royalblue")
grid_col_target =c(
  "Endothelial" = "tomato")

grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_target = tibble(target_type = grid_col_target %>% names(), color_target_type = grid_col_target)

circos_links = circos_links %>% mutate(receptor = paste(receptor," ")) # extra space: make a difference between a gene as ligand and a gene as target!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
links_circle = circos_links %>% select(ligand, receptor, weight)

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color = circos_links %>% distinct(receptor,color_target_type)
grid_target_color = target_color$color_target_type %>% set_names(target_color$receptor)

grid_col =c(grid_ligand_color,grid_target_color)

# give the option that links in the circos plot will be transparant ~ ligand-target potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 

# order ligands and targets
target_order = circos_links$receptor %>% unique()
# ligand_order = c(CAF_specific_ligands,general_ligands,endothelial_specific_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
ligand_order = circos_links$ligand
order = c(ligand_order,target_order)

# prepare gaps between different segments
width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_target = 15
width_same_cell_same_target_type = 0.5

# rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Langerhans") %>% distinct(ligand) %>% nrow() -1)),
# width_different_cell,
gaps = c(
  # width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Macrophage") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Keratinocyte") %>% distinct(ligand) %>% nrow() -1)), 
  width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Fibroblast") %>% distinct(ligand) %>% nrow() -1)), 
  width_ligand_target,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_type == "Endothelial") %>% distinct(receptor) %>% nrow() -1)),
  width_ligand_target
)

# render the circos plot
# par(mar = c(15, 4, 4, 2) + 0.1)
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow",annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #
