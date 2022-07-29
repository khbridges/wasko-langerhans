# reading in DE results from csv (Python generated)
# annotating wrt GO, KEGG databases; viz IPA results
# makes sure to reset working directory
de_res_ <- read.csv("/Users/katebridges/Downloads/LC_DEres_20200519.csv")

# set up kegg database
# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
BiocManager::install(c("gage","gageData"))
BiocManager::install("org.Mm.eg.db")

library(gage)
# load in libraries to annotate data
library(AnnotationDbi)
library(org.Mm.eg.db)
library(stringr)

run_gsea <- function(de_res) {
  kg.hsa <- kegg.gsets(species="mmu")
  
  # separation of signaling and metabolism pathways
  kegg.sig.gs <- kg.hsa$kg.sets[kg.hsa$sig.idx]
  kegg.met.gs <- kg.hsa$kg.sets[kg.hsa$met.idx]
  
  # set up go database
  go.hs <- go.gsets(species="mouse")
  go.bp.gs <- go.hs$go.sets[go.hs$go.subs$BP]
  go.mf.gs <- go.hs$go.sets[go.hs$go.subs$MF]
  go.cc.gs <- go.hs$go.sets[go.hs$go.subs$CC]
  
  de_res$entrez <- mapIds(org.Mm.eg.db, keys=as.character(de_res$gene), column="ENTREZID", keytype="SYMBOL", multiVals="first")
  de_res$symbol <- as.character(de_res$gene)
  de_res$genename <- mapIds(org.Mm.eg.db, keys=as.character(de_res$gene), column="GENENAME", keytype="SYMBOL", multiVals="first")
  
  de_res.fc <- de_res$log2fc
  names(de_res.fc) <- de_res$entrez
  
  # Run enrichment analysis on all log fc
  fc.kegg.sig.p <- gage(de_res.fc, gsets = kegg.sig.gs)
  fc.kegg.met.p <- gage(de_res.fc, gsets = kegg.met.gs)
  
  fc.go.bp.p <- gage(de_res.fc, gsets = go.bp.gs)
  fc.go.mf.p <- gage(de_res.fc, gsets = go.mf.gs)
  fc.go.cc.p <- gage(de_res.fc, gsets = go.cc.gs)
  
  # convert the kegg results to data frames
  fc.kegg.sig.p.up <- as.data.frame(fc.kegg.sig.p$greater)
  fc.kegg.sig.p.up <- na.omit(fc.kegg.sig.p.up)
  
  fc.kegg.sig.p.down <- as.data.frame(fc.kegg.sig.p$less)
  fc.kegg.sig.p.down <- na.omit(fc.kegg.sig.p.down)
  
  fc.kegg.met.p.up <- as.data.frame(fc.kegg.met.p$greater)
  fc.kegg.met.p.up <- na.omit(fc.kegg.met.p.up)
  
  fc.kegg.met.p.down <- as.data.frame(fc.kegg.met.p$less)
  fc.kegg.met.p.down <- na.omit(fc.kegg.met.p.down)
  
  # convert the go results to data frames
  fc.go.bp.p.up <- as.data.frame(fc.go.bp.p$greater)
  fc.go.bp.p.up <- na.omit(fc.go.bp.p.up)
  fc.go.mf.p.up <- as.data.frame(fc.go.mf.p$greater)
  fc.go.mf.p.up <- na.omit(fc.go.mf.p.up)
  fc.go.cc.p.up <- as.data.frame(fc.go.cc.p$greater)
  fc.go.cc.p.up <- na.omit(fc.go.cc.p.up)
  
  fc.go.bp.p.down <- as.data.frame(fc.go.bp.p$less)
  fc.go.bp.p.down <- na.omit(fc.go.bp.p.down)
  fc.go.mf.p.down <- as.data.frame(fc.go.mf.p$less)
  fc.go.mf.p.down <- na.omit(fc.go.mf.p.down)
  fc.go.cc.p.down <- as.data.frame(fc.go.cc.p$less)
  fc.go.cc.p.down <- na.omit(fc.go.cc.p.down)
  
  fc.go.bp <- rbind(fc.go.bp.p.up[fc.go.bp.p.up$p.val < 0.05,], 
                    fc.go.bp.p.down[fc.go.bp.p.down$p.val < 0.05,])
  fc.go.bp <- fc.go.bp[order(fc.go.bp$p.val),]
  split <- str_split_fixed(rownames(fc.go.bp), " ", 2)
  fc.go.bp$goID <- split[,1]
  fc.go.bp$goTerm <- split[,2]
  
  fc.go.mf <- rbind(fc.go.mf.p.up[fc.go.mf.p.up$p.val < 0.05,], 
                    fc.go.mf.p.down[fc.go.mf.p.down$p.val < 0.05,])
  fc.go.mf <- fc.go.mf[order(fc.go.mf$p.val),]
  split_ <- str_split_fixed(rownames(fc.go.mf), " ", 2)
  fc.go.mf$goID <- split_[,1]
  fc.go.mf$goTerm <- split_[,2]
  
  fc.kegg.sig <- rbind(fc.kegg.sig.p.up[fc.kegg.sig.p.up$p.val < 0.05,], 
                       fc.kegg.sig.p.down[fc.kegg.sig.p.down$p.val < 0.05,])
  fc.kegg.sig <- fc.kegg.sig[order(fc.kegg.sig$p.val),]
  split_0 <- str_split_fixed(rownames(fc.kegg.sig), " ", 2)
  fc.kegg.sig$goID <- split_0[,1]
  fc.kegg.sig$goTerm <- split_0[,2]
  
  fc.kegg.met <- rbind(fc.kegg.met.p.up[fc.kegg.met.p.up$p.val < 0.05,], 
                       fc.kegg.met.p.down[fc.kegg.met.p.down$p.val < 0.05,])
  fc.kegg.met <- fc.kegg.met[order(fc.kegg.met$p.val),]
  split_1 <- str_split_fixed(rownames(fc.kegg.met), " ", 2)
  fc.kegg.met$goID <- split_1[,1]
  fc.kegg.met$goTerm <- split_1[,2]
  
  return(list(fc.go.bp, fc.go.mf, fc.kegg.sig, fc.kegg.met))
}

# generate results LC DEGs
lc_test <- run_gsea(de_res_)
lc_bp <- lc_test[[1]] # GO BP
lc_mf <- lc_test[[2]]
lc_sig <- lc_test[[3]]
lc_met <- lc_test[[4]]


# visualization of GO/KEGG analysis
# bubble plot to summarize GO BP results
library(ggplot2)
library(dplyr)

plot_gsea <- function(output, threshold) {
  output$ind <- seq(-1, -dim(output)[1], by=-1)
  mid <- 0
  ggplot(output[output$p.val < threshold,], 
        aes(x=-log(p.val), y=ind, size=set.size, color=stat.mean)) +
   geom_point(alpha=0.7)+
   scale_size(range = c(1, 12), name="Gene set size")+
    scale_y_continuous(breaks = output[output$p.val < threshold,]$ind, 
                       labels = str_to_sentence(output[output$p.val < threshold,]$goTerm),"") + 
    scale_x_continuous('-log(p-value)') + 
    scale_color_gradient2(midpoint=mid, low="blue", mid="white",
                          high="red", space ="Lab", name="Z-score")
}

plot_gsea(lc_bp, 0.005)

# mapping back to GO terms of interest
term <- go.bp.gs$'GO:0071621 granulocyte chemotaxis'
term_full <- mapIds(org.Mm.eg.db, keys=as.character(term), column="SYMBOL", keytype="ENTREZID", multiVals="first")

# now interested in overlap with LC DEGs
overlap <- intersect(de_res_[de_res_$pval < 0.05,]$gene, term_full)
