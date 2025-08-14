# Figure 2 ####
library(Seurat)
library(tidyr)
library(ggplot2)
library(Matrix)
library(reticulate)
library(parallel)
library(dplyr)
library(zoo)
library(ggpubr)
library(pheatmap)
library(WebGestaltR)

# Isolate tumor B cells ####
seu <- readRDS("data/PT3.rds")
DimPlot(seu,group.by = 'tmp',cols=c("lightgrey","steelblue"),
        reduction = 'umap',)+
  labs(title="Isolate Tumor B cells")+NoLegend()

b_seu <- readRDS("data/PT3_b_seu.rds")
DimPlot(b_seu, group.by = 'Treatment_days',
        cols=c("#99D8C9","#41AE76","#238B45","#00441B"))
DimPlot(b_seu, group.by = 'PBN_sensitivity',
        cols=c("chocolate2","peachpuff4"))
b_seu$Phase <- factor(b_seu$Phase,levels=c("G1","S","G2M"))
DimPlot(b_seu, group.by = 'Phase',cols = c("grey", "orange", "red"))+
  labs(title="Cell cycle phase")
DimPlot(b_seu, group.by = 'sample',label=F,repel=T)
DimPlot(b_seu, group.by = 'seurat_clusters',label=T,repel=T)

# Heatmap
library(pheatmap)
m <- do.call(rbind,lapply(c(3,6,5,1,4,7),function(x){
  tmp <- FindMarkers(b_seu,assay="SCT",
                     only.pos = T,
                     ident.1 = x,
                     ident.2 = c(0,2),
                     logfc.threshold = 0.25,
                     min.diff.pct = 0.1)
  tmp$gene <- rownames(tmp)
  tmp$cluster <- x
  tmp
}))
m <- m[m$p_val_adj<0.05,]
subsample <- b_seu[,sample(colnames(b_seu),2000)]
good <- c(0,2,3,6,5,1,4,7)
subsample <- subsample[,subsample$seurat_clusters %in% good]
subsample$seurat_clusters <- as.factor(as.character(subsample$seurat_clusters))

top <- m %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup()
top <- arrange(top,factor(cluster,
                          levels=c(3,6,5,1,7,4)))
gene.use <- c(top$gene,"RAD21")
rna <- subsample@assays$SCT@scale.data[gene.use,]
asplit_cells <- split(colnames(subsample), subsample$seurat_clusters)
n <- 2
lvl <- names(asplit_cells)
means <- do.call(cbind, lapply(lvl, function(x){
  df <- rna[,asplit_cells[[x]]]
  t(apply(df,1,function(x){rollapply(x,n,mean,by=n)}))
}))
seurat_cluster <- unlist(lapply(lvl, function(x) rep(x, length(asplit_cells[[x]])%/%n)))
anno_col <- data.frame(seurat_cluster)
rownames(anno_col) <- colnames(means) <- paste0(seq(1:ncol(means)),colnames(means))
anno_col <- arrange(anno_col,factor(seurat_cluster,
                                    levels=c(0,2,3,6,5,1,7,4)))
meta <- subsample@meta.data
pheatmap(means[,rownames(anno_col)],
         cluster_rows = F,cluster_cols = F,
         show_colnames = F,
         breaks = seq(-2, 2, length.out =101),
         gaps_row = c(10,20,30,40,50),
         scale = "row",annotation_col = anno_col,
         main='DEGs in each cluster',
         fontsize_row = 7)

