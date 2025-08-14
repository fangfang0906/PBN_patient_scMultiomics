# Figure 1B ####
setwd("/Users/fyan1/Downloads/Projects/Yang_PBN_patient_scRNAseq/")
library(Seurat)
library(tidyr)
library(ggplot2)
library(Matrix)
library(reticulate)
library(parallel)
library(dplyr)
use_python('/Users/fyan1/Library/r-miniconda-arm64/envs/r-reticulate/bin/python')
set.seed(1234)

# load data ####
seu <- readRDS("data/seu.rds")
DimPlot(seu, group.by = 'seurat_clusters',label=T,repel = T,raster=FALSE)
seu$Phase <- factor(seu$Phase,levels=c("G1","S","G2M"))
DimPlot(seu, group.by = 'Phase',cols = c("grey", "orange", "red"))+
  labs(title="Cell cycle phase")
DimPlot(seu, group.by = 'patient',label=F,repel = T)
DimPlot(seu, group.by = 'sample',label=F,repel = T,reduction='umap')
DimPlot(seu, group.by = 'PBN_sensitivity',cols=c("chocolate2","peachpuff4"))

# Figure1C, Barplot showing cell type props ####
library(dplyr)
meta <- seu@meta.data
t <- table(meta$cell_type)
good <- names(t)[t>20]
meta <- meta[meta$cell_type %in% good,]
df <- meta %>% group_by(sample, cell_type) %>% 
  summarise(n = n(), .groups = "drop")
df <- as.data.frame(df)
library(RColorBrewer)  
colors <- c("darkgrey", "#8DD3C7", "#FFFFB3", "#BEBADA",
            "#80B1D3", "#CCEBC5", "#B3DE69", "#BC80BD",
            "#FDB462", "#FB8072", "#FCCDE5")
lvls <- c("UP_R4","UP_R3",
          "UP_S2","UP_S1",
          "P_R10","P_S10",
          "P_R9","P_S9",
          "P_R8","P_S8",
          "P_R6","P_S6",
          "P_R5","P_S5",
          "P_R4","P_S4",
          "P_R7","P_R3",
          "P_S7","P_S3",
          "P_R2","P_S2",
          "P_R1","P_S1")
ggplot(df, aes(x = factor(sample, levels = lvls), y = n, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  ggtitle("Cell Type Counts per Sample") +
  xlab("Sample") + ylab("Number of Cells") +
  scale_fill_manual(values = colors) +
  theme_bw() +
  coord_flip()

# Figure 1D, Dot plot ####
tmp <- seu[,sample(colnames(seu),5000)]
tmp <- tmp[,!(tmp$cell_type %in% c("doublets","Mono",
                                   "other T","other"))]
Idents(tmp) <- tmp$cell_type
levels(tmp) <- c("CD14 Mono","CD16 Mono","DC","NK",'CD8 T',
                 "CD4 T","stem-like B","B")
tmp <- PrepSCTFindMarkers(tmp)
m <- FindAllMarkers(tmp,assay = "SCT", only.pos = T,
                    min.pct = 0.2,min.diff.pct=0.2,
                    logfc.threshold = 1,
                    max.cells.per.ident = 500)
top10 <- m %>% group_by(cluster) %>% top_n(3, avg_log2FC)
lst <- unique(top10$gene)
lst[1:6] <- c("S100A8","S100A12",'CD14','FCGR3A', 'MS4A7','CDKN1C')
lst[22:24] <- c("CCND1","CD79A","MS4A1")
DotPlot(tmp, features = lst) + 
  RotatedAxis()+
  theme(axis.title.y = element_blank(),axis.title.x=element_blank(),
        axis.text.x = element_text(size = 8,angle = 30, vjust = 0.8, 
                                   hjust=1),
        axis.text.y = element_text(size = 10))+
  scale_color_viridis_c()

top10 <- m %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup()
tmp <- seu[,sample(colnames(seu),3000)]
DoHeatmap(tmp, features = top10$gene) + NoLegend()

# Figure 1F ####
library(corrplot)
library(dplyr)
pca <- t(b_seu@reductions$pca@cell.embeddings)
asplit <- split(colnames(b_seu), b_seu$sample)
aggr <- do.call(cbind, lapply(names(asplit), function(x) {
  rowMeans(pca[,asplit[[x]]])
}))
colnames(aggr) <- names(asplit)
tmp <- aggr[,c("P_S8","P_R8","P_S10","P_R10","P_S3","P_S7","P_R3","P_R7")]
cor <- cor(tmp)
corrplot(cor, tl.col = "black")
corrplot(cor, method = "color",
         tl.col = "black", tl.srt = 45, 
         addCoef.col = "black", number.cex = 0.7)