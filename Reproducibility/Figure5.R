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
library(ggvenn)
library(WebGestaltR)

# Validate in unpaired samples ####
b_seu <- readRDS("data/b_seu.rds")
paired <- b_seu[,b_seu$patient %in% c("P_9","P_3","P_7")]
unpaired <- b_seu[,!(b_seu$patient %in% c("P_9","P_3","P_7"))]
t <- table(unpaired$sample)
good <- names(t)[t>20]
unpaired <- unpaired[,unpaired$sample %in% good]
Idents(unpaired) <- unpaired$PBN_sensitivity
unpaired <- PrepSCTFindMarkers(unpaired)
m <- FindMarkers(unpaired,assay="SCT",
                 ident.1 = "R",only.pos = F,
                 logfc.threshold = 0,
                 min.diff.pct = 0)
o <- read.delim("results/overlapped_3patients_sig_genes.csv",header = T,
                sep=",")
o <- rownames(o)
summary(m[o,'avg_log2FC']>0)

deg <- FindMarkers(unpaired,assay="SCT",
                   ident.1 = "R",only.pos = T,
                   logfc.threshold = 0.25,
                   min.diff.pct = 0.1)
deg <- rownames(deg)[deg$avg_log2FC>0.5]
deg_res <- runEnrich(gene=deg,thr=0.05,name="unpaired")

m$sig <- rownames(m) %in% o
library(ggrepel)
highlight <- m[m$sig=="TRUE" & m$avg_log2FC>0.5 &
                           -log10(m$p_val_adj)>15,]
highlight$gene <- rownames(highlight)

ggplot(m, aes(x = avg_log2FC, y = -log10(p_val_adj),
              color=sig))+
  
  geom_point(data = subset(m, sig == "FALSE"), color = "lightgrey") +
  geom_point(data = subset(m, sig == "TRUE"), color = "steelblue") +
  theme_bw()+
  geom_vline(xintercept=c(0),linetype=2)+
  ggtitle("Volcano plot")+
  geom_label_repel(data = highlight,
                  aes(label = gene),
                  size = 3,color = "black",
                  box.padding = 0.5,max.overlaps = 10,
                  point.padding = 0.5)
