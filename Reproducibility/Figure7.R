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
use_python('/Users/fyan1/Library/r-miniconda-arm64/envs/r-reticulate/bin/python')
set.seed(1234)

tmp <- readRDS("data/PT3_HSPC.rds")
Idents(tmp) <- tmp$new
m <- FindMarkers(tmp,ident.1 = 'stem_C1',
                 ident.2 = "stem_C2",only.pos=T,
                 logfc.threshold = 0.25,
                 max.cells.per.ident =3000,
                 min.pct=0.1
)
m <- m[m$p_val_adj<0.05,]
m2 <- FindMarkers(tmp,ident.1 = 'stem_C2',,
                  only.pos=T,ident.2 = 'stem_C1',
                  logfc.threshold = 0.25,
                  max.cells.per.ident =3000,
                  min.pct=0.1)
m2 <- m2[m2$p_val_adj<0.05,]

res <- runEnrich(gene=rownames(m),thr=0.05,
                 name="up_in_5",
                 db=database)
res2 <- runEnrich(gene=rownames(m2),thr=0.05,
                  name="up_in_6",
                  db=database)

Idents(tmp) <- tmp$new
m <- FindMarkers(tmp,ident.1 = c('stem_C1',"stem_C2"),
                 ident.2 = c('R',"S"),,only.pos=T,
                 logfc.threshold = 0.25,
                 max.cells.per.ident =3000,
                 min.pct=0.1
)
m <- m[m$p_val_adj<0.05,]
m2 <- FindMarkers(tmp,ident.1 = c('R',"S"),
                  only.pos=T,ident.2 = c('stem_C1',"stem_C2"),
                  logfc.threshold = 0.25,
                  max.cells.per.ident =3000,
                  min.pct=0.1)
m2 <- m2[m2$p_val_adj<0.05,]
res <- runEnrich(gene=rownames(m),thr=0.05,
                 name="up_in_stem",
                 db=database)
res2 <- runEnrich(gene=rownames(m2),thr=0.05,
                  name="up_in_committed",
                  db=database)


VlnPlot(tmp,c("CD44","SOX4","TM4SF1","CD24","CD19","CCND1",'MS4A1','CD79A'),ncol = 4,
        group.by = 'new')
FeaturePlot(tmp,c("HDAC3","HDAC4","HDAC6"),ncol = 3)
FeaturePlot(tmp,c("CD44","SOX4","CCND1",'MS4A1'),ncol = 4)
VlnPlot(tmp,c("KRT8","KRT18","EZR","COL1A1", "COL1A2", "COL3A1"),ncol = 3,
        group.by = 'new')

