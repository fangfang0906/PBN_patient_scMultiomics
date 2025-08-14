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

seu <- readRDS("data/seu.rds")
tmp <- seu[,seu$patient=="P_3" & 
             seu$predicted.celltype.l2 %in% c('HSPC',"B naive","B memory","B intermediate","Plasmablast")]
VlnPlot(tmp,'percent.mt',pt.size=0)

rm(seu)
DefaultAssay(tmp) <- 'prediction.score.celltype.l2'
VlnPlot(tmp,c("HSPC","B intermediate","B naive","B memory"),
        group.by = 'predicted.celltype.l2',pt.size=0.01,ncol=4)

library("WebGestaltR")
database <- c(
  "geneontology_Biological_Process_noRedundant",
  "pathway_KEGG",
  "community-contributed_Hallmark50")
runEnrich <- function(gene,thr,name,db){
  WebGestaltR(enrichMethod="ORA", organism="hsapiens",
              enrichDatabase=db,
              enrichDatabaseType="genesymbol",
              interestGene=gene,
              minNum=10,maxNum = 1000,
              interestGeneType="genesymbol",
              fdrThr = thr,
              referenceSet ="genome_protein-coding",
              referenceGeneType="genesymbol",isOutput = T,
              outputDirectory = "./",
              projectName = name)  
}
DefaultAssay(tmp) <- "RNA"
tmp <- JoinLayers(tmp,layers = 'counts')
tmp <- SCTransform(tmp,vars.to.regress = "percent.mt",
                   vst.flavor = "v1",return.only.var.genes = F)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
tmp <- CellCycleScoring(tmp, s.features = s.genes, 
                        g2m.features = g2m.genes, 
                        set.ident = F)
tmp <- RunPCA(tmp)
tmp <- RunUMAP(tmp, reduction = 'pca', dims = 1:20)
tmp <- FindNeighbors(tmp, dims = 1:20)
tmp <- FindClusters(tmp,resolution = 0.2)
DimPlot(tmp,label=T)

tmp$new <- "stem_C1"
tmp$new[tmp$seurat_clusters %in% c(6)] <- "stem_C2"
tmp$new[tmp$seurat_clusters %in% c(0,2,7,8)] <- "S"
tmp$new[tmp$seurat_clusters %in% c(1,4,3)] <- "R"
saveRDS(tmp,"PT3_HSPC.rds")

tmp <- readRDS("PT3_HSPC.rds")
tmp$new2 <- tmp$new


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
# write.table(m,"up_in_HSPC.csv",row.names = T,col.names = T,quote=F,
#             sep=",")
# write.table(m2,"down_in_HSPC.csv",row.names = T,col.names = T,quote=F,
#             sep=",")
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

plot_10cell <- function(x,y=1.5){
  meta <- tmp@meta.data
  expr <- tmp@assays$SCT@scale.data[x,]
  meta$new2 <- "stem-like B"
  meta$new2[meta$new %in% c("R","S")] <- "committed B"
  df <- do.call(rbind,lapply(unique(meta$new2),function(x){
    sample <- rownames(meta)[meta$new2==x]
    data.frame(sample=x,
               expr=rollapply(expr[sample],10,by=10,mean))
  }))
  df <- data.frame(df)
  # df$sample <- factor(df$sample,levels=c("S","R","HSPC"))
  my_comparisons <- list(c("stem-like B","committed B"))
  ggplot(df, aes(sample, expr, fill = sample)) +
    geom_violin(adjust =1,trim=TRUE, scale = "width")+
    geom_boxplot(width = 0.2)+
    stat_compare_means(comparisons = my_comparisons)+
    geom_jitter(size=0.1)+
    ylab("Relative expression") +xlab("")+theme_bw()+
    ggtitle(x)+NoLegend()
}

plot_10cell("CD44")|plot_10cell("SOX4")|plot_10cell("CCND1")|plot_10cell("MS4A1")

p1 <- plot_10cell("CD44")|plot_10cell("SOX4")|plot_10cell("TM4SF1")|plot_10cell("CD24")
p2 <- plot_10cell("CD19")|plot_10cell("CCND1")|plot_10cell("MS4A1")|plot_10cell("CD79A")
p1/p2

a <- read.delim("Project_up_in_HSPC/enrichment_results_up_in_HSPC.txt")
a_top10 <- a[order(a$FDR), ][1:15, ]
a_top_2 <- a[order(a$enrichmentRatio,decreasing = T), ][1:15, ]
a_top_2 <- a_top_2[!(a_top_2$geneSet %in% a_top10$geneSet),]
p1 <- ggplot(a, aes(x=enrichmentRatio, y=-log10(FDR))) +
  geom_point(aes(color=-log10(FDR), size=enrichmentRatio)) +
  geom_text(data=a_top10, 
            aes(label=description), 
            hjust=0, vjust=0, size=3) +
  geom_text(data=a_top_2, 
            aes(label=description), 
            hjust=0, vjust=0, size=3) +
  scale_color_gradient(low = "grey",high = "purple4") +
  theme_bw() +
  labs(x = "Enrichment Ratio", y = "-log10(FDR)", 
       title = "Upregulated pathways in stem-like malignant population",
       color = "-log10(FDR)")
a <- read.delim("Project_down_in_HSPC/enrichment_results_down_in_HSPC.txt")
a$FDR[a$FDR==0] <- min(a$FDR[a$FDR!=0])*0.01
a_top10 <- a[order(a$FDR), ][1:20, ]
a_top_2 <- a[order(a$enrichmentRatio,decreasing = T), ][1:20, ]
a_top_2 <- a_top_2[!(a_top_2$geneSet %in% a_top10$geneSet),]
p2 <- ggplot(a, aes(x=enrichmentRatio, y=-log10(FDR))) +
  geom_point(aes(color=-log10(FDR), size=enrichmentRatio)) +
  geom_text_repel(data=a_top10, 
                  aes(label=description), 
                  hjust=0, vjust=0, size=3,max.overlaps = 30) +
  geom_text_repel(data=a_top_2, 
                  aes(label=description), 
                  hjust=0, vjust=0, size=3) +
  scale_color_gradient(low = "grey",high = "purple4") +
  theme_bw() +
  labs(x = "Enrichment Ratio", y = "-log10(FDR)", 
       title = "Downregulated pathways in stem-like malignant population",
       color = "-log10(FDR)")
p1|p2

# RPL vs MRPL ####
cyto_ribo <- rownames(m)[grep("^(RPS|RPL)",rownames(m))]
cyto_ribo2 <- rownames(m2)[grep("^(RPS|RPL)",rownames(m2))]
mito_ribo <- rownames(m2)[grep("^(MRPS|MRPL)",rownames(m2))]

gene.use <- c(cyto_ribo,cyto_ribo2,mito_ribo)
rna <- tmp@assays$SCT@scale.data[gene.use,]
asplit_cells <- split(colnames(tmp), tmp$new)
n <- 2
lvl <- names(asplit_cells)
means <- do.call(cbind, lapply(lvl, function(x){
  df <- rna[,asplit_cells[[x]]]
  t(apply(df,1,function(x){rollapply(x,n,mean,by=n)}))
}))
seurat_cluster <- unlist(lapply(lvl, function(x) rep(x, length(asplit_cells[[x]])%/%n)))
anno_col <- data.frame(seurat_cluster)
rownames(anno_col) <- colnames(means) <- paste0(seq(1:ncol(means)),colnames(means))
meta <- tmp@meta.data
pheatmap(means[,rownames(anno_col)],
         cluster_rows = F,cluster_cols = F,
         show_colnames = F,
         breaks = seq(-2, 2, length.out =101),
         gaps_row = c(15,53+15),gaps_col= c(791),
         scale = "row",annotation_col = anno_col,
         main='Cytoplasmic ribosomal genes vs Mitochondrial ribosomal genes',
         fontsize_row = 7)

gene.use <- a[a$description=="Citrate cycle (TCA cycle)",'userId']
gene.use <- unlist(strsplit(gene.use,split=";",fixed=T)[[1]])
rna <- tmp@assays$SCT@scale.data[gene.use,]
asplit_cells <- split(colnames(tmp), tmp$new)
n <- 2
lvl <- names(asplit_cells)
means <- do.call(cbind, lapply(lvl, function(x){
  df <- rna[,asplit_cells[[x]]]
  t(apply(df,1,function(x){rollapply(x,n,mean,by=n)}))
}))
seurat_cluster <- unlist(lapply(lvl, function(x) rep(x, length(asplit_cells[[x]])%/%n)))
anno_col <- data.frame(seurat_cluster)
rownames(anno_col) <- colnames(means) <- paste0(seq(1:ncol(means)),colnames(means))
pheatmap(means[,rownames(anno_col)],
         cluster_rows = F,cluster_cols = F,
         show_colnames = F,
         breaks = seq(-2, 2, length.out =101),
         scale = "row",annotation_col = anno_col,
         main='Citrate cycle (TCA cycle)',
         fontsize_row = 7)

gene.use <- a[a$description=="B cell receptor signaling pathway",'userId']
gene.use <- unlist(strsplit(gene.use,split=";",fixed=T)[[1]])
rna <- tmp@assays$SCT@scale.data[gene.use,]
asplit_cells <- split(colnames(tmp), tmp$new)
n <- 2
lvl <- names(asplit_cells)
means <- do.call(cbind, lapply(lvl, function(x){
  df <- rna[,asplit_cells[[x]]]
  t(apply(df,1,function(x){rollapply(x,n,mean,by=n)}))
}))
seurat_cluster <- unlist(lapply(lvl, function(x) rep(x, length(asplit_cells[[x]])%/%n)))
anno_col <- data.frame(seurat_cluster)
rownames(anno_col) <- colnames(means) <- paste0(seq(1:ncol(means)),colnames(means))
anno_col$B_cell_pathway <- colMeans(means[,rownames(anno_col)])
my_comparisons <- list(c("stem-like B","committed B"))
ggplot(anno_col, aes(seurat_cluster, B_cell_pathway, fill = seurat_cluster)) +
  geom_violin(adjust =1,trim=TRUE, scale = "width")+
  geom_boxplot(width = 0.2)+
  stat_compare_means(comparisons = my_comparisons)+
  geom_jitter(size=0.1)+
  ylab("Relative expression") +xlab("")+theme_bw()+
  ggtitle('B cell signaliing pathway')+NoLegend()

t_cell <- seu

seu <- readRDS("data/seu.rds")
t_cell <- seu[,seu$patient=="P_3" & 
             seu$predicted.celltype.l1 %in% c('CD4 T')]
merge <- merge(tmp,y=t_cell)                      

# inferCNV ####
library(Seurat)
library(dplyr)
library(infercnv)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library("RColorBrewer")
df <- tmp@meta.data
df <- df[,c('seurat_clusters'),drop=FALSE]
write.table(df,file="annotation.txt",quote=F,
            row.names = T,col.names = F,sep="\t")
# Gene_order ####
gtf <- read.delim("infercnv/gtf.txt",header = F)
colnames(gtf) <- c("gene_name","seqnames","start","end")
gtf <- gtf[gtf$gene_name %in% rownames(tmp),]
write.table(gtf,file="gene_order.txt",quote=F,row.names = F,
            col.names = F,sep="\t")
# InferCNV ####
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=tmp@assays$SCT$counts,
                                    annotations_file="annotation.txt",
                                    delim="\t",
                                    gene_order_file="gene_order.txt",
                                    ref_group_names=c("0"))

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,
                             out_dir="PT3_S_control_add_HSPC",
                             cluster_by_groups=T,
                             denoise=TRUE,
                             HMM=F,
                             output_format='pdf')

infercnv_obj = readRDS('P3_7_HMM/run.final.infercnv_obj')
m <- c("1","4","7","5","6","3")
infercnv_obj@observation_grouped_cell_indices <- infercnv_obj@observation_grouped_cell_indices[m]
plot_cnv(infercnv_obj, 
         output_filename='P_3_7', 
         # x.center=1, 
         # plot_chr_scale = T,
         output_format='pdf',
         #useRaster=T,
         title = "InferCNV")


