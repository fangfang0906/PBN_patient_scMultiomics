library(Seurat)
library(monocle)
library(RColorBrewer)

#Extract data, phenotype data, and feature data from the SeuratObject
cds <- readRDS("data/P3.cds")
p1 <- plot_cell_trajectory(cds,color_by = "sample")+
  scale_color_manual(values=c("chocolate1","chocolate4","steelblue","steelblue4"))
p2 <- plot_cell_trajectory(cds,color_by = "Pseudotime")+
  scale_color_gradient(low = "darkblue", high = "yellow")
p3 <- plot_cell_trajectory(cds,color_by = "State")
pData(cds)[,'Treatment_days'] <- as.character(pData(cds)[,'days'])
p4 <- plot_cell_trajectory(cds,color_by = "Treatment_days")+
  scale_color_manual(values=c("#99D8C9","#41AE76","#238B45","#00441B"))
(p1|p2)/(p3|p4)

# Real time vs pseudotime ####
df <- data.frame(pData(cds),"pseudotime"=cds@phenoData@data$Pseudotime)
df$days <- factor(df$days,
                            levels=c("-46","-3","376","394"))
ggplot(df, aes(x = pseudotime, color = days, fill = days)) +
  geom_density(alpha = 0.3) +
  labs(title = "Density of Cells Along Pseudotime by Sample",
       x = "Pseudotime",
       y = "Density") +
  theme_minimal() +
  scale_color_manual(values = c("#99D8C9", "#41AE76", "#238B45", "#00441B")) +
  scale_fill_manual(values = c("#99D8C9", "#41AE76", "#238B45", "#00441B")) +
  theme(legend.title = element_blank())

# pseudotime_heatmap ####
library(viridis)
cds_subset <- cds[sig_gene_names,]
newdata <- data.frame(Pseudotime = seq(min(cds@phenoData@data$Pseudotime), 
                                       max(cds@phenoData@data$Pseudotime),
                                       length.out = 100)) 
tmp <- plot_pseudotime_heatmap(cds[sig_gene_names,],
                               num_clusters = 3,
                               cores = 8,
                               add_annotation_col = newdata,
                               show_rownames = F,
                               hmcols = viridis(256),
                               return_heatmap = T)
#dev.off()

gene_group <- cutree(tmp$tree_row,3)
gene_group <- data.frame(gene_group)
idx <- match(rownames(gene_group),rownames(diff_test_res))
gene_group$pval <- diff_test_res[idx,'pval']
gene_group$qval <- diff_test_res[idx,'qval']

test_gene <- rownames(gene_group)[gene_group$gene_group==3][1:3]
test_gene <- c("HSPA5","DDIT4","NOP56")
plot_genes_in_pseudotime(cds[test_gene,], 
                         color_by = "PBN_sensitivity")
