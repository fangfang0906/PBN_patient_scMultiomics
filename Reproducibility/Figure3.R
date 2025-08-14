library(copykit)
library(BiocParallel)
library(ComplexHeatmap)
library(circlize) 
register(MulticoreParam(progressbar = T, workers = 8), default = T)

# Process
tumor <- readRDS("data/scDNA/P3_scDNA.rds")
tumor <- findAneuploidCells(tumor)
# Mark low-quality cells for filtering
tumor <- findOutliers(tumor,resolution = 0.2,k=10)

# Visualize cells labeled by filter and aneuploid status
tumor_r7 <- tumor[, colData(tumor)$sample_info == "P_R3"]
plotHeatmap(tumor_r7, label = c('outlier', 'is_aneuploid','sample_info'), 
            row_split = 'outlier')

plotHeatmap(tumor, label = c( 'is_aneuploid','sample_info'), 
            row_split = 'sample_info')
plotHeatmap(
  tumor,
  label = c('sample_info'),
  row_split = 'sample_info',
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
)

# Remove cells marked as low-quality from the copykit object
tumor <- tumor[,SummarizedExperiment::colData(tumor)$outlier == FALSE]

# kNN smooth profiles
tumor <- knnSmooth(tumor)

# Create a umap embedding 
tumor <- runUmap(tumor,n_neighbors=30,min_dist=0,seed=42)
plotUmap(tumor,label = 'sample_info')

# Search for the K value that maximizes jaccard similarity for clustering of subclones
# This step is optional. A fixed K value can be provided later to findClusters()
tumor <- findSuggestedK(tumor)
plotSuggestedK(tumor)

# Find clusters of similar copy number profiles and plot the results
# If no k_subclones value is provided, automatically detect it from findSuggestedK()
tumor  <- findClusters(tumor)
my_colors <- c("c1"="darkred", "c3"="steelblue", "c2"="grey")
plotUmap(tumor, label = 'subclones')+ scale_fill_manual(values = my_colors)

# Infer absolute copy numbers
tumor <- calcInteger(tumor, method = 'scquantum', assay = 'smoothed_bincounts')
plotMetrics(tumor, metric = 'ploidy', label = 'ploidy_score')
df <- colData(tumor)
p1 <- plotRatio(tumor,sample_name=rownames(df)[df$sample_info=='P_S3'][1])
p2 <- plotRatio(tumor,sample_name=rownames(df)[df$sample_info=='P_R3'][1])
p1/p2





