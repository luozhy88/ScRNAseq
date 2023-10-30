
library(scater)
library(SingleCellExperiment)
library(org.Hs.eg.db)
library(singleCellTK)
library(AnnotationHub)
library(Seurat)
library(edgeR)
library(scRNAseq)
sce <- readRDS("../01_data_prep/scrna_19sample_sce.RDS")
rownames(sce) = rowData(sce)[,"gene_name"]

# Quality control (using mitochondrial genes).
is.mito <- grepl("^MT-", rownames(sce) )
qcstats <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))
filtered <- quickPerCellQC(qcstats, percent_subsets="subsets_Mito_percent")
sce2 <- sce[, !( filtered$discard) ] # | qcstats$total <500





# Normalization.
sce2 <- logNormCounts(sce2)
meta=colData(sce2)  %>% data.frame()

# Feature selection.
library(scran)
dec <- modelGeneVar(sce2)
hvg <- getTopHVGs(dec, prop=0.1)

# PCA.
library(scater)
set.seed(1234)
sce2 <- runPCA(sce2, ncomponents=10, subset_row=hvg)

# Clustering.
library(bluster)
colLabels(sce2) <- clusterCells(sce2, use.dimred='PCA',
                               BLUSPARAM=NNGraphParam(cluster.fun="louvain"))    

# # Visualization.
# sce2 <- runUMAP(sce2, dimred = 'PCA')
# plotUMAP(sce2, colour_by="label")


# Visualization.
sce2 <- runUMAP(sce2)
plotUMAP(sce2,colorBy="label")


