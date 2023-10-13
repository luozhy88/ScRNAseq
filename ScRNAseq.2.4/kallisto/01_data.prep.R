

library(SingleCellExperiment)
library(rtracklayer)
library(Seurat)
library(SeuratDisk)
library(tibble)
library(dplyr)


## format convert
# Convert("combined_matrix.h5ad", dest = "h5seurat", overwrite = TRUE)
## Seurat format
sixhosp_21subjects_star <- LoadH5Seurat("/data3/zhiyu/pipelines/scRNA/nfcore_scRNA_sixhosp_21subjects_kallisto/output/kallisto/mtx_conversions/combined_matrix.h5seurat")
## SCE format, the best format for single cell data analysis
sce <- as.SingleCellExperiment(sixhosp_21subjects_star)

## we need to annotate the gene code using the gtf for alignment (better to use the same GTF here as used for alignment)
gene.data <- rtracklayer::import("~/data/Projects/YNYK/sixhop/sixhosp_nc_scd1_amci_ad_shjw_ad_version2/Singlecell_Sixhop/scRNA/20231009/nfcore_scRNA_sixhosp_21subjects_kallisto_ge/input/Gencode/gencode.v44.primary_assembly.annotation.gtf")
# Cleaning up the object.
gene.data <- gene.data[gene.data$type=="gene"]
names(gene.data) <- gene.data$gene_id
is.gene.related <- grep("gene_", colnames(mcols(gene.data)))
mcols(gene.data) <- mcols(gene.data)[,is.gene.related]
rowRanges(sce) <- gene.data[rownames(sce)]
## add meta info
meta <- read.csv("../input/meta/bulkrna_scrna_meta/sixhosp_singlecell_rna_seq_meta.csv")


coldata <- colData(sce) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "row") %>% 
  dplyr::left_join(meta) %>% 
  tibble::column_to_rownames(var = "row") %>%
  DataFrame()

colData(sce) <- coldata
## now sce is ready for downstream analysis
saveRDS(sce, file = "scrna_19sample_sce.RDS")

# 
# 
# ## meta
# ```{r include=FALSE}
# library(readxl)
# library(dplyr)
# library(SeuratDisk)
# ```
# 
# 
# ```{r}
# sixhosp_bulk_rna_seq_meta <- read_excel("../input/meta/bulkrna_scrna_meta/sixhosp_bulk_rna_seq_meta.xlsx") %>% data.frame() %>% dplyr::select(sample_id,age,sex,diagnosis) %>% dplyr::distinct()
# colnames(sixhosp_bulk_rna_seq_meta)[1]="sample"
# ```
# 
# 
# ## construct singlecellexperiment
# ```{r}
# sixhosp_21subjects_kallisto <- LoadH5Seurat("/data3/zhiyu/pipelines/scRNA/nfcore_scRNA_sixhosp_21subjects_kallisto/output/kallisto/mtx_conversions/combined_matrix.h5seurat")
# sce <- as.SingleCellExperiment(sixhosp_21subjects_kallisto)
# 
# col_dat <- data.frame(colData(sce))
# common_id=intersect(col_dat$sample,sixhosp_bulk_rna_seq_meta$sample)
# meta=merge(col_dat,sixhosp_bulk_rna_seq_meta,by="sample",all.x=TRUE)
# meta$sample=as.character(meta$sample)
# meta$diagnosis=as.character(meta$diagnosis)
# 
# colData(sce)=meta %>% DataFrame() # 把新的meta表融合到sce
# sce=sce[,!is.na(sce$diagnosis) ]
# 
# 
# saveRDS(sce, file = "scrna_19sample_sce.RDS")
# 
# ```
# 
# 
# 
