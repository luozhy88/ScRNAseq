if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")

library(Seurat)
library(SeuratDisk)

Convert("combined_matrix.h5ad", dest = "h5seurat", overwrite = TRUE)
sixhosp_21subjects_alevin <- LoadH5Seurat("combined_matrix.h5seurat")