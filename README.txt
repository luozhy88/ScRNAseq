# ScRNAseq
https://bioconductor.org/books/3.13/OSCA.intro/the-singlecellexperiment-class.html#overview-2
## Only keeping the first two assays
assays(sce) <- assays(sce)[1:2]

# scFlow
https://combiz.github.io/scFlow/articles/scFlow.html
devtools::install_github("neurogenomics/scFlowExamples")
BiocManager::install("uniftest")
BiocManager::install("DropletUtils")
install.packages("ids")
devtools::install_github("NathanSkene/One2One")
devtools::install_github(repo = "hhoeflin/hdf5r")
devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
devtools::install_github("vqf/nVennR")
devtools::install_github("combiz/scFlow")
devtools::install_github("combiz/scFlowData")#https://api.github.com/repos/combiz/scFlowData/tarball/HEAD
