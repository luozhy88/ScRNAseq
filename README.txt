Barcode是每个凝胶微珠的身份证号码，区分细胞
UMI是每个DNA标签分子的身份证号码，区分分子
Poly(dT)VN作用是与mRNA的Ploy(A)尾巴相结合，作为逆转录的引物，逆转录出cDNA
https://www.youtube.com/watch?v=dbE1UlpxzHQ

# 质控
https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/advanced/barcode-rank-plot#understanding-data


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

# BPCells
https://bnprks.github.io/BPCells/articles/pbmc3k.html
