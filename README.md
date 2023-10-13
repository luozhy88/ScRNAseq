# ScRNAseq
https://bioconductor.org/books/3.13/OSCA.intro/the-singlecellexperiment-class.html#overview-2
## Only keeping the first two assays
assays(sce) <- assays(sce)[1:2]
