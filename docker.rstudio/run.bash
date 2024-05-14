
#https://nbisweden.github.io/workshop-scRNAseq/other/containers.html
docker run --platform=linux/amd64 --rm -p 8789:8787 -e PASSWORD=scrnaseq -v ${PWD}:/home/rstudio/workdir ghcr.io/nbisweden/workshop-scrnaseq:2024-bioconductor-r4.3.0
