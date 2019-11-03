#!/usr/bin/env Rscript
#Usage: Rscript geneCounts_PCA.r mergedCounts_file_transposed.csv
#Usage Ex: Rscript geneCounts_PCA.r mergedCounts_legacy_transposed.csv
#R script to generate PCA of gene count matrices
#{ggfortify} lets {ggplot2} know how to interpret PCA objects
library(ggfortify)
#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)
#Test if there are two input arguments
if (length(args)!=1) {
  stop("One file name must be supplied.n", call.=FALSE)
}
#Retrieve gene count tables with gene IDs from input file
gCount0 = read.csv(args[1], sep=",", row.names=1)
#Plot principal componants of PCA performed with prcomp
autoplot(prcomp(gCount0))