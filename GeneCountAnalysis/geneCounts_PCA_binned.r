#!/usr/bin/env Rscript
#Usage: Rscript geneCounts_PCA_binned.r mergedCounts_file_transposed.csv bin
#Usage Ex: Rscript geneCounts_PCA_binned.r mergedCounts_legacy_transposed.csv treatment
#R script to generate binned PCA of gene count matrices
#{ggfortify} lets {ggplot2} know how to interpret PCA objects
library(ggfortify)
#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)
#Test if there are three input arguments
if (length(args)!=2) {
  stop("One file name and one column name must be supplied.n", call.=FALSE)
}
#Retrieve gene count tables with gene IDs
gCount0 = read.csv(args[1], sep=",", row.names=1)
#Plot principal componants of PCA performed with prcomp
autoplot(prcomp(gCount0[, names(gCount0)!=args[2]]), data=gCount0, colour=args[2])