#!/usr/bin/env Rscript
#Usage: Rscript geneCounts_kMeans_binned.r mergedCounts_file_transposed.csv bin k
#Usage Ex: Rscript geneCounts_kMeans_binned.r mergedCounts_legacy_transposed.csv treatment 3
#R script to generate binned kMeans clustering of gene count matrices
#{ggfortify} lets {ggplot2} know how to interpret kMeans objects
library(ggfortify)
#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)
#Test if there are four input arguments
if (length(args)!=3) {
  stop("One file name, one column name, and the value for k must be supplied.n", call.=FALSE)
}
#Retrieve gene count tables with gene IDs from input file
gCount0 = read.csv(args[1], sep=",", row.names=1)
#Plot principal componants of kMeans clustering performed with prcomp
autoplot(kmeans(gCount0[, names(gCount0)!=args[2]], args[3]), data=gCount0, label=TRUE, label.size=3)