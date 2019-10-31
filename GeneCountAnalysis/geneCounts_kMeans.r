#!/usr/bin/env Rscript
#Usage: Rscript geneCounts_kMeans.r mergedCounts_file_transposed.csv k SCALE
#Usage Ex: Rscript geneCounts_kMeans.r mergedCounts_legacy_transposed.csv 3 FALSE
#R script to generate kMeans clustering of gene count matrices
#{ggfortify} lets {ggplot2} know how to interpret kMeans objects
library(ggfortify)
#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)
#Test if there are three input arguments
if (length(args)!=3) {
  stop("One file name must be supplied.n", call.=FALSE)
}
#Retrieve gene count tables with gene IDs from input file
gCount0 = read.csv(args[1], sep=",", row.names=1)
#Plot principal componants of kMeans clustering performed with prcomp
autoplot(kmeans(gCount0, args[2], scale.=args[3]), data=gCount0, label=TRUE, label.size=3)