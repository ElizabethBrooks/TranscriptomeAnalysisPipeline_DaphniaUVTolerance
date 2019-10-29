#!/usr/bin/env Rscript
#R script to generate kMeans clustering of gene count matrices
#{ggfortify} lets {ggplot2} know how to interpret kMeans objects
library(ggfortify)
#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)
#Test if there is one input argument
if (length(args)!=1) {
  stop("One file name must be supplied.n", call.=FALSE)
}
#Retrieve gene count tables with gene IDs from input file
gCount0 = read.csv(args[1], sep=",", row.names=1)
#Plot principal componants of kMeans clustering performed with prcomp
autoplot(kmeans(gCount0, 3), data = gCount0, label = TRUE, label.size = 3)