#!/usr/bin/env Rscript
#Usage: Rscript geneCounts_PCA_eigens_binned.r mergedCounts_file_transposed.csv bin
#Usage Ex: Rscript geneCounts_PCA_eigens_binned.r mergedCounts_legacy_transposed.csv treatment
#R script to generate binned PCA with loading vectors of gene count matrices
#{ggfortify} lets {ggplot2} know how to interpret PCA objects
library(ggfortify)
#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)
#Test if there is one input argument
if (length(args)!=2) {
  stop("One file name and one column name must be supplied.n", call.=FALSE)
}
#Retrieve gene count tables with gene IDs from input file
gCount0 = read.csv(args[1], sep=",", row.names=1)
#Plot principal componants and loading vectors of PCA performed with prcomp
autoplot(prcomp(gCount0[, names(gCount0)!=args[2]]), data=gCount0, colour=args[2], loadings=TRUE, 
	loadings.colour='green', loadings.label=TRUE, loadings.label.size=3)