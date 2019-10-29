#!/usr/bin/env Rscript
#R script to generate scree plot for binned PCA of gene count matrices
#Library {tidyverse} used to generate scree plots
library(tidyverse)
#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)
#Test if there is one input argument
if (length(args)!=2) {
  stop("One file name must be supplied.n", call.=FALSE)
}
#Retrieve gene count tables with gene IDs from input file
gCount0 = read.csv(args[1], sep=",", row.names=1)
#Generate principal componants of PCA performed with prcomp
gene.pca <- prcomp(gCount0[ , names(gCount0) != args[2]])
# Extract eigenvalues/variances
get_eig(gene.pca)
#Default plot
#fviz_eig(gene.pca, addlabels = TRUE, ylim = c(0, 85))
#Scree plot - Eigenvalues
fviz_eig(gene.pca, choice = "eigenvalue", addlabels=TRUE)
#Use only bar or line plot: geom = "bar" or geom = "line"
#fviz_eig(gene.pca, geom="line")