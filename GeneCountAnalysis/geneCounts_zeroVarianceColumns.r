#!/usr/bin/env Rscript
#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)
#Test if there is one input argument
if (length(args)!=1) {
  stop("One file name must be supplied.n", call.=FALSE)
}
#R script to identify columns of matrix with zero-variance
#Retrieve gene count table with gene IDs
gCount0 = read.csv(args[1], sep=",", row.names=1)
#Convert gene count table to matrix
gCount0M <- as.matrix(gCount0)
#Identify columns with zero-variance
which(apply(gCount0M, 2, var)==0)