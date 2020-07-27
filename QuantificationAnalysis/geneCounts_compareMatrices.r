#!/usr/bin/env Rscript
#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)
#Test if there is one input argument
if (length(args)!=2) {
  stop("One file name must be supplied.n", call.=FALSE)
}
#R script to compare two gene count matrices
#Retrieve gene count tables with gene IDs
gCount1 = read.csv(args[1], sep=",", row.names=1)
gCount2 = read.csv(args[2], sep=",", row.names=1)
#Convert gene count tables to matrices
gCount1M <- as.matrix(gCount1)
gCount2M <- as.matrix(gCount2)
#Determine matrix equality
eq <- gCount1M==gCount2M
#Get the percentage of non equal values
round(sum(!eq)/length(eq)*100, 2)