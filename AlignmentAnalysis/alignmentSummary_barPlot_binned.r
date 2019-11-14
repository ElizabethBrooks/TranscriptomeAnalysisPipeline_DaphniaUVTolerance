#!/usr/bin/env Rscript
#Usage: Rscript alignmentSummary_barPlot_binned.r countsFile.csv colName
#Usage Ex: Rscript alignmentSummary_barPlot_binned.r ../AlignmentStats_Analysis/alignmentSummarized_merged_differences.csv overall
#R script to generate grouped and colored bar plots
#Import library
library(ggplot2)
#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)
#Test if there are two input arguments
if (length(args)!=2) {
  stop("One file name and one column name must be supplied.n", call.=FALSE)
}
#Retrieve gene count tables with gene IDs
aStats = read.csv(args[1] sep=",", row.names=1)
counts <- table(aStats$sample, aStats$args[2], aStats$method)
#Generate grouped and colored bar plot
ggplot(counts, aes(factor(sample), args[2], fill=method)) + 
  geom_bar(stat="identity", position="dodge") + 
  scale_fill_brewer(palette="Set1")