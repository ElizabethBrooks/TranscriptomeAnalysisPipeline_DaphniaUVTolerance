#!/usr/bin/env Rscript
#Usage: Rscript alignmentSummary_barPlot_binned.r alignmentSummaryFiles
#Usage Ex: Rscript alignmentSummary_barPlot_binned.r alignmentSummarized_hisat2_run1_formatted.csv alignmentSummarized_tophat2_run2_formatted.csv
#R script to generate grouped and colored bar plots

#Import library
library(ggplot2)
#Retrieve input file name of gene counts
args=commandArgs(trailingOnly=TRUE)
numArgs=length(args)
#Test if there is one input argument
if (length(args)!=1) {
  stop("One file name must be supplied.n", call.=FALSE)
}
#Retrieve alignment stats
aStats <- do.call(rbind,lapply(args,read.csv))
counts <- data.frame(aStats$sample, aStats$overall, aStats$concordant, aStats$method)
#Create matrix for multiple plots
par(mfrow=c(2,1))
#Generate grouped and colored bar plot
ggplot(counts, aes(factor(aStats.sample), aStats.overall, fill=aStats.method)) + 
  geom_bar(stat="identity", position="dodge") + 
  scale_fill_brewer(palette="Set1")
#Generate second grouped and colored bar plot
ggplot(counts, aes(factor(aStats.sample), aStats.concordant, fill=aStats.method)) + 
  geom_bar(stat="identity", position="dodge") + 
  scale_fill_brewer(palette="Set1")
