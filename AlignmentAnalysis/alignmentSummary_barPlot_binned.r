#!/usr/bin/env Rscript
#Usage: Rscript alignmentSummary_barPlot_binned.r percentsFile.csv
#Usage Ex: Rscript alignmentSummary_barPlot_binned.r alignmentSummarized_legacyTophat2Hisat2_merged.csv
#R script to generate grouped and colored bar plots
#Import library
library(ggplot2)
#Retrieve input file name of gene counts
args=commandArgs(trailingOnly=TRUE)
#Test if there are two input arguments
if (length(args)!=1) {
  stop("One file name must be supplied.n", call.=FALSE)
}
#Retrieve alignment stats
aStats = read.csv(args[1], sep=",")
counts <- data.frame(aStats$sampleNum, aStats$overall, aStats$concordant, aStats$method)
#Create matrix for multiple plots
par(mfrow=c(2,1))
#Generate grouped and colored bar plot
ggplot(counts, aes(factor(aStats.sampleNum), aStats.overall, fill=aStats.method)) + 
  geom_bar(stat="identity", position="dodge") + 
  scale_fill_brewer(palette="Set1")
  #Generate grouped and colored bar plot
ggplot(counts, aes(factor(aStats.sampleNum), aStats.concordant, fill=aStats.method)) + 
  geom_bar(stat="identity", position="dodge") + 
  scale_fill_brewer(palette="Set1")