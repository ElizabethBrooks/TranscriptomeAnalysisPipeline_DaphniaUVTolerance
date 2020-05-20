#!/usr/bin/env Rscript
#Usage: Rscript alignmentSummary_barPlotMedians_binnedCompareGenotype alignmentSummaryFiles
#Usage Ex: Rscript alignmentSummary_barPlotMedians_binnedCompareGenotype alignmentSummarized_hisat2_run1_formatted.csv alignmentSummarized_tophat2_run2_formatted.csv
#R script to generate grouped and colored bar plots of genotype medians

#Installations need to be performed once
#The easiest way to get ggplot2 is to install the whole tidyverse
#install.packages("tidyverse")
#install.packages("stringr")
#Import librarys
library(ggplot2)
library(stringr)
#Retrieve input file name of gene counts
args=commandArgs(trailingOnly=TRUE)
numArgs=length(args)
#Retrieve alignment stats
aStats <- do.call(rbind,lapply(args[2:numArgs],read.csv))
#Create data frame of compined alignment stats
counts <- data.frame(aStats$genotype, aStats$overallMedian, aStats$concordantMedian, aStats$overallSd, aStats$concordantSd, aStats$targetGenotype)
#Create matrix for multiple plots
par(mfrow=c(2,1))
#Set the plot titles
plotTitle <- args[1]
plotTitle <- gsub("_", " ", plotTitle)
#Generate grouped and colored bar plot
plotOverall <- ggplot(counts, aes(x=aStats.genotype, y=aStats.overallMedian, fill=aStats.targetGenotype)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(x=aStats.genotype, ymin=aStats.overallMedian-aStats.overallSd, ymax=aStats.overallMedian+aStats.overallSd), width=.2, position=position_dodge(.9)) +
  ggtitle(plotTitle) +
  xlab("Sample Genotype") +
  ylab("Median Overall Percent") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Set1")
plotOverall <- plotOverall + guides(fill=guide_legend(title="Target Genotype"))
#Save overall percentages plot as a jpg
outFile <- paste(normalizePath(dirname(args[2])), "plotMedianOverallPercentages.jpg", sep="/")
ggsave(outFile)
#Generate second grouped and colored bar plot
plotConc <- ggplot(counts, aes(x=aStats.genotype, y=aStats.concordantMedian, fill=aStats.targetGenotype)) +
  geom_bar(stat="identity", position="dodge") + 
  geom_errorbar(aes(x=aStats.genotype, ymin=aStats.concordantMedian-aStats.concordantSd, ymax=aStats.concordantMedian+aStats.concordantSd), width=.2, position=position_dodge(.9)) +
  ggtitle(plotTitle) +
  xlab("Sample Genotype") +
  ylab("Median Concordant Percent") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Set1")
plotConc <- plotConc + guides(fill=guide_legend(title="Target Genotype"))
#Save concordant percentages plot as a jpg
outFile <- paste(normalizePath(dirname(args[2])), "plotMedianConcordantPercentages.jpg", sep="/")
ggsave(outFile)
