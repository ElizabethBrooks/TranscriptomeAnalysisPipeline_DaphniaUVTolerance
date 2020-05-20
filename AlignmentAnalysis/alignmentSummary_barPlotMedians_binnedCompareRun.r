#!/usr/bin/env Rscript
#Usage: Rscript alignmentSummary_barPlotMedians_binnedCompareRun.r alignmentSummaryFiles
#Usage Ex: Rscript alignmentSummary_barPlotMedians_binnedCompareRun.r alignmentSummarized_hisat2_run1_formatted.csv alignmentSummarized_tophat2_run2_formatted.csv
#R script to generate grouped and colored bar plots of genotype medians

#Installations need to be performed once
#The easiest way to get ggplot2 is to install the whole tidyverse
#install.packages("tidyverse")
#install.packages("stringr")
#Import librarys
library(ggplot2)
library(stringr)
library(matrixStats)
#Retrieve input file name of gene counts
args=commandArgs(trailingOnly=TRUE)
numArgs=length(args)
#Test if there is one input argument
if (length(args)!=2) {
  stop("Two file names must be supplied.n", call.=FALSE)
}
#Retrieve alignment stats
aStats <- do.call(rbind,lapply(args,read.csv))
#Create data frame of compined alignment stats
counts <- data.frame(aStats$genotype, aStats$overallMedian, aStats$concordantMedian, aStats$overallSd, aStats$concordantSd, aStats$run)
#Create matrix for multiple plots
par(mfrow=c(2,1))
#Set the plot titles
plotTitle1 <- basename(args[1])
plotTitle1 <- str_remove(plotTitle1, "alignmentSummarized_")
plotTitle1 <- str_remove(plotTitle1, "_medians.csv")
plotTitle2 <- basename(args[2])
plotTitle2 <- str_remove(plotTitle2, "alignmentSummarized_")
plotTitle2 <- str_remove(plotTitle2, "_medians.csv")
plotTitle <- paste(plotTitle1, plotTitle2, sep=" vs ")
#Generate grouped and colored bar plot
plotOverall <- ggplot(counts, aes(x=aStats.genotype, y=aStats.overallMedian, fill=aStats.run)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(x=aStats.genotype, ymin=aStats.overallMedian-aStats.overallSd, ymax=aStats.overallMedian+aStats.overallSd), width=.2, position=position_dodge(.9)) +
  ggtitle(plotTitle) +
  xlab("Genotype") +
  ylab("Overall Median Percent") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Set1")
plotOverall <- plotOverall + guides(fill=guide_legend(title="Run Number"))
#Save overall percentages plot as a jpg
outFile <- paste(normalizePath(dirname(args[1])), "plotOverallMedianPercentages.jpg", sep="/")
ggsave(outFile)
#Generate second grouped and colored bar plot
plotConc <- ggplot(counts, aes(x=aStats.genotype, y=aStats.concordantMedian, fill=aStats.run)) +
  geom_bar(stat="identity", position="dodge") + 
  geom_errorbar(aes(x=aStats.genotype, ymin=aStats.concordantMedian-aStats.concordantSd, ymax=aStats.concordantMedian+aStats.concordantSd), width=.2, position=position_dodge(.9)) +
  ggtitle(plotTitle) +
  xlab("Genotype") +
  ylab("Overall Median Percent") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Set1")
plotConc <- plotConc + guides(fill=guide_legend(title="Run Number"))
#Save concordant percentages plot as a jpg
outFile <- paste(normalizePath(dirname(args[1])), "plotConcordantMedianPercentages.jpg", sep="/")
ggsave(outFile)
