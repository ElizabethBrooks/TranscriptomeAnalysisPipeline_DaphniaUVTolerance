#!/usr/bin/env Rscript
#Usage: Rscript alignmentSummary_barPlot_binnedCompareSoftware.r alignmentSummaryFiles
#Usage Ex: Rscript alignmentSummary_barPlot_binnedCompareSoftware.r alignmentSummarized_hisat2_run1_formatted.csv alignmentSummarized_tophat2_run2_formatted.csv
#R script to generate grouped and colored bar plots

#The easiest way to get ggplot2 is to install the whole tidyverse
#install.packages("tidyverse")
#Import library
library(ggplot2)
#Retrieve input file name of gene counts
args=commandArgs(trailingOnly=TRUE)
numArgs=length(args)
#Test if there is one input argument
if (length(args)!=2) {
  stop("Two file names must be supplied.n", call.=FALSE)
}
#Retrieve alignment stats
aStats <- do.call(rbind,lapply(args,read.csv))
#Re-set row names to match samples
subsetLength <- length(rownames(aStats))/2
subsetNames <- rownames(aStats)[1:subsetLength]
fullsetNames <- c(subsetNames,subsetNames)
fullsetNames <- as.numeric(fullsetNames)
#Create data frame of compined alignment stats
counts <- data.frame(fullsetNames, aStats$overall, aStats$concordant, aStats$software)
#Create matrix for multiple plots
par(mfrow=c(2,1))
#Generate grouped and colored bar plot
plotTitle <- basename(args[1])
plotTitle <- str_remove(plotTitle, "alignmentSummarized_")
plotTitle <- str_remove(plotTitle, "_formatted.csv")
plotOverall <- ggplot(counts, aes(factor(fullsetNames), aStats.overall, fill=aStats.software)) + 
  geom_bar(stat="identity", position="stack") +
  ggtitle(plotTitle) +
  xlab("Sample Number") +
  ylab("Overall Percent") +
  scale_fill_brewer(palette="Set1")
plotOverall <- plotOverall + guides(fill=guide_legend(title="Software"))
#Save overall percentages plot as a jpg
outFile <- paste(normalizePath(dirname(args[1])), "plotOverallPercentages.jpg", sep="/")
ggsave(outFile)
#Generate second grouped and colored bar plot
plotTitle <- basename(args[2])
plotTitle <- str_remove(plotTitle, "alignmentSummarized_")
plotTitle <- str_remove(plotTitle, "_formatted.csv")
plotConc <- ggplot(counts, aes(factor(fullsetNames), aStats.concordant, fill=aStats.software)) + 
  geom_bar(stat="identity", position="stack") + 
  ggtitle(plotTitle) +
  xlab("Sample Number") +
  ylab("Concordant Percent") +
  scale_fill_brewer(palette="Set1")
plotConc <- plotConc + guides(fill=guide_legend(title="Software"))
#Save concordant percentages plot as a jpg
outFile <- paste(normalizePath(dirname(args[1])), "plotConcordantPercentages.jpg", sep="/")
ggsave(outFile)
