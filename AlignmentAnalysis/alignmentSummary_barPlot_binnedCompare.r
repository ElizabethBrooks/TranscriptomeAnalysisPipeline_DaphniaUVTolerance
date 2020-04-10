#!/usr/bin/env Rscript
#Usage: Rscript alignmentSummary_barPlot_binned.r alignmentSummaryFiles inputFilesPath
#Usage Ex: Rscript alignmentSummary_barPlot_binned.r alignmentSummarized_hisat2_run1_formatted.csv alignmentSummarized_tophat2_run2_formatted.csv inputFilesPath
#R script to generate grouped and colored bar plots

#The easiest way to get ggplot2 is to install the whole tidyverse
#install.packages("tidyverse")
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
#Re-set row names to match samples
subsetLength <- length(rownames(aStats))/2
subsetNames <- rownames(aStats)[1:subsetLength]
fullsetNames <- c(subsetNames,subsetNames)
fullsetNames <- as.numeric(fullsetNames)
#Create data frame of compined alignment stats
counts <- data.frame(fullsetNames, aStats$overall, aStats$concordant, aStats$method)
#Create matrix for multiple plots
par(mfrow=c(2,1))
#Generate grouped and colored bar plot
plotOverall <- ggplot(counts, aes(factor(fullsetNames), aStats.overall, fill=aStats.method)) + 
  geom_bar(stat="identity", position="stack") +
  xlab("Sample Number") +
  ylab("Overall Percent") +
  scale_fill_brewer(palette="Set1")
plotOverall <- plotOverall + guides(fill=guide_legend(title="Software"))
#Save overall percentages plot as a jpg
outFile <- paste(normalizePath(dirname(args[1])), "plotOverallPercentages.jpg", sep="/")
jpeg(outFile)
grid.newpage()
grid.draw(plotOverall)
dev.off()
#Generate second grouped and colored bar plot
plotConc <- ggplot(counts, aes(factor(fullsetNames), aStats.concordant, fill=aStats.method)) + 
  geom_bar(stat="identity", position="stack",by=2) + 
  xlab("Sample Number") +
  ylab("Concordant Percent") +
  scale_fill_brewer(palette="Set1")
plotConc <- plotConc + guides(fill=guide_legend(title="Software"))
#Save concordant percentages plot as a jpg
outFile <- paste(normalizePath(dirname(args[1])), "plotOverallPercentages.jpg", sep="/")
jpeg(outFile)
grid.newpage()
grid.draw(plotConc)
dev.off()
