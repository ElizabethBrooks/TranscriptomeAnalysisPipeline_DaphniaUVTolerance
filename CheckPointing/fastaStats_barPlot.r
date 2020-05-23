#!/usr/bin/env Rscript
#Usage: Rscript fastaStats_barPlot.r statsSummaryFile
#Usage Ex: Rscript fastaStats_barPlot.r /afs/crc.nd.edu/group/pfrenderlab/mendel/ebrooks/rnaseq/trimmed_run1_assemblyTrinity_mergedFasta/trimmed_run1_assemblyTrinity_mergedFasta_summary.csv
#R script to generate bar plots of file stats

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
#Test if there is one input argument
if (numArgs!=1) {
  stop("A file name must be supplied.n", call.=FALSE)
}

#Retrieve file stats
aStats <- read.csv(file=args[1])
#Create data frame of combined file stats
counts <- data.frame(aStats$file, aStats$sequences, aStats$lines, aStats$bytes)

#Create matrix for multiple plots
par(mfrow=c(3,1))

#Set the plot titles
plotTitle <- "File Statistics"

#Generate grouped bar plot
plotSeqs <- ggplot(counts, aes(x=aStats.file, y=aStats.sequences)) + 
  geom_bar(stat="identity", fill="steelblue") +
  theme_minimal() +
  ggtitle(plotTitle) +
  xlab("File Name") +
  ylab("Number of Sequences")
#Save sequence stats plot as a jpg
outFile <- paste(normalizePath(dirname(args[1])), "plotSequences.jpg", sep="/")
print(outFile)
ggsave(outFile)

#Generate second grouped bar plot
plotLines <- ggplot(counts, aes(x=aStats.file, y=aStats.lines) + 
  geom_bar(stat="identity", fill="steelblue") +
  theme_minimal() + 
  ggtitle(plotTitle) +
  xlab("File Name") +
  ylab("Number of Lines")
#Save line stats plot as a jpg
outFile <- paste(normalizePath(dirname(args[1])), "plotLines.jpg", sep="/")
ggsave(outFile)

#Generate second grouped bar plot
plotBytes <- ggplot(counts, aes(x=aStats.file, y=aStats.bytes) + 
  geom_bar(stat="identity", fill="steelblue")+
  theme_minimal() + 
  ggtitle(plotTitle) +
  xlab("File Name") +
  ylab("Number of Bytes")
#Save bytes stats plot as a jpg
outFile <- paste(normalizePath(dirname(args[1])), "plotBytes.jpg", sep="/")
ggsave(outFile)
