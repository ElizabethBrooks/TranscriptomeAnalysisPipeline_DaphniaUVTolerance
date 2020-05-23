#!/usr/bin/env Rscript
#Usage: Rscript fastaStats_barPlot.r statsSummaryFile
#Usage Ex: Rscript fastaStats_barPlot.r mergedFasta_summary.csv
#R script to generate grouped and colored bar plots

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
  stop("A file name must be supplied... exiting.n", call.=FALSE)
}

#Retrieve file stats
aStats <- read.csv(file=args[1])
#Create data frame of combined file stats
counts <- data.frame(aStats$file, aStats$sequences, aStats$lines, aStats$bytes)

#Create matrix for multiple plots
par(mfrow=c(2,1))

#Set the plot titles
plotTitle <- "File Statistics"

#Generate grouped bar plot
plotOverall <- ggplot(counts, aes(x=aStats.file, y=aStats.sequences)) + 
  geom_bar(stat="identity", position="dodge") +
  ggtitle(plotTitle) +
  xlab("File Name") +
  ylab("Number of Sequences") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Set1")
#Save sequence stats plot as a jpg
outFile <- paste(normalizePath(dirname(args[1])), "plotSequences.jpg", sep="/")
ggsave(outFile)

#Generate second grouped bar plot
plotConc <- ggplot(counts, aes(x=aStats.file, y=aStats.lines) + 
  geom_bar(stat="identity", position="dodge") + 
  ggtitle(plotTitle) +
  xlab("File Name") +
  ylab("Number of Lines") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Set1")
#Save line stats plot as a jpg
outFile <- paste(normalizePath(dirname(args[1])), "plotLines.jpg", sep="/")
ggsave(outFile)

#Generate second grouped bar plot
plotConc <- ggplot(counts, aes(x=aStats.file, y=aStats.bytes) + 
  geom_bar(stat="identity", position="dodge") + 
  ggtitle(plotTitle) +
  xlab("File Name") +
  ylab("Number of Bytes") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Set1")
#Save bytes stats plot as a jpg
outFile <- paste(normalizePath(dirname(args[1])), "plotBytes.jpg", sep="/")
ggsave(outFile)
