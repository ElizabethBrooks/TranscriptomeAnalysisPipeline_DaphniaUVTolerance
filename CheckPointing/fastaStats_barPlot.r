#!/usr/bin/env Rscript
#Usage: Rscript fastaStats_barPlot.r statsSummaryFile
#Usage Ex: Rscript fastaStats_barPlot.r /afs/crc.nd.edu/group/pfrenderlab/mendel/ebrooks/rnaseq/trimmed_run1_assemblyTrinity_mergedFasta/trimmed_run1_assemblyTrinity_mergedFasta_summary.csv
#R script to generate bar plots of file stats

#Installations need to be performed once
#The easiest way to get ggplot2 is to install the whole tidyverse
#install.packages("tidyverse")
#install.packages("ggpubr")

#Import librarys
library(ggplot2)
library(ggpubr)

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
  geom_text(aes(label=aStats.sequences), vjust=1.6, color="white", size=3.5) +
  xlab("File") +
  ylab("Sequences")

#Generate second grouped bar plot
plotLines <- ggplot(counts, aes(x=aStats.file, y=aStats.lines)) + 
  geom_bar(stat="identity", fill="steelblue") +
  theme_minimal() + 
  geom_text(aes(label=aStats.lines), vjust=1.6, color="white", size=3.5) +
  xlab("File") +
  ylab("Lines")

#Generate second grouped bar plot
plotMBytes <- ggplot(counts, aes(x=aStats.file, y=aStats.bytes)) + 
  geom_bar(stat="identity", fill="steelblue")+
  theme_minimal() +
  geom_text(aes(label=aStats.bytes), vjust=1.6, color="white", size=3.5) +
  xlab("File") +
  ylab("MB")

#Arrange stats plots on grid
finalPlot <- ggarrange(plotSeqs, plotLines, plotMBytes, nrow=3)
#Add plot title
finalPlot <- annotate_figure(finalPlot,
  top=text_grob(plotTitle, color="black", face="bold", size=14))

#Save file stats plots as a jpg
outFile <- paste(normalizePath(dirname(args[1])), "fileStats.jpg", sep="/")
ggexport(finalPlot, filename=outFile)
