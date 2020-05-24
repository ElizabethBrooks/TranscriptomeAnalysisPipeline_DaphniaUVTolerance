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
if (numArgs!=2) {
  stop("A file name and plot title must be supplied.n", call.=FALSE)
}

#Retrieve file stats
aStats <- read.csv(file=args[1])
#Create data frame of combined file stats
counts <- data.frame(aStats$file, aStats$sequences, aStats$lines, aStats$MB)

#Create matrix for multiple plots
par(mfrow=c(3,1))

#Set the plot titles
plotTitle1 <- args[2]
plotTitle2 <- "File Stats"
plotTitle <- paste(plotTitle1, plotTitle2, sep=" ")

#Generate grouped bar plot
plotSeqs <- ggplot(counts, aes(x=aStats.file, y=aStats.sequences)) + 
  geom_bar(stat="identity", fill="steelblue") +
  theme_minimal() +
  geom_text(aes(label=aStats.sequences), vjust=1.6, color="white", size=3.5) +
  ylab("Sequences") +
  axis.title.x = element_blank()

#Generate second grouped bar plot
plotLines <- ggplot(counts, aes(x=aStats.file, y=aStats.lines)) + 
  geom_bar(stat="identity", fill="steelblue") +
  theme_minimal() + 
  geom_text(aes(label=aStats.lines), vjust=1.6, color="white", size=3.5) +
  ylab("Lines") +
  axis.title.x = element_blank()

#Generate second grouped bar plot
plotMB <- ggplot(counts, aes(x=aStats.file, y=aStats.MB)) + 
  geom_bar(stat="identity", fill="steelblue")+
  theme_minimal() +
  geom_text(aes(label=aStats.MB), vjust=1.6, color="white", size=3.5) +
  ylab("MB") +
  axis.title.x = element_blank()

#Arrange stats plots on grid
finalPlot <- ggarrange(plotSeqs, plotLines, plotMB, nrow=3)
#Add plot title
finalPlot <- annotate_figure(finalPlot,
  top=text_grob(plotTitle, color="black", face="bold", size=14),
  bottom=text_grob("File Number", color="black", face="bold", size=12))

#Save file stats plots as a jpg
outFile <- paste(normalizePath(dirname(args[1])), "fileStats.jpg", sep="/")
ggexport(finalPlot, filename=outFile)
