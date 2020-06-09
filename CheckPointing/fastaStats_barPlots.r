#!/usr/bin/env Rscript
#Usage: Rscript fastaStats_barPlots.r title statsSummaryFile
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
aStats <- read.csv(file=args[2])
totalIndex <- length(rownames(aStats))
endFiles <- totalIndex-2

#Create data frame of combined file stats
statsMerged <- data.frame(aStats$genotype, aStats$sequences, aStats$lines, aStats$MB)
statsFiles <- statsMerged[1:endFiles,]

#Create matrix for multiple plots
par(mfrow=c(3,1))

#Set the plot titles and output file
plotTitle1 <- args[1]
plotTitle2 <- "File Stats"
plotTitle <- paste(plotTitle1, plotTitle2, sep=" ")

#Generate sequences bar plots
plotSeqsMerged <- ggplot(statsMerged, aes(x=aStats.genotype, y=aStats.sequences)) + 
  geom_bar(stat="identity", fill="steelblue") +
  theme_minimal() +
  geom_text(aes(label=aStats.sequences), vjust=1.6, color="white", size=3.5) +
  ylab("Sequences") +
  theme(axis.title.x = element_blank())
plotSeqsFiles <- ggplot(statsFiles, aes(x=aStats.genotype, y=aStats.sequences)) + 
  geom_bar(stat="identity", fill="steelblue") +
  theme_minimal() +
  geom_text(aes(label=aStats.sequences), vjust=1.6, color="white", size=3.5) +
  ylab("Sequences") +
  theme(axis.title.x = element_blank())

#Generate lines bar plots
plotLinesMerged <- ggplot(statsMerged, aes(x=aStats.genotype, y=aStats.lines)) + 
  geom_bar(stat="identity", fill="steelblue") +
  theme_minimal() + 
  geom_text(aes(label=aStats.lines), vjust=1.6, color="white", size=3.5) +
  ylab("Lines") +
  theme(axis.title.x = element_blank())
plotLinesFiles <- ggplot(statsFiles, aes(x=aStats.genotype, y=aStats.lines)) + 
  geom_bar(stat="identity", fill="steelblue") +
  theme_minimal() + 
  geom_text(aes(label=aStats.lines), vjust=1.6, color="white", size=3.5) +
  ylab("Lines") +
  theme(axis.title.x = element_blank())

#Generate MB bar plots
plotMBMerged <- ggplot(statsMerged, aes(x=aStats.genotype, y=aStats.MB)) + 
  geom_bar(stat="identity", fill="steelblue")+
  theme_minimal() +
  geom_text(aes(label=aStats.MB), vjust=1.6, color="white", size=3.5) +
  ylab("MB") +
  theme(axis.title.x = element_blank())
plotMBFiles <- ggplot(statsFiles, aes(x=aStats.genotype, y=aStats.MB)) + 
  geom_bar(stat="identity", fill="steelblue")+
  theme_minimal() +
  geom_text(aes(label=aStats.MB), vjust=1.6, color="white", size=3.5) +
  ylab("MB") +
  theme(axis.title.x = element_blank())

#Arrange stats plots on grid
finalMergedPlot <- ggarrange(plotSeqsMerged, plotLinesMerged, plotMBMerged, nrow=3)
finalFilesPlot <- ggarrange(plotSeqsFiles, plotLinesFiles, plotMBFiles, nrow=3)
#Add plot title and x-axis label
finalMergedPlot <- annotate_figure(finalMergedPlot,
  top=text_grob(plotTitle, color="black", face="bold", size=14),
  bottom=text_grob("File Number", color="black", size=12))
finalFilesPlot <- annotate_figure(finalFilesPlot,
  top=text_grob(plotTitle, color="black", face="bold", size=14),
  bottom=text_grob("File Number", color="black", size=12))

#Set output file names
plotMergedOut <- "mergedStats.jpg"
plotFilesOut <- "fileStats.jpg"
plotMergedFile <- paste(plotTitle1, plotMergedOut, sep="_")
plotFilesFile <- paste(plotTitle1, plotFilesOut, sep="_")
#Set output path to the input file path
outMergedFile <- paste(normalizePath(dirname(args[2])), plotMergedFile, sep="/")
outFilesFile <- paste(normalizePath(dirname(args[2])), plotFilesFile, sep="/")

#Save file stats plots as a jpg
ggexport(finalMergedPlot, filename=outMergedFile)
ggexport(finalFilesPlot, filename=outFilesFile)
