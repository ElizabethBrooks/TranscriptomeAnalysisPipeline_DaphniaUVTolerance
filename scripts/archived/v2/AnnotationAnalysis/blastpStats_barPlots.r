#!/usr/bin/env Rscript
#Usage: Rscript blastpStats_barPlots.r title blastpSummaryFile
#Usage Ex: Rscript blastpStats_barPlots.r trimmed_run1 trimmed_run1_blastp_summary.txt
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
statsMerged <- data.frame(aStats$query, aStats$queryHits, aStats$dbHits, aStats$bestHits)

#Create matrix for multiple plots
par(mfrow=c(3,1))

#Set the plot titles and output file
plotTitle1 <- "Transcriptome"
plotTitle2 <- "Blastp Hits to"
plotTitle3 <- aStats$db[1]
plotTitle4 <- paste(plotTitle1, plotTitle2, sep=" ")
plotTitle <- paste(plotTitle4, plotTitle3, sep=" ")

#Generate bar plot of best hits
plotBest <- ggplot(statsMerged, aes(x=aStats.query, y=aStats.bestHits)) + 
  geom_bar(stat="identity", fill="steelblue") +
  theme_minimal() +
  geom_text(aes(label=aStats.bestHits), vjust=1.6, color="white", size=3.5) +
  ylab("Best Hits") +
  theme(axis.title.x = element_blank())

#Generate bar plot of query hits
plotQuery <- ggplot(statsMerged, aes(x=aStats.query, y=aStats.queryHits)) + 
  geom_bar(stat="identity", fill="steelblue") +
  theme_minimal() + 
  geom_text(aes(label=aStats.queryHits), vjust=1.6, color="white", size=3.5) +
  ylab("Query Hits") +
  theme(axis.title.x = element_blank())

#Generate bar plot of DB hits
plotDB <- ggplot(statsMerged, aes(x=aStats.query, y=aStats.dbHits)) + 
  geom_bar(stat="identity", fill="steelblue") +
  theme_minimal() + 
  geom_text(aes(label=aStats.dbHits), vjust=1.6, color="white", size=3.5) +
  ylab("DB Hits") +
  theme(axis.title.x = element_blank())


#Arrange stats plots on grid
finalPlot <- ggarrange(plotBest, plotQuery, plotDB, nrow=3)
#Add plot title and x-axis label
finalPlot <- annotate_figure(finalPlot,
  top=text_grob(plotTitle, color="black", face="bold", size=14),
  bottom=text_grob("Transcriptome Genotype", color="black", size=12))

#Set output file names
plotOut1 <- args[1]
plotOut2 <- "blastpStats.jpg"
plotFile1 <- paste(plotOut1, plotTitle3, sep="_")
plotFile <- paste(plotFile1, plotOut2, sep="_")
#Set output path to the input file path
outFile <- paste(normalizePath(dirname(args[2])), plotFile, sep="/")

#Save file stats plots as a jpg
ggexport(finalPlot, filename=outFile)
