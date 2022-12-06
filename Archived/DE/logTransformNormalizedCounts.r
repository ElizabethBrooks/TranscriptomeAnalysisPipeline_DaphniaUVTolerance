#!/usr/bin/env Rscript

# R script to log transform normalized gene counts

# Load libraries
library(ggplot2)
library(dplyr)

#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)

#Set working directory
#workingDir = args[1];
workingDir = "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Tolerance"
setwd(workingDir); 

#Turn off scientific notation
options(scipen = 999)

#Import gene counts
inputCounts <- read.csv(file="glmQLF_normalizedCounts.csv", row.names="gene")

#Perform log transformations of log2(x+1)
outputCounts <- log2(inputCounts+1)

#Write log transformed counts to a file
write.table(outputCounts, file="normalizedCounts_logTransformed.csv", sep=",", row.names=TRUE, quote=FALSE)

#Test plots
jpeg("glmQLF_normalizedCounts_scatterPlot.jpg")
ggplot(inputCounts, aes(x=seq_along(Y05_VIS_Pool1), y=Y05_VIS_Pool1)) + 
  geom_point()
dev.off()

jpeg("glmQLF_normalizedCounts_logTransformed_scatterPlot.jpg")
ggplot(outputCounts, aes(x=seq_along(Y05_VIS_Pool1), y=Y05_VIS_Pool1)) + 
  geom_point()
dev.off()
