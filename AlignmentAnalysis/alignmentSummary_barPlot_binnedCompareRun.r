#!/usr/bin/env Rscript
#Usage: Rscript alignmentSummary_barPlot_binnedCompareRun.r alignmentSummaryFiles
#Usage Ex: Rscript alignmentSummary_barPlot_binnedCompareRun.r alignmentSummarized_hisat2_run1_formatted.csv alignmentSummarized_tophat2_run2_formatted.csv
#R script to generate grouped and colored bar plots

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
#Re-set row names to match samples
subsetLength <- length(rownames(aStats))/2
subsetNames <- rownames(aStats)[1:subsetLength]
fullsetNames <- c(subsetNames,subsetNames)
fullsetNames <- as.numeric(fullsetNames)
#Create data frame of compined alignment stats
counts <- data.frame(fullsetNames, aStats$overall, aStats$concordant, aStats$run)
#Calculate row median values for each genotype
curCols <- 2:3
for (i in nrow(counts)) {
  #Use a sliding window of 6 for my data set
  max <- i+5
  curRows <- i:max
  curRows
  counts[curRows,curCols]
  #Calulate row medians for each genotype column set
  #curMedians <- rowMedians(counts[curCols], rows=curRows, na.rm=FALSE, dim.=dim(counts))
  #curMedians
  #Add latest medians to final matrix
  #finalMedians <- Merge(finalMedians, curMedians, by=NULL)
  #finalMedians
  #Move the sliding window
  i <- i+5
}

