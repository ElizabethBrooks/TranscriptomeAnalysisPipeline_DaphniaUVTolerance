# Load libraries
library(ggplot2)
library("dplyr")

#Set working directory
#workingDir = args[1];
workingDir="/Users/bamflappy/PfrenderLab/dMelUV/DEA_PA42_v4.1/glmQLFAnalysis/"
setwd(workingDir); 

#Turn off scientific notation
options(scipen = 999)

#Import gene counts
inputCounts <- read.csv(file="glmQLF_normalizedCounts.csv")

#Set output file name
outFile <- "glmQLF_normalizedCounts_logTransformed.csv"

#Perform log transformations of log2(x+1)
outputCounts <- inputCounts
for(i in 2:ncol(inputCounts)) {
  colName <- colnames(inputCounts)[i]
  outputCounts[,i] <- log2(inputCounts[,i])
}

#Write log transformed counts to a file
write.csv(outputCounts, file=outFile)
