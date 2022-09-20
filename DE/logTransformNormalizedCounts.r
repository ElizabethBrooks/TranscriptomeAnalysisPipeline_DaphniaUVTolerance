# Load libraries
library(ggplot2)
library(dplyr)

#Set working directory
#workingDir = args[1];
workingDir="/Users/bamflappy/PfrenderLab/dMelUV/DEA_PA42_v4.1/glmQLFAnalysis/"
setwd(workingDir); 

#Turn off scientific notation
options(scipen = 999)

#Import gene counts
inputTable <- read.csv(file="glmQLF_normalizedCounts.csv")

#Trim the data table
inputCounts <- head(inputTable, - 5)

#Set output file name
outFile <- "glmQLF_normalizedCounts_logTransformed.csv"

#Perform log transformations of log2(x+1)
outputCounts <- inputCounts
for(i in 2:ncol(inputCounts)) {
  outputCounts[,i] <- log2(inputCounts[,i]+1)
}

#Write log transformed counts to a file
write.csv(outputCounts, file=outFile)

#Test plots
ggplot(inputCounts, aes(x=seq_along(Y05_VIS_Pool1), y=Y05_VIS_Pool1)) + 
  geom_point()

ggplot(outputCounts, aes(x=seq_along(Y05_VIS_Pool1), y=Y05_VIS_Pool1)) + 
  geom_point()
