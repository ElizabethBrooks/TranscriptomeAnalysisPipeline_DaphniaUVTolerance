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

#Perform log transformations of log2(x+1)
for(i in 2:24) {
  var <- "logCol"
  colName <- paste(var,i, sep = "", collapse = "")
  inputCounts[[colName]] <- log2(inputCounts[,i])
}
