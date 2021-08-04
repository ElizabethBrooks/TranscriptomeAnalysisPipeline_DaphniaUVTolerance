#Set working directory
#workingDir = args[1];
workingDir="~/PfrenderLab/WGCNA_PA42_v4.1/SelectionAnalysis"
setwd(workingDir); 

#Load libraries
library(ggplot2)

#Import ka ks results
inputKaks <- read.csv(file="~/PfrenderLab/PA42_v4.1/PA42_v4.1_Olympics_kaksResults.txt")
cleanKaks <- na.omit(inputKaks)
filtKaks <- cleanKaks[which(cleanKaks$dNdS < 10),]

#Generate a scatter plot of dN and dS values
ggplot(filtKaks, aes(x=dS, y=dN)) +
  geom_point(aes(size=dNdS))
