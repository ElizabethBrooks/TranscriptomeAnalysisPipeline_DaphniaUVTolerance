#Set working directory
#workingDir = args[1];
workingDir="~/PfrenderLab/dMelUV/WGCNA_PA42_v4.1"
setwd(workingDir); 

# Load libraries
library(ggplot2)

# Load the expression and trait data saved in the first part
lnames1 = load(file = "PA42_v4.1_dataInputTol.RData");

# Load network data saved in the second part.
lnames2 = load(file = "PA42_v4.1_networkConstructionTol_auto_threshold8_signedNowick.RData");
numRow = length(moduleColors)

#Import ka ks results
inputKaks <- read.csv(file="~/PfrenderLab/PA42_v4.1/PA42_v4.1_Olympics_kaksResults.csv")
cleanKaks <- na.omit(inputKaks)
filtKaks <- cleanKaks[which(cleanKaks$dNdS <= 5),]
#numRow = nrow(filtKaks)

