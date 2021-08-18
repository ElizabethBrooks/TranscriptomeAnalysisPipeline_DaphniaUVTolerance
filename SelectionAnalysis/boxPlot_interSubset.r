#Set working directory
#workingDir = args[1];
workingDir="~/PfrenderLab/WGCNA_PA42_v4.1"
setwd(workingDir); 

# Load libraries
library(ggplot2)

# Load the expression and trait data saved in the first part
lnames1 = load(file = "PA42_v4.1_dataInputInter.RData");

# Load network data saved in the second part.
lnames2 = load(file = "PA42_v4.1_networkConstructionInter_auto_threshold8_signed.RData");

#Import ka ks results
inputKaks <- read.csv(file="~/PfrenderLab/PA42_v4.1/PA42_v4.1_Olympics_kaksResults.csv")
cleanKaks <- na.omit(inputKaks)
numRow = nrow(cleanKaks)

#Retrieve the dNdS of genes in each module
for(var in 1:numRow)
{
  curGene <- cleanKaks$geneID[var]
  cleanKaks$Color <- ifelse(curGene %in% names(datExprInter[moduleColors,]), moduleColors[which(names(datExprInter[moduleColors,])==curGene)], "None")
}

#Box plot of ka/ks
boxplot(dNdS~Color,data=cleanKaks, main="Module dNdS Ratios",
        xlab="Module Color", ylab="dNdS")