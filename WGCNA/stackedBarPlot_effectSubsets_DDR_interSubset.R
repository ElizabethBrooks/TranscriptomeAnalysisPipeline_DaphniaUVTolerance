#Set working directory
#workingDir = args[1];
workingDir="/home/mae/Documents/RNASeq_Workshop_ND/WGCNA_PA42_v4.1/effectSubsets"
setwd(workingDir); 

# Load libraries
library(ggplot2)

# Load the expression and trait data saved in the first part
lnames1 = load(file = "PA42_v4.1_dataInputInter.RData");

# Load network data saved in the second part.
lnames2 = load(file = "PA42_v4.1_networkConstructionInter_auto_threshold8_signed.RData");

ddr <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/DDRGOTF_Dmel_PA42_v4.1_combined_geneIDs_uniq.csv")
SETDDR <- ddr[,1]

#Get module color list
colorList = unique(moduleColors)
numRow = length(colorList)
colorSets <- data.frame(matrix(ncol = 2, nrow = numRow))

#Retrieve the percent of genes in each module
for(var in 1:length(colorList))
{
  #Print the color of the module
  #print(colorList[var])
  #Add DDR data for the current module
  numDDR <- which(names(datExprInter)[moduleColors==colorList[var]] %in% SETDDR)
  colorSets[var,1] <- colorList[var]
  colorSets[var,2] <- length(numDDR)
  #Print the number of DDR genes in the current module
  #print(length(numDDR))
}

#Set column names
names(colorSets) = c("Color","Genes")

#Create stacked bar plot
jpeg("stackedBarPlotInter_DDR.jpg", width = 844, height = 596)
colorPlot <- ggplot(colorSets, aes(y=Genes, x=Color)) + 
  geom_bar(position="stack", stat="identity", fill="steelblue")
colorPlot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
