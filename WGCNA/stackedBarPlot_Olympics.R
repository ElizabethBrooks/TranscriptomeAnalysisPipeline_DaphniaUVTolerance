#Set working directory
#workingDir = args[1];
workingDir="/home/mae/Documents/RNASeq_Workshop_ND/WGCNA_PA42_v4.1/allGenes"
setwd(workingDir); 

# Load libraries
library(WGCNA)
library(ggplot2)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

# Load the expression and trait data saved in the first part
lnames = load(file = "PA42_v4.1_dataInput.RData");

# Load network data saved in the second part.
lnames = load(file = "PA42_v4.1_networkConstruction_auto_threshold8.RData");

ddr <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/DDRGOTF_Dmel_PA42_v4.1_combined_geneIDs_uniq.csv")
SETDDR <- ddr[,1]

geneCountsInter <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/geneCounts_mergedHisat2_PA42_v4.1/glmQLFAnalysis_FDR0.10/glmQLF_2WayANOVA_interaction_topTags_filtered.csv")
geneCountsTreat <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/geneCounts_mergedHisat2_PA42_v4.1/glmQLFAnalysis_FDR0.10/glmQLF_2WayANOVA_UVvsVIS_topTags_filtered.csv")
geneCountsTol <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/geneCounts_mergedHisat2_PA42_v4.1/glmQLFAnalysis_FDR0.10/glmQLF_2WayANOVA_TvsN_topTags_filtered.csv")
SETInterIn <- geneCountsInter[,1]
SETTreatIn <- geneCountsTreat[,1]
SETTolIn <- geneCountsTol[,1]


#Get module color list
colorList = unique(moduleColors)
numRow = length(colorList)*3
colorSets <- data.frame(matrix(ncol = 3, nrow = numRow))

#Retrieve the percent of genes in each module
for(var in 1:length(colorList))
{
  #Add interaction data to first block
  rowNum = var
  colorSets[rowNum,1] <- colorList[var]
  colorSets[rowNum,2] <- length(which(names(datExpr)[moduleColors==colorList[var]] %in% SETInterIn))/length(names(datExpr)[moduleColors==colorList[var]])*100
  colorSets[rowNum,3] <- "Inter"
    
  #Add treatment data to second block
  rowNum = length(colorList)+var
  colorSets[rowNum,1] <- colorList[var]
  colorSets[rowNum,2] <- length(which(names(datExpr)[moduleColors==colorList[var]] %in% SETTreatIn))/length(names(datExpr)[moduleColors==colorList[var]])*100
  colorSets[rowNum,3] <- "Treat"
  
  #Add tolerance data to third block
  rowNum = length(colorList)*2+var
  colorSets[rowNum,1] <- colorList[var]
  colorSets[rowNum,2] <- length(which(names(datExpr)[moduleColors==colorList[var]] %in% SETTolIn))/length(names(datExpr)[moduleColors==colorList[var]])*100
  colorSets[rowNum,3] <- "Tol"
}

#Set column names
names(colorSets) = c("color","percent","effect")

#Create stacked bar plot
colorPlot <- ggplot(colorSets, aes(fill=effect, y=percent, x=color)) + 
  geom_bar(position="stack", stat="identity")
colorPlot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

