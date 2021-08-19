#Set working directory
#workingDir = args[1];
workingDir="/home/mae/Documents/RNASeq_Workshop_ND/WGCNA_PA42_v4.1/effectSubsets"
setwd(workingDir); 

# Load libraries
library(ggplot2)

# Load the expression and trait data saved in the first part
lnames1 = load(file = "PA42_v4.1_dataInputTreat.RData");

# Load network data saved in the second part.
lnames2 = load(file = "PA42_v4.1_networkConstructionTreat_auto_threshold8_signed.RData");

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
  #Determine the sets and intersections of the effects
  SETInter <- which(names(datExprTreat)[moduleColors==colorList[var]] %in% SETInterIn)
  SETTreat <- which(names(datExprTreat)[moduleColors==colorList[var]] %in% SETTreatIn)
  SETTol <- which(names(datExprTreat)[moduleColors==colorList[var]] %in% SETTolIn)
  SETIntersectAll <- intersect(intersect(SETInter,SETTreat),SETTol)
  SETInterTreat <- setdiff(intersect(SETInter,SETTreat),SETIntersectAll)
  SETInterTol <- setdiff(intersect(SETInter,SETTol),SETIntersectAll)
  SETTreatTol <- setdiff(intersect(SETTreat,SETTol),SETIntersectAll)
  #Filter intersections from the effects sets
  SETInter <- setdiff(setdiff(setdiff(SETInter,SETInterTreat),SETInterTol),SETIntersectAll)
  SETTreat <- setdiff(setdiff(setdiff(SETTreat,SETInterTreat),SETTreatTol),SETIntersectAll)
  SETTol <- setdiff(setdiff(setdiff(SETTol,SETTreatTol),SETInterTol),SETIntersectAll)
  
  #Add interaction data to first block
  rowNum = var
  colorSets[rowNum,1] <- colorList[var]
  colorSets[rowNum,2] <- length(SETInter)/length(names(datExprTreat)[moduleColors==colorList[var]])*100
  colorSets[rowNum,3] <- "Interaction"
  
  #Add treatment data to second block
  rowNum = length(colorList)+var
  colorSets[rowNum,1] <- colorList[var]
  colorSets[rowNum,2] <- length(SETTreat)/length(names(datExprTreat)[moduleColors==colorList[var]])*100
  colorSets[rowNum,3] <- "Treatment"
  
  #Add tolerance data to third block
  rowNum = length(colorList)*2+var
  colorSets[rowNum,1] <- colorList[var]
  colorSets[rowNum,2] <- length(SETTol)/length(names(datExprTreat)[moduleColors==colorList[var]])*100
  colorSets[rowNum,3] <- "Tolerance"
  
  #Add intersection of all data to third block
  rowNum = length(colorList)*3+var
  colorSets[rowNum,1] <- colorList[var]
  colorSets[rowNum,2] <- length(SETIntersectAll)/length(names(datExprTreat)[moduleColors==colorList[var]])*100
  colorSets[rowNum,3] <- "IntersectAll"
  
  #Add intersection of interaction and treatment data to third block
  rowNum = length(colorList)*4+var
  colorSets[rowNum,1] <- colorList[var]
  colorSets[rowNum,2] <- length(SETInterTreat)/length(names(datExprTreat)[moduleColors==colorList[var]])*100
  colorSets[rowNum,3] <- "Treatment&Interaction"
  
  #Add intersection of interaction and tolerance data to third block
  rowNum = length(colorList)*5+var
  colorSets[rowNum,1] <- colorList[var]
  colorSets[rowNum,2] <- length(SETInterTol)/length(names(datExprTreat)[moduleColors==colorList[var]])*100
  colorSets[rowNum,3] <- "Tolerance&Interaction"
  
  #Add intersection of treatment and tolerance data to third block
  rowNum = length(colorList)*6+var
  colorSets[rowNum,1] <- colorList[var]
  colorSets[rowNum,2] <- length(SETTreatTol)/length(names(datExprTreat)[moduleColors==colorList[var]])*100
  colorSets[rowNum,3] <- "Treatment&Tolerance"
}

#Set column names
names(colorSets) = c("Color","Percent","Effect")

#Create stacked bar plot
jpeg("stackedBarPlotTreat_moduleEffectSubsets_filtered.jpg", width = 844, height = 596)
colorPlot <- ggplot(colorSets, aes(fill=Effect, y=Percent, x=Color)) + 
  geom_bar(position="stack", stat="identity")
colorPlot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
