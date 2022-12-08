#!/usr/bin/env Rscript

# Load libraries
library(ggplot2)

#The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)

#Set working directory
workingDir = args[1]
setwd(workingDir);

# set counts directory
deDir = args[2]

# retrieve subset tag tag
set <- args[3]

# set the minimum module size
minModSize <- args[4]

# set the full subset tag name
tag <- paste(set, minModSize, sep="_")

# Load the expression and trait data saved in the first part
importFile <- paste(set, "dataInput.RData", sep="-")
lnames1 = load(file = importFile)

# Load network data saved in the second part
importFile <- paste(tag, "networkConstruction-stepByStep.RData", sep="-")
lnames2 = load(file = importFile)

# import DE data for the Olympics
# interaction
inFile <- paste(deDir, "glmQLF_2WayANOVA_interaction_topTags_LFC1.2.csv", sep="/")
geneCountsInter <- read.csv(file=inFile)
SETInterIn <- geneCountsInter[,1]
# treatment
inFile <- paste(deDir, "glmQLF_2WayANOVA_UVvsVIS_topTags_LFC1.2.csv", sep="/")
geneCountsTreat <- read.csv(file=inFile)
SETTreatIn <- geneCountsTreat[,1]
# tolerance
inFile <- paste(deDir, "glmQLF_2WayANOVA_TvsN_topTags_LFC1.2.csv", sep="/")
geneCountsTol <- read.csv(file=inFile)
SETTolIn <- geneCountsTol[,1]

#Get module color list
colorList = unique(moduleColors)
numRow = length(colorList)*7
colorSets <- data.frame(matrix(ncol = 3, nrow = numRow))

#Retrieve the percent of genes in each module
for(var in 1:length(colorList))
{
  #Determine the sets and intersections of the effects
  SETInter <- which(names(datExpr)[moduleColors==colorList[var]] %in% SETInterIn)
  SETTreat <- which(names(datExpr)[moduleColors==colorList[var]] %in% SETTreatIn)
  SETTol <- which(names(datExpr)[moduleColors==colorList[var]] %in% SETTolIn)
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
  colorSets[rowNum,2] <- length(SETInter)/length(names(datExpr)[moduleColors==colorList[var]])*100
  colorSets[rowNum,3] <- "Interaction"
  #Add treatment data to second block
  rowNum = length(colorList)+var
  colorSets[rowNum,1] <- colorList[var]
  colorSets[rowNum,2] <- length(SETTreat)/length(names(datExpr)[moduleColors==colorList[var]])*100
  colorSets[rowNum,3] <- "Treatment"
  #Add tolerance data to third block
  rowNum = length(colorList)*2+var
  colorSets[rowNum,1] <- colorList[var]
  colorSets[rowNum,2] <- length(SETTol)/length(names(datExpr)[moduleColors==colorList[var]])*100
  colorSets[rowNum,3] <- "Tolerance"
  #Add intersection of all data to third block
  rowNum = length(colorList)*3+var
  colorSets[rowNum,1] <- colorList[var]
  colorSets[rowNum,2] <- length(SETIntersectAll)/length(names(datExpr)[moduleColors==colorList[var]])*100
  colorSets[rowNum,3] <- "IntersectAll"
  #Add intersection of interaction and treatment data to third block
  rowNum = length(colorList)*4+var
  colorSets[rowNum,1] <- colorList[var]
  colorSets[rowNum,2] <- length(SETInterTreat)/length(names(datExpr)[moduleColors==colorList[var]])*100
  colorSets[rowNum,3] <- "Treatment&Interaction"
  #Add intersection of interaction and tolerance data to third block
  rowNum = length(colorList)*5+var
  colorSets[rowNum,1] <- colorList[var]
  colorSets[rowNum,2] <- length(SETInterTol)/length(names(datExpr)[moduleColors==colorList[var]])*100
  colorSets[rowNum,3] <- "Tolerance&Interaction"
  #Add intersection of treatment and tolerance data to third block
  rowNum = length(colorList)*6+var
  colorSets[rowNum,1] <- colorList[var]
  colorSets[rowNum,2] <- length(SETTreatTol)/length(names(datExpr)[moduleColors==colorList[var]])*100
  colorSets[rowNum,3] <- "Treatment&Tolerance"
}

#Set column names
names(colorSets) = c("Color","Percent","Effect")

#Create stacked bar plot
exportFile <- paste(tag, "stackedBarPlot_moduleEffectSubsets.jpg", sep="_")
jpeg(file = exportFile, width = 844, height = 596)
colorPlot <- ggplot(colorSets, aes(fill=Effect, y=Percent, x=Color)) + 
  geom_bar(position="stack", stat="identity")
colorPlot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
