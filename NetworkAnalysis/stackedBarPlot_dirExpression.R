#!/usr/bin/env Rscript

#if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
#BiocManager::install('_______')

#Load the libraries
library(ggplot2)

#The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)

#Set working directory
workingDir = args[1];
setwd(workingDir)

# set counts directory
deDir = args[2]

# retrieve subset tag tag
set <- args[3]

# set the minimum module size
minModSize <- args[4]

# retrieve DE set tag
deSetTag <- args[5]

# set the full subset tag name
tag <- paste(set, minModSize, sep="_")

# Load the expression and trait data saved in the first part
importFile <- paste(set, "dataInput.RData", sep="-")
lnames1 = load(file = importFile)

# Load network data saved in the second part
importFile <- paste(tag, "networkConstruction-stepByStep.RData", sep="-")
lnames2 = load(file = importFile)

#Import DEGs
#geneCountsInter <- read.csv(file="~/PfrenderLab/DEA_PA42_v4.1/glmQLFAnalysis_FDR0.10/glmQLF_2WayANOVA_interaction_topTags.csv")
#SETInterIn <- geneCountsInter[,1:2]
#Get sign of logFC
#SETInterIn$dirExpr <- sign(SETInterIn[,2])

# import DE data for the Olympics
# interaction
inFile <- paste("glmQLF_2WayANOVA", deSetTag, sep="_")
inFile <- paste(inFile, "topTags_LFC1.2.csv", sep="_")
inFile <- paste(deDir, inFile, sep="/")
geneCounts <- read.csv(file=inFile)
#SETIn <- geneCounts[geneCounts$FDR<0.05,1:2]
SETIn <- geneCounts[,1:2]
# get sign of logFC
SETIn$dirExpr <- sign(SETIn[,2])

#Get module color list
colorList = unique(moduleColors)
numRow = length(colorList)*2
colorSets <- data.frame(matrix(ncol = 3, nrow = numRow))

#Retrieve the percent of genes in each module
for(var in 1:length(colorList))
{
  #Get table of logFC signs
  subSet <- SETIn[which(SETIn[,1] %in% names(datExpr)[moduleColors==colorList[var]]),]
  dirTable <- table(subSet[,3])
  #Count number of down regulated genes
  rowNum = var
  numDown <- ifelse(-1 %in% names(dirTable), dirTable[names(dirTable)==-1], 0)
  colorSets[rowNum,1] <- ifelse(numDown==0, 0, numDown/nrow(subSet))
  colorSets[rowNum,2] <- colorList[var]
  colorSets[rowNum,3] <- "Down"
  #Count number of up regulated genes
  rowNum = length(colorList)+var
  numUp <- ifelse(1 %in% names(dirTable), dirTable[names(dirTable)==1], 0)
  colorSets[rowNum,1] <- ifelse(numUp==0, 0, numUp/nrow(subSet))
  colorSets[rowNum,2] <- colorList[var]
  colorSets[rowNum,3] <- "Up"
  #Number of non DE genes, should be zero
  #rowNum = length(colorList)*2+var
  #numZero <- ifelse(0 %in% names(dirTable), dirTable[names(dirTable)==0], 0)
  #colorSets[rowNum,1] <- ifelse(numZero==0, 0, numZero/nrow(subSet))
  #colorSets[rowNum,2] <- colorList[var]
  #colorSets[rowNum,3] <- "Zero"
}

#Set column names
names(colorSets) = c("Percent","Color","Direction")

#Create stacked bar plot
exportFile <- paste("stackedBarPlot", deSetTag, sep="_")
exportFile <- paste(exportFile, "moduleDirExpr.jpg", sep="_")
exportFile <- paste(tag, exportFile, sep="_")
jpeg(file = exportFile, width = 844, height = 596)
colorPlot <- ggplot(colorSets, aes(fill=Direction, y=Percent, x=Color)) + 
  geom_bar(position="stack", stat="identity")
colorPlot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
