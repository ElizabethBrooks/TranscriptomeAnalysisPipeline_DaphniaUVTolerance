#if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
#BiocManager::install('_______')

#Set working directory
workingDir="~/PfrenderLab/WGCNA_PA42_v4.1"
setwd(workingDir); 

#Load the libraries
library(ggplot2)

# Load the expression and trait data saved in the first part
lnames1 = load(file = "PA42_v4.1_dataInputTol.RData");

# Load network data saved in the second part.
lnames2 = load(file = "PA42_v4.1_networkConstructionTol_auto_threshold8_signedNowick.RData");

#Import DEGs
#geneCountsTol <- read.csv(file=args[1])
geneCountsTol <- read.csv(file="~/PfrenderLab/DEA_PA42_v4.1/glmQLFAnalysis_FDR0.10/glmQLF_2WayANOVA_TvsN_topTags.csv")
SETTolIn <- geneCountsTol[,1:2]
#Get sign of logFC
SETTolIn$dirExpr <- sign(SETTolIn[,2])

#Get module color list
colorList = unique(moduleColors)
numRow = length(colorList)*2
colorSets <- data.frame(matrix(ncol = 3, nrow = numRow))

#Retrieve the percent of genes in each module
for(var in 1:length(colorList))
{
  #Get table of logFC signs
  subSet <- SETTolIn[which(SETTolIn[,1] %in% names(datExprTol)[moduleColors==colorList[var]]),]
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
jpeg("stackedBarPlotTol_moduleDirExpression_signedNowick.jpg", width = 844, height = 596)
colorPlot <- ggplot(colorSets, aes(fill=Direction, y=Percent, x=Color)) + 
  geom_bar(position="stack", stat="identity")
colorPlot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
