#if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
#BiocManager::install('_______')

#Set working directory
workingDir="~/PfrenderLab/DEA_PA42_v4.1"
setwd(workingDir); 

#Load the libraries
library(ggplot2)

#Import DEGs
#geneCountsInter <- read.csv(file=args[1])
#geneCountsTreat <- read.csv(file=args[2])
#geneCountsTol <- read.csv(file=args[3])
geneCountsInter <- read.csv(file="glmQLFAnalysis_FDR0.10/glmQLF_2WayANOVA_interaction_topTags.csv")
geneCountsTreat <- read.csv(file="glmQLFAnalysis_FDR0.10/glmQLF_2WayANOVA_UVvsVIS_topTags.csv")
geneCountsTol <- read.csv(file="glmQLFAnalysis_FDR0.10/glmQLF_2WayANOVA_TvsN_topTags.csv")
SETInterIn <- geneCountsInter[,1]
SETTreatIn <- geneCountsTreat[,1]
SETTolIn <- geneCountsTol[,1]

#Get module color list
colorList = unique(moduleColors)
numRow = length(colorList)
colorSets <- data.frame(matrix(ncol = 3, nrow = numRow))

#Get sign of logFC
sign(x)

#Retrieve the percent of genes in each module
for(var in 1:length(colorList))
{
  
}

#Set column names
names(colorSets) = c("Positive","Negative","Zero")

#Create stacked bar plot
jpeg("stackedBarPlotInter_moduleEffectSubsets.jpg", width = 844, height = 596)
colorPlot <- ggplot(colorSets, aes(fill=Effect, y=Percent, x=Color)) + 
  geom_bar(position="stack", stat="identity")
colorPlot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()