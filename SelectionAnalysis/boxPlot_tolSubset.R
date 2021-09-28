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

#Setup results table
resultsTbl <- as.data.frame(matrix(nrow=numRow,ncol=3))
colnames(resultsTbl) <- c("geneID","Color","dNdS")

#Retrieve the dNdS of genes in each module
for(var in 1:numRow)
{
  #curGene <- filtKaks$geneID[var]
  curGene <- names(datExprTol[moduleColors,])[var]
  #filtKaks$Color[var] <- ifelse(curGene %in% names(datExprTol[moduleColors,]), moduleColors[which(names(datExprTol[moduleColors,])==curGene)], "None")
  resultsTbl$geneID[var] <- curGene
  resultsTbl$Color[var] <- ifelse(curGene %in% filtKaks$geneID, moduleColors[which(names(datExprTol[moduleColors,])==curGene)], "NA")
  resultsTbl$dNdS[var] <- ifelse(curGene %in% filtKaks$geneID, filtKaks[which(filtKaks$geneID==curGene),5], "NA")
}

#Remove NA entries
cleanResults <- transform(resultsTbl, dNdS = as.numeric(dNdS))
cleanResults <- na.omit(cleanResults)

#Box plot of ka/ks
#boxplot(dNdS~Color,data=cleanResults, main="Module dNdS Ratios",
#        xlab="Module Color", ylab="dNdS")
jpeg("SelectionAnalysis/boxPlotTol_filtered_dNdS_signedNowick.jpg", width = 844, height = 596)
p <- ggplot(cleanResults, aes(x=Color, y=dNdS)) + 
  geom_boxplot(notch=TRUE)
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p + geom_hline(yintercept=1, linetype="dashed", color = "red")
dev.off()

