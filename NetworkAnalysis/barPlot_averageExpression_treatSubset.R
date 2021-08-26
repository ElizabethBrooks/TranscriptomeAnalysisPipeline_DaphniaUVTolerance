#if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
#BiocManager::install('_______')

#Set working directory
workingDir="~/PfrenderLab/WGCNA_PA42_v4.1"
setwd(workingDir); 

#Load the libraries
library(edgeR)
library(ggplot2)

# Load the expression and trait data saved in the first part
lnames1 = load(file = "PA42_v4.1_dataInputTreat.RData");

# Load network data saved in the second part.
lnames2 = load(file = "PA42_v4.1_networkConstructionTreat_auto_threshold8_signed.RData");

#Import gene count data for the Olympics
countsTable <- read.csv(file="~/PfrenderLab/PA42_v4.1/geneCounts_cleaned_PA42_v4.1.csv", row.names="gene")[ ,1:24]
#Import grouping factor
targets <- read.csv(file="~/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_Olympics.csv", row.names="sample")

#Setup a design matrix
group <- factor(paste(targets$treatment,targets$genotype,sep="."))
#Create DGE list object
list <- DGEList(counts=countsTable,group=group)
colnames(list) <- rownames(targets)

#Retain genes only if it is expressed at a minimum level
keep <- filterByExpr(list)
list <- list[keep, , keep.lib.sizes=FALSE]

#Use TMM normalization to eliminate composition biases between libraries
list <- calcNormFactors(list)
normList <- cpm(list, normalized.lib.sizes=TRUE)

#Subset the treatment and control samples
SETVIS <- normList[,grepl("VIS",colnames(normList))]
SETUV <- normList[,grepl("UV",colnames(normList))]
#Subset the tolerant and non-tolerant samples
SETTolVIS <- SETVIS[,grepl("Y05|Y023",colnames(SETVIS))]
SETTolUV <- SETUV[,grepl("Y05|Y023",colnames(SETUV))]
SETNTolVIS <- SETVIS[,grepl("E05|R2",colnames(SETVIS))]
SETNTolUV <- SETUV[,grepl("E05|R2",colnames(SETUV))]

#Calculate the means of each gene across samples
meanSETTolVIS <- rowMeans(SETTolVIS)
meanSETTolUV <- rowMeans(SETTolUV)
meanSETNTolVIS <- rowMeans(SETNTolVIS)
meanSETNTolUV <- rowMeans(SETNTolUV)

#Get module color list
colorList = unique(moduleColors)
#Set up data frame for mean values
colorSets <- data.frame(Mean=double(), Color=character(), Sample=character(), stringsAsFactors=FALSE)

#Retrieve the percent of genes in each module
for(var in 1:length(colorList))
{
  #Prepare and add tolerant VIS data
  meanTolVIS <- as.data.frame(unname(meanSETTolVIS[names(datExprTreat[which(names(datExprTreat)[moduleColors==colorList[var]] %in% names(meanSETTolVIS))])]))
  if(length(meanTolVIS) > 0){
    namesTolVIS <- as.data.frame(rep(colorList[var], length(meanTolVIS)))
    sampleTolVIS <- as.data.frame(rep("Tolerant&VIS", length(meanTolVIS)))
    tagsTolVIS <- cbind(namesTolVIS, sampleTolVIS)
    colorSetTolVIS <- cbind(meanTolVIS, tagsTolVIS)
    names(colorSetTolVIS) <- c("Mean", "Color", "Sample")
    colorSets <- rbind(colorSets,colorSetTolVIS)
  }
  #Prepare and add tolerant UV data
  meanTolUV <- as.data.frame(unname(meanSETTolUV[names(datExprTreat[which(names(datExprTreat)[moduleColors==colorList[var]] %in% names(meanSETTolUV))])]))
  if(length(meanTolUV) > 0){
    namesTolUV <- as.data.frame(rep(colorList[var], length(meanTolUV)))
    sampleTolUV <- as.data.frame(rep("Tolerant&UV", length(meanTolUV)))
    tagsTolUV <- cbind(namesTolUV, sampleTolUV)
    colorSetTolUV <- cbind(meanTolUV, tagsTolUV)
    names(colorSetTolUV) <- c("Mean", "Color", "Sample")
    colorSets <- rbind(colorSets,colorSetTolUV)
  }
  #Prepare and add tolerant VIS genotype data
  meanNTolVIS <- as.data.frame(unname(meanSETNTolVIS[names(datExprTreat[which(names(datExprTreat)[moduleColors==colorList[var]] %in% names(meanSETNTolVIS))])]))
  if(length(meanNTolVIS) > 0){
    namesNTolVIS <- as.data.frame(rep(colorList[var], length(meanNTolVIS)))
    sampleNTolVIS <- as.data.frame(rep("Sensitive&VIS", length(meanNTolVIS)))
    tagsNTolVIS <- cbind(namesNTolVIS, sampleNTolVIS)
    colorSetNTolVIS <- cbind(meanNTolVIS, tagsNTolVIS)
    names(colorSetNTolVIS) <- c("Mean", "Color", "Sample")
    colorSets <- rbind(colorSets,colorSetNTolVIS)
  }
  #Prepare and add not tolerant UV genotype data
  meanNTolUV <- as.data.frame(unname(meanSETNTolUV[names(datExprTreat[which(names(datExprTreat)[moduleColors==colorList[var]] %in% names(meanSETNTolUV))])]))
  if(length(meanNTolUV) > 0){
    namesNTolUV <- as.data.frame(rep(colorList[var], length(meanNTolUV)))
    sampleNTolUV <- as.data.frame(rep("Sensitive&UV", length(meanNTolUV)))
    tagsNTolUV <- cbind(namesNTolUV, sampleNTolUV)
    colorSetNTolUV <- cbind(meanNTolUV, tagsNTolUV)
    names(colorSetNTolUV) <- c("Mean", "Color", "Sample")
    colorSets <- rbind(colorSets,colorSetNTolUV)
  }
}

#Subset by sample type
colorSetTolVIS <- colorSets[colorSets$Sample=="Tolerant&VIS",]
colorSetTolUV <- colorSets[colorSets$Sample=="Tolerant&UV",]
colorSetNTolVIS <- colorSets[colorSets$Sample=="Sensitive&UV",]
colorSetNTolUV <- colorSets[colorSets$Sample=="Sensitive&VIS",]

#Generate box plot of median treatment values
#p<-ggplot(rbind(colorSetTolVIS,colorSetTolUV), aes(x=Color, y=Mean, fill=Sample)) +
#  geom_boxplot(position=position_dodge(1))
#p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#Generate box plot of median tolerance values
#p<-ggplot(rbind(colorSetNTolVIS,colorSetNTolUV), aes(x=Color, y=Mean, fill=Sample)) +
#  geom_boxplot(position=position_dodge(1))
#p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#Generate a line plot of median treatment values
#ggplot(rbind(colorSetTolVIS,colorSetTolUV), aes(x=Sample, y=Mean, colour=Color, group=Color)) + geom_line()

#Grouped barplot of median tolerant values
p <- ggplot(rbind(colorSetTolVIS,colorSetTolUV), aes(x=Color, y=Mean, fill=Sample)) +
  geom_bar(stat = "identity", position = "dodge")
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#Grouped barplot of median not tolerants values
p <- ggplot(rbind(colorSetNTolVIS,colorSetNTolUV), aes(x=Color, y=Mean, fill=Sample)) +
  geom_bar(stat = "identity", position = "dodge")
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#df1  = transform(df, mean=rowMeans(df[cols]), sd=apply(df[cols],1, sd))
#ggplot(df1, aes(x=as.factor(Gene), y=mean, fill=Species)) +
#geom_bar(position=position_dodge(), stat="identity", colour='black') +
#  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,position=position_dodge(.9))
