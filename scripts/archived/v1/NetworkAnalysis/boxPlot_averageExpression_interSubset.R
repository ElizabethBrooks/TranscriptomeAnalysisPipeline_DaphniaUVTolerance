#if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
#BiocManager::install('_______')

#Load the libraries
library(edgeR)
library(ggplot2)

# Load the expression and trait data saved in the first part
lnames1 = load(file = "PA42_v4.1_dataInputTol.RData");

# Load network data saved in the second part.
lnames2 = load(file = "PA42_v4.1_networkConstructionTol_auto_threshold8_signed.RData");

#Import gene count data for the Olympics
countsTable <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/geneCounts_cleaned_PA42_v4.1.csv", row.names="gene")[ ,1:24]
#Import grouping factor
targets <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/expDesign_Olympics.csv", row.names="sample")

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
SETTol <- normList[,grepl("Y05|Y023",colnames(normList))]
SETNTol <- normList[,grepl("E05|R2",colnames(normList))]

#Calculate the means of each gene across samples
meanSETVIS <- rowMeans(SETVIS)
meanSETUV <- rowMeans(SETUV)
#Calculate the means of each gene across samples
meanSETTol <- rowMeans(SETTol)
meanSETNTol <- rowMeans(SETNTol)

#Get module color list
colorList = unique(moduleColors)
#Set up data frame for mean values
colorSets <- data.frame(Mean=double(), Color=character(), Sample=character(), stringsAsFactors=FALSE)

#Retrieve the percent of genes in each module
for(var in 1:length(colorList))
{
  #Prepare and add VIS data
  meanVIS <- as.data.frame(unname(meanSETVIS[names(datExprTol[which(names(datExprTol)[moduleColors==colorList[var]] %in% names(meanSETVIS))])]))
  if(length(meanVIS) > 0){
    namesVIS <- as.data.frame(rep(colorList[var], length(meanVIS)))
    sampleVIS <- as.data.frame(rep("VIS", length(meanVIS)))
    tagsVIS <- cbind(namesVIS, sampleVIS)
    colorSetVIS <- cbind(meanVIS, tagsVIS)
    names(colorSetVIS) <- c("Mean", "Color", "Sample")
    colorSets <- rbind(colorSets,colorSetVIS)
  }
  #Prepare and add UV data
  meanUV <- as.data.frame(unname(meanSETUV[names(datExprTol[which(names(datExprTol)[moduleColors==colorList[var]] %in% names(meanSETUV))])]))
  if(length(meanUV) > 0){
    namesUV <- as.data.frame(rep(colorList[var], length(meanUV)))
    sampleUV <- as.data.frame(rep("UV", length(meanUV)))
    tagsUV <- cbind(namesUV, sampleUV)
    colorSetUV <- cbind(meanUV, tagsUV)
    names(colorSetUV) <- c("Mean", "Color", "Sample")
    colorSets <- rbind(colorSets,colorSetUV)
  }
  #Prepare and add tolerant genotype data
  meanTol <- as.data.frame(unname(meanSETTol[names(datExprTol[which(names(datExprTol)[moduleColors==colorList[var]] %in% names(meanSETTol))])]))
  if(length(meanTol) > 0){
    namesTol <- as.data.frame(rep(colorList[var], length(meanTol)))
    sampleTol <- as.data.frame(rep("Tol", length(meanTol)))
    tagsTol <- cbind(namesTol, sampleTol)
    colorSetTol <- cbind(meanTol, tagsTol)
    names(colorSetTol) <- c("Mean", "Color", "Sample")
    colorSets <- rbind(colorSets,colorSetTol)
  }
  #Prepare and add not tolerant genotype data
  meanNTol <- as.data.frame(unname(meanSETNTol[names(datExprTol[which(names(datExprTol)[moduleColors==colorList[var]] %in% names(meanSETNTol))])]))
  if(length(meanNTol) > 0){
    namesNTol <- as.data.frame(rep(colorList[var], length(meanNTol)))
    sampleNTol <- as.data.frame(rep("NTol", length(meanNTol)))
    tagsNTol <- cbind(namesNTol, sampleNTol)
    colorSetNTol <- cbind(meanNTol, tagsNTol)
    names(colorSetNTol) <- c("Mean", "Color", "Sample")
    colorSets <- rbind(colorSets,colorSetNTol)
  }
}

#Subset by sample type
colorSetVIS <- colorSets[colorSets$Sample=="VIS",]
colorSetUV <- colorSets[colorSets$Sample=="UV",]
colorSetTol <- colorSets[colorSets$Sample=="Tol",]
colorSetNTol <- colorSets[colorSets$Sample=="NTol",]

#Generate box plot of median treatment values
p<-ggplot(rbind(colorSetVIS,colorSetUV), aes(x=Color, y=Mean, fill=Sample)) +
  geom_boxplot(position=position_dodge(1))
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#Generate box plot of median tolerance values
p<-ggplot(rbind(colorSetTol,colorSetNTol), aes(x=Color, y=Mean, fill=Sample)) +
  geom_boxplot(position=position_dodge(1))
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
