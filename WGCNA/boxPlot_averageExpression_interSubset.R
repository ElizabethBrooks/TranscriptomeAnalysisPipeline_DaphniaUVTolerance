#if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
#BiocManager::install('_______')

#Load the libraries
library(edgeR)
library(ggplot2)

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

#Calculate the means of each gene across samples
meanSETVIS <- rowMeans(SETVIS)
meanSETUV <- rowMeans(SETUV)
