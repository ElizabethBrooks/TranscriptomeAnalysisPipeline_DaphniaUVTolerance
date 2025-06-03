#!/usr/bin/env Rscript
#Usage: Rscript normalizeFilteredCounts.r
#Usage Ex: Rscript normalizeFilteredCounts.r
#R script to normalized and log transform filtered gene counts

# Load libraries
library("edgeR")
library("limma")

#Import gene count data
#countsTable <- read.csv(file=args[1], row.names="gene")[ ,args[2]:args[3]]
countsTable <- read.csv(file="/Users/bamflappy/PfrenderLab/PA42_v4.1/geneCounts_cleaned_PA42_v4.1.csv", row.names="gene")[ ,1:24]
#Import grouping factor
targets <- read.csv(file="/Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_binned_Olympics.csv", row.names="sample")

#Setup a design matrix
group <- factor(paste(targets$treatment,targets$genotype,sep="."))
#Create DGE list object
list <- DGEList(counts=countsTable,group=group)
colnames(list) <- rownames(targets)

#Retain genes only if it is expressed at a minimum level
keep <- filterByExpr(list)
summary(keep)
list <- list[keep, , keep.lib.sizes=FALSE]

#Use TMM normalization to eliminate composition biases
# between libraries
list <- calcNormFactors(list)
normList <- cpm(list, normalized.lib.sizes=TRUE)

#Write normalized counts to file
write.table(normList, file="/Users/bamflappy/PfrenderLab/PA42_v4.1/PA42_v4.1_normalizedCountsOlympics.csv", sep=",", row.names=TRUE)

#Log transform normalized counts
#Use TMM normalization to eliminate composition biases
# between libraries
list <- calcNormFactors(list)
normList <- cpm(list, log=TRUE, normalized.lib.sizes=TRUE)

#Write log transformed normalized counts to file
write.table(normList, file="/Users/bamflappy/PfrenderLab/PA42_v4.1/PA42_v4.1_normalizedLogCountsOlympics.csv", sep=",", row.names=TRUE)
