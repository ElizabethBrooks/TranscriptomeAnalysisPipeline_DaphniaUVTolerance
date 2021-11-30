#!/usr/bin/env Rscript
#Usage: Rscript geneSetTest_camera.r countsFile startColumn endColumn factorGroupingFile
#Usage Ex: Rscript geneSetTest_camera.r cleaned.csv 1 24 expDesign_binned_Olympics.csv
#R script to perform gene set enrichment testing using camera

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

#The experimental design is parametrized with a one-way layout, 
# where one coefficient is assigned to each group
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# First test set
index1 <- 1:20

# Second test set
index2 <- 21:40

#Individual gene set tests
#Gene set enrichment test for the first set
camera(normList, index1, design)
#Gene set enrichment test for the second set
camera(normList, index2, design)

#Batch gene set tests
#Inter-gene correlation will be estimated for each tested set
camera(normList, list(set1=index1,set2=index2), design, inter.gene.cor=NA)
#Specify an inter-gene correlation of 0.01
camera(normList, list(set1=index1,set2=index2), design, inter.gene.cor=0.01)

# Pre-ranked version
fit <- eBayes(lmFit(normList, design))
cameraPR(fit$t[,2], list(set1=index1,set2=index2))
#cameraPR(fit$p.value[,2], list(set1=index1,set2=index2))
cameraPR(fit$F, list(set1=index1,set2=index2))
#cameraPR(fit$F.p.value, list(set1=index1,set2=index2))
