#!/usr/bin/env Rscript
#Usage: Rscript statistics_edgeR.r countsFile.csv startColPos endColPos
#Usage Ex: Rscript statistics_edgeR.r ../GeneCounts_Merged/merged_counts_legacy_cleaned.csv 1 6
#R script to perform statistical analysis of gene count tables using edgeR
#Install edgeR, this should only need to be done once
#Since edgeR is already installed on the CRC this can be skipped if using the module
#bioLite("edgeR")
#Load the edgeR library

#GLM Method
#Import gene count data
countsTable <- read.csv(file=args[1], row.names="gene")[ ,args[2]:args[3]]
head(countsTable)
#Add grouping factor
group <- factor(c(rep("ctrl",3),rep("treat",3)))
#Create DGE list object
list <- DGEList(counts=countsTable,group=group)
#Calculate normalized factors
list <- calcNormFactors(list)
list$samples
#Design the model
design <- model.matrix(~group)
#Estimate common dispersion, trended dispersions, and tagwise dispersions
list <- estimateDisp(list,design)

#Perform quasi-likelihood F-tests
fit <- glmQLFit(list,design)
qlf <- glmQLFTest(fit,coef=2)
topTags(qlf)