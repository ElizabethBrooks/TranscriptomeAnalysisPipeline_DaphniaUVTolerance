#!/usr/bin/env Rscript
#Usage: Rscript statistics_edgeR.r countsFile.csv startColPos endColPos
#Usage Ex: Rscript statistics_edgeR.r ../GeneCounts_Merged/merged_counts_legacy_cleaned.csv 1 6
#R script to perform statistical analysis of gene count tables using edgeR
#Install edgeR, this should only need to be done once
#Since edgeR is already installed on the CRC this can be skipped if using the module
#bioLite("edgeR")
#Load the edgeR library

#Import gene count data
x <- read.delim("TableOfCounts.txt",row.names="Symbol")
#Add grouping factor
group <- factor(c(1,1,2,2))
#Generate list of DE genes
y <- DGEList(counts=x,group=group)
#Calculate normalized factors
y <- calcNormFactors(y)
y$samples
#Design the model
design <- model.matrix(~group)
#Estimate pairwise disperssion
y <- estimateDisp(y,design)

#Perform quasi-likelihood F-tests
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
topTags(qlf)