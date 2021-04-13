#!/usr/bin/env Rscript
#Usage: Rscript exactTest_edgeR.r countsFile.csv startColPos endColPos FDR genotype
#Usage Ex: Rscript exactTest_edgeR.r genome_sortedName_samtoolsHisat2_run1_counted_htseq_run1 31 36 0.10 Sierra
#R script to perform statistical analysis of gene count tables using edgeR exact test

#Install edgeR, this should only need to be done once
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("edgeR")

#Load the edgeR library
library("edgeR")

#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)

#Import gene count data
countsTable <- read.csv(file=args[1], row.names="gene")[ ,args[2]:args[3]]
#Add grouping factor
group <- factor(c(rep("ctrl",3),rep("treat",3)))
#Create DGE list object
list <- DGEList(counts=countsTable,group=group)
#Retrieve input genotype
genotype <- args[4]

#There is no purpose in analysing genes that are not expressed in either 
# experimental condition, so genes are first filtered on expression levels
keep <- filterByExpr(list)
list <- list[keep, , keep.lib.sizes=FALSE]
#Calculate normalized factors
list <- calcNormFactors(list)
#Generate normalized counts
normList <- cpm(list, normalized.lib.sizes=TRUE)

#Produce a matrix of pseudo-counts
#Estimate common dispersion and tagwise dispersions
list <- estimateDisp(list)

#Perform an exact test for treat vs ctrl
tested <- exactTest(list, pair=c("ctrl", "treat"))
write.table(tested, file=paste("exactTest", genotype, "pvalues.csv", sep = "_"), sep=",", row.names=TRUE)