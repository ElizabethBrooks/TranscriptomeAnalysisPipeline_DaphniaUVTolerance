#!/usr/bin/env Rscript
#Usage: Rscript statistics_edgeR.r countsFile.csv startColPos endColPos
#Usage Ex: Rscript statistics_exactTest_edgeR.r geneCounts_merged_counted_htseqTophat2_run1_fullset_cleaned.csv 31 36
#R script to perform statistical analysis of gene count tables using
# a DESeq2 two-group comparison 
#Install DESeq2, this should only need to be done once
#bioLite("DESeq2")

#Load the DESeq2 library
library("DESeq2")
#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)
#Test if there is one input argument
if (length(args)!=3) {
  stop("One file name and a range of columns must be supplied.n", call.=FALSE)
}

#Import gene count data
countsTable <- read.csv(file=args[1], row.names="gene")[ ,args[2]:args[3]]
head(countsTable)
#Add grouping factor
group <- factor(c(rep("ctrl",3),rep("treat",3)))

#TO DO
## Example 1: two-group comparison
dds <- makeExampleDESeqDataSet(m=4)
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition","B","A"))