#!/usr/bin/env Rscript
#Usage: Rscript normalize_edgeR.r countsFile factorGroupingFile
#Usage Ex: Rscript normalize_edgeR.r cleaned.csv expDesign_fullSet.csv
#R script to perform statistical analysis of gene count tables using edgeR two way ANOVA

#Install edgeR and statmod, this should only need to be done once
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("edgeR")
#install.packages("statmod")

#Load the edgeR library
library("edgeR")

#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)
#Test if there is one input argument
if (length(args)!=2) {
  stop("Two file names must be supplied.n", call.=FALSE)
}

#Import gene count data
countsTable <- read.csv(file=args[1], row.names="gene")
head(countsTable)
#Import grouping factor
targets <- read.csv(file=args[2], row.names="sample")

#Setup a design matrix
group <- factor(paste(targets$treatment,targets$genotype,sep="."))
#cbind(targets,Group=group)
#Create DGE list object
list <- DGEList(counts=countsTable,group=group)
colnames(list) <- targets$sample

#Retain genes only if it is expressed at a minimum level
keep <- filterByExpr(list)
summary(keep)
list <- list[keep, , keep.lib.sizes=FALSE]

#Use TMM normalization to eliminate composition biases
# between libraries
list <- calcNormFactors(list)
#list$samples
#Write normalized counts to file
normList <- cpm(list, normalized.lib.sizes=TRUE)
write.table(normList, file="normalizedCounts.csv", sep=",", row.names=TRUE)