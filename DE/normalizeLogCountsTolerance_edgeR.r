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

#Set working directory
workingDir = args[1];
setwd(workingDir)

#Import gene count data
inputTable <- read.csv(file=args[2], row.names="gene")[ ,args[3]:args[4]]

#Trim the data table
countsTable <- head(inputTable, - 5)

#Import grouping factor
targets <- read.csv(file=args[5], row.names="sample")

#Setup a design matrix
group <- factor(paste(targets$treatment,targets$tolerance,sep="."))
#cbind(targets,Group=group)
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
#list$samples
#Write log transformed normalized counts to file
normList <- cpm(list, normalized.lib.sizes=TRUE, log=TRUE)
write.table(normList, file="normalizedCounts_logTransformed.csv", sep=",", row.names=TRUE, quote=FALSE)
