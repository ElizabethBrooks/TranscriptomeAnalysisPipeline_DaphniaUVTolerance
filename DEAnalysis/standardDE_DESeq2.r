#!/usr/bin/env Rscript
#Usage: Rscript standardDE_DESeq2.r countsFile factorGroupingFile
#Usage Ex: Rscript standardDE_DESeq2.r rawCounts.csv expDesign.csv
#R script to perform DE analysis of raw gene counts using DESeq2

#Install packages, this should only need to be done once
install.packages("tidyverse")
install.packages("ggplot2")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("apeglm")
BiocManager::install("vsn")
BiocManager::install("hexbin")

#after you've installed packages:
library("tidyverse")
library("DESeq2")
library("vsn")
library("hexbin")
library("ggplot2")

#Import gene count data
countsTable <- read.csv(file="/home/mae/Documents/BIOS60132/ND_ICG_2020_Practical_Eight_Data/rawCounts.csv", row.names="gene")
#head(countsTable)
#Import grouping factor
targets <- read.csv(file="/home/mae/Documents/BIOS60132/ND_ICG_2020_Practical_Eight_Data/expDesign.csv", row.names="sample")

#Construct a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData=countsTable,
                              colData=targets,
                              design= ~ genotype + treatment)

#Pre-filter to reduce data set size
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#Specify the reference level
dds$condition <- relevel(dds$treatment, ref="VIS")

#Perform standard DESeq2 DE analysis
dds <- DESeq(dds)

#List the coefficients
resultsNames(dds)

#Retrieve the log2 fold change and Wald test
# p value for the last variable in the design
results <- results(dds, name="treatment_VIS_vs_UV")
#Order results table by the smallest p value
resultsOrdered <- results[order(results$pvalue),]
#View top 10 unfiltered results
head(resultsOrdered, 10)
#Generate MA plot of unfiltered results
jpeg("plotMA_results.jpg")
plotMA(results, main="UVR Treatment Effect - E. Brooks")
dev.off()

#Shrink LFCs associated with treatment
resultsShrunk <- lfcShrink(dds, coef="treatment_VIS_vs_UV", type="apeglm")
#Order shrunk results table by the smallest p value
resultsShrunkOrdered <- resultsShrunk[order(resultsShrunk$pvalue),]
#View top 10 shrunk results
head(resultsShrunkOrdered, 10)
#Generate MA plot of shrunk LFC estimates
jpeg("plotMA_resultsShrunk.jpg")
plotMA(resultsShrunk, main="Shrunk LFC Estimates for Treatment - E. Brooks")
dev.off()

#Filter results for a p-value of < 0.05
resultsPV <- subset(results, padj < 0.05)
#Further filter results for a log2FC > 1
resultsLFC <- subset(resultsPV, abs(log2FoldChange) > 1) 
#View top 10 filtered results
head(resultsLFC, 10)

#Plot effects of transformations on the variance
# this gives log2(n + 1)
ntd <- normTransform(dds)
msd <- meanSdPlot(assay(ntd))
#Output mean SD plot with title
jpeg("plotMeanSD.jpg")
msd$gg + ggtitle("Effects of Transformations on Variance - E. Brooks")
dev.off()
