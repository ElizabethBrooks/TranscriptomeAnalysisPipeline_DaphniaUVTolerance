#!/usr/bin/env Rscript
# usage: Rscript ANOVA_OlympicsGenotypes_aov.r workingDir pValues
# R script to perform statistical analysis of gene count tables using aov
# note: https://www.r-bloggers.com/2022/05/two-way-anova-example-in-r-quick-guide/

# turn off scientific notation
options(scipen=999)

# retrieve input file name of gene counts
#args = commandArgs(trailingOnly=TRUE)

# set working directory
#workingDir = args[1];
workingDir = "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCNA_DEGenotypes"
setwd(workingDir)

# import p values
#inputTable <- read.csv(args[2], row.names="module")
inputTable <- read.csv("/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCNA_DEGenotypes/OLYM_WGCNA_aov_pairwise_summary_pValues.csv", row.names="module")

# generate FDR adjusted p values for each module
outputTable <- inputTable
for (row in 1:nrow(inputTable)) {
  outputTable[row,1:3] <- p.adjust(inputTable[row,1:3], method = "fdr", n = length(inputTable[row,1:3]))
}

# write test statistics to a file
write.csv(outputTable, file="OLYM_WGCNA_aov_pairwise_summary_fdr.csv", row.names=TRUE)

