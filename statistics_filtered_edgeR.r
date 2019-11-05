#!/usr/bin/env Rscript
#R script to perform statistical analysis of filtered gene count tables using edgeR
#Usage: Rscript statistics_edgeR.r countsFile.csv startColPos endColPos
#Usage Ex: Rscript statistics_edgeR.r ../GeneCounts_Merged/merged_counts_legacy_cleaned.csv 1 6

#Install edgeR, this should only need to be done once
#Since edgeR is already installed on the CRC this can be skipped if using the module
#bioLite("edgeR")
#Load the edgeR library
library("edgeR")
#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)
#Test if there is one input argument
if (length(args)!=3) {
  stop("One file name and a range of columns must be supplied.n", call.=FALSE)
}
#Read input gene count table
countsTable <- read.csv(file=args[1], row.names="gene")[ ,args[2]:args[3]]
head(countsTable)
#Set control and treatment order
conds <- c(rep("ctrl",3),rep("treat",3))
#Generate list of DE genes
dgList <- DGEList(counts=countsTable, group=conds)
#Filter out counts less than 1 per million
countsPerMillion <- cpm(dgList)
keep <- which(rowSums(countsPerMillion) >= 1)
dgList <- dgList[keep,]
#View summary of filtered counts
summary(cpm(dgList))
#Normalize data
dgList <- calcNormFactors(dgList, method="TMM")
#We can examine inter-sample relationships by producing a plot
# based on mutlidimensional scaling
plotMDS(dgList)