#!/usr/bin/env Rscript
#R script to perform statistical analysis of gene count tables using edgeR exact tests
# and goseq GO analysis
#Usage: Rscript geneOntology_exactTest_edgeRGoseq.r countsFile.csv startColPos endColPos
#Usage Ex: Rscript geneOntology_exactTest_edgeRGoseq.r geneCounts_merged_counted_htseqTophat2_run1_fullset_cleaned.csv 31 36

#Install necessary libraries, this should only need to be done once
#Since edgeR is already installed on the CRC this can be skipped if using the module
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("edgeR")
#BiocManager::install("goseq")
#BiocManager::install("AnnotationForge")
#Load the necessary libraries
library("edgeR")
library("goseq")
library("AnnotationForge")
#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)
#Test if there is one input argument
if (length(args)!=3) {
  stop("One file name and a range of columns must be supplied.n", call.=FALSE)
}

#Import GO annotations for Daphnia pulex from NCBI
makeOrgPackageFromNCBI(version="", outputDir="/home/mae/Downloads/", 
                       maintainer="", author="", tax_id="6669", 
                       genus="Daphnia", species="pulex")

#Import gene count data
countsTable <- read.csv(file=args[1], row.names="gene")[ ,args[2]:args[3]]
head(countsTable)
#Add grouping factor
group <- factor(c(rep("ctrl",3),rep("treat",3)))
#Create DGE list object
list <- DGEList(counts=countsTable,group=group)

#There is no purpose in analysing genes that are not expressed in either 
# experimental condition, so genes are first filtered on expression levels
keep <- filterByExpr(list)
table(keep)
list <- list[keep, , keep.lib.sizes=FALSE]
#Calculate normalized factors
list <- calcNormFactors(list)
#Write normalized counts to file
normList <- cpm(list, normalized.lib.sizes=TRUE)
write.table(normList, file="stats_normalizedCounts.csv", sep=",", row.names=TRUE)
#View normalization factors
list$samples
dim(list)

#Produce a matrix of pseudo-counts
#Estimate common dispersion and tagwise dispersions
list <- estimateDisp(list)
list$common.dispersion

#Perform an exact test for treat vs ctrl
tested <- exactTest(list, pair=c("ctrl", "treat"))
topTags(tested)
#Create results table of DE genes
resultsTbl <- topTags(tested, n=nrow(tested$table))$table
#Output resulting table
write.table(resultsTbl, file="stats_exactTest.csv", sep=",", row.names=TRUE)

#Begin goseq analysis
#First, format the DE genes into a vector suitable for use with goseq
genes=as.integer(p.adjust(tested$table$PValue[tested$table$logFC!=0],
	method="BH")<.05)
names(genes)=row.names(tested$table[tested$table$logFC!=0,])
#Output resulting table
write.table(genes, file="DEGenes.csv", sep=",", row.names=TRUE)
