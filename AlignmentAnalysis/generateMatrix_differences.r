#!/usr/bin/env Rscript
#Usage: Rscript generateMatrix_differences.r countsFile_method1.csv countsFile_method2.csv colName
#Usage Ex: Rscript generateMatrix_differences.r ../../AlignmentStats_Analysis/alignmentSummarized_legacy_subset.csv ../../AlignmentStats_Analysis/alignmentSummarized_hisat2.csv overall
#R script to generate matrix of column diferences
#Test if there are three input arguments
if (length(args)!=3) {
  stop("Two file names and one column name must be supplied.n", call.=FALSE)
}
#Retrieve gene count tables with gene IDs
aStats0 = read.csv(args[1] sep=",", row.names=1)
aStats1 = read.csv(args[2] sep=",", row.names=1)
#Create table of selected alignment stats column
dataTable <- data.frame(A=aStats0$args[3], B=aStats1$args[3])
#Generate column of differences and add to table
dataTable$C <- (dataTable$A - dataTable$B)