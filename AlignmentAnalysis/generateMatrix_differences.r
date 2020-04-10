#!/usr/bin/env Rscript
#Usage: Rscript generateMatrix_differences.r percentsFile_legacy.csv percentsFile_newMethod.csv
#Usage Ex: Rscript generateMatrix_differences.r ../../AlignmentStats_Analysis/alignmentSummarized_legacy_subset_trimmed.csv ../../AlignmentStats_Analysis/alignmentSummarized_hisat2_trimmed.csv
#R script to generate matrix of column diferences

#Retrieve input file name of gene counts
args=commandArgs(trailingOnly=TRUE)
#Test if there are three input arguments
if (length(args)!=2) {
  stop("Two file names and one column name must be supplied.n", call.=FALSE)
}
#Retrieve overall alignment percents columns
aStatsOA = read.csv(args[1], sep=",")[ , 1]
aStatsOB = read.csv(args[2], sep=",")[ , 1]
#Retrieve concordant alignment percents columns
aStatsCA = read.csv(args[1], sep=",")[ , 2]
aStatsCB = read.csv(args[2], sep=",")[ , 2]
#Create table of overall alignment stats column
dataTableO <- data.frame(OA=aStatsOA, OB=aStatsOB)
#Create table of concordant alignment stats column
dataTableC <- data.frame(CA=aStatsCA, CB=aStatsCB)
#Generate column of overall squared differences and add to table
dataTableOD <- (dataTableO$OA - dataTableO$OB)*(dataTableO$OA - dataTableO$OB)
#Generate column of concordant squared differences and add to table
dataTableCD <- (dataTableC$CA - dataTableC$CB)*(dataTableC$CA - dataTableC$CB)
#Write generated differences matrices to csv files
write.csv(dataTableOD, file="alignmentSummarized_differences_overall.csv", row.names=FALSE)
write.csv(dataTableCD, file="alignmentSummarized_differences_concordant.csv", row.names=FALSE)