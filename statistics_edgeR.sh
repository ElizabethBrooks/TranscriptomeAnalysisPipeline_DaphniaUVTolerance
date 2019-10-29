#!/usr/bin/env Rscript
#R script to perform statistical analysis of gene count tables using edgeR
#Install edgeR, this should only need to be done once
#Since edgeR is already installed on the CRC this can be skipped if using the module
#bioLite("edgeR")
#Load the edgeR library
library("edgeR")
#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)
#Test if there is one input argument
if (length(args)!=1) {
  stop("One file name must be supplied.n", call.=FALSE)
}
#Read input gene count table
countsTable <- read.delim( file=arg[1], row.names="gene" ) head(countsTable)
#Set control and treatment order
conds <- factor( c("ctrl","ctrl","ctrl","treat","treat","treat"  ) )
#Generate list of DE genes
cds<- DGEList( counts=countsTable, group=conds )
d <- calcNormFactors( cds )
d <- estimateCommonDisp( d )
d <- estimateTagwiseDisp( d)
de<- exactTest( d , pair = c( "ctrl" , "treat" ) )
#Create results table of DE genes
resultsTbl <- topTags( de, n = nrow( de$table ) )$table
#Output resulting table
write.table(resultsTbl, file=arg[1]".out.csv", sep=",", row.names=TRUE)