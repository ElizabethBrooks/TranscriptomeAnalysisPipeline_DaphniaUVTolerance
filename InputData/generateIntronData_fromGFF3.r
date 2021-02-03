#!/usr/bin/env Rscript
#R script to generate intron data from a gff3 file
#Usage: Rscript generateIntronData_fromGFF3.r gffPath
#Usage Ex: Rscript generateIntronData_fromGFF3.r /afs/crc.nd.edu/group/pfrenderlab/bateson/ebrooks/rnaseq/genomicResources_PA42_v4.1/PA42.4.1.gff
#Usage Ex: Rscript generateIntronData_fromGFF3.r /afs/crc.nd.edu/group/pfrenderlab/bateson/ebrooks/rnaseq/genomicResources_PA42_v3.0/PA42.3.0.annotation.18440.gff

#Install Rsubread using BiocManager, this should only need to be done once
#BiocManager::install("GenomicFeatures")
#Load the Rsubread library for featureCounts
library(GenomicFeatures)
#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)
#Test if there is one input argument
if (length(args)!=2) {
  stop("One genome file path and an aligned/sorted read folder path must be supplied.n", 
  	call.=FALSE)
}
#Store gff3 file as TxDb object
geneDB <- makeTxDbFromGFF(file=args[1])
#Retrieve intron data
intronData <- intronsByTranscript(geneDB)
#View intron data
intronData <- unlist(intronData)
summary(width(intronData))