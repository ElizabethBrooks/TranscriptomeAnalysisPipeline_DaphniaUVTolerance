#!/usr/bin/env Rscript
#R script to translate cds to protein sequences
#Usage: Rscript translateCDS_seqinr.r cdsPath

#if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
#BiocManager::install("seqinr")

#Load the libraries
library(seqinr)

#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)

realcds <- read.fasta(file=args[1])[[1]]
translate(seq=realcds)