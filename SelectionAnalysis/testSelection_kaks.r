#if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
#BiocManager::install("seqinr")

#Load the libraries
library(seqinr)

#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)

#Import file with protein alignments
s <- read.alignment(file=args[1], format = "fasta")

#Generate ka ks values
kaks(s)