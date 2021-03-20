#if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
#BiocManager::install("seqinr")

#Load the libraries
#library(seqinr)

#Import file with protein alignments
s <- read.alignment(file = "/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/gene4487_cds_allDaphnia_aligned.fasta", format = "fasta")

#Generate ka ks values
kaks(s)
