#!/usr/bin/env Rscript
#R script to translate cds to protein sequences
#Usage: Rscript translateCDS_seqinr.r cdsPath

#if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
#BiocManager::install("seqinr")

#Load the libraries
library(seqinr)
library(stringr)

#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)

#Import fasta sequence
realcds <- read.fasta(file=args[1])[[1]]

#Find longest ORF
curORF=0
longORF=0
lenfORF=0
lenrORF=0
longf=0
longs="F"
for (i in 0:2) {
	#Translate to proteins for the current forward ORF
	realpep <- translate(seq=realcds, frame=i, sens="F")
	#Check length of pep sequence
	len=length(realpep)
	#Check occurence of stop codons
	numStops=str_count(realpep, "\\*")
	#Loop over each pep in current ORF
	for (j in 1:len) {
		if(numStops[j] == 1){ #Stop codon
			if(curORF > lenfORF){ #Longest current ORF
				lenfORF=curORF
			}
			curORF=0
		}else{
			curORF=curORF+1
		}
	}
	#Translate to proteins for the current reverse ORF
	realpep <- translate(seq=realcds, frame=i, sens="R")
	#Check length of pep sequence
	len=length(realpep)
	#Check occurence of stop codons
	numStops=str_count(realpep, "\\*")
	#Loop over each pep in current ORF
	for (j in len:1) {
		if(numStops[j] == 1){ #Stop codon
			if(curORF > lenrORF){ #Longest current ORF
				lenrORF=curORF
			}
			curORF=0
		}else{
			curORF=curORF+1
		}
	}
	#Keep longest ORF translation
	if(lenrORF > longORF){
		longORF=lenrORF
		longf=i
		longs="R"
	}else if(lenfORF > longORF){
		longORF=lenfORF
		longf=i
		longs="F"
	}
}

#Output longest ORF
finalpep <- translate(seq=realcds, frame=longf, sens=longs)
cat(finalpep, sep="")
