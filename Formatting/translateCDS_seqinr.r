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

#Import fasta sequence
realcds <- read.fasta(file=args[1])[[1]]

#Find longest ORF
curORF=0
longORF=0
lenORF=0
longf=0
longs="F"
for (i in 1:6) {
	#Check frame
	if(i == 1 | 4){
		fNum=0
	}elseif(i == 2 | 5){
		fNum=1
	}else{
		fNum=2
	}
	#Check sense
	if(i < 4){
		sTag="F"
	}else{
		sTag="R"
	}
	#Translate to proteins for the current ORF
	realpep <- translate(seq=realcds, frame=fNum, sens=sTag)
	#Check length of pep sequence
	len=length(realpep)
	#Loop over each pep in current ORF
	for (i in 1:len) {
		if(realpep[i] == "\\*"){
			curORF=curORF+1
		}else{
			if(curORF > lenORF){ 
				lenORF=curORF
			}
			curORF=0
		}
	}
	#Keep longest ORF translation
	if(lenORF > longORF){
		longORF=lenORF
		longf=fNum
		longs=sTag
	}
}

#Output longest ORF
finalpep <- translate(seq=realcds, frame=longf, sens=longs)
cat(finalpep, sep="")
