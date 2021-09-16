#!/usr/bin/env Rscript
#R script to translate cds to protein sequences
#Usage: Rscript translateCDS_longestORF_seqinr.r cdsPath

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

#Initialize variables
longORF=0
longf=0
longs="F"

#Find longest ORF
for (i in 0:2) {
	#Translate to proteins for the current forward ORF
	realpep <- translate(seq=realcds, frame=i, sens="F")
	#Check length of pep sequence
	len=length(realpep)
	#Check occurence of start codons
	markedpep <- sub("M", "+", realpep)
	#Loop over each pep in current ORF
	curORF=0
	startFlag=0
	stopFlag=0
	for (j in 1:len) {
		if(markedpep[j] == "+"){ #Start codon
		  curORF=curORF+1
		  startFlag=1
		}
	  if (startFlag == 1 && stopFlag == 0){ #Start codon was found
			if(markedpep[j] == "*"){ #Stop codon
			  stopFlag=1
			  if(curORF > longORF){ #Longest current ORF
			    longORF=curORF
			    longf=i
			    longs="F"
			  }
			  curORF=0
			}else{
			  curORF=curORF+1
			}
	  }
	}
	#Translate to proteins for the current reverse ORF
	#realpep <- translate(seq=realcds, frame=i, sens="R")
	#Check length of pep sequence
	#len=length(realpep)
	#Check occurence of stop codons
	#markedpep <- sub("M", "+", realpep)
	#Loop over each pep in current ORF
	#curORF=0
	#startFlag=0
	#stopFlag=0
	#for (j in 1:len) {
	#  if(markedpep[j] == "+"){ #Start codon
	#    curORF=curORF+1
	#    startFlag=1
	#  }
	#  if (startFlag == 1  && stopFlag == 0){ #Start codon was found
	#    if(markedpep[j] == "*"){ #Stop codon
	#      stopFlag=1
	#      if(curORF > longORF){ #Longest current ORF
	#        longORF=curORF
	#        longf=i
	#        longs="R"
	#      }
	#      curORF=0
	#    }else{
	#      curORF=curORF+1
	#    }
	#  }
	#}
}

#Output longest ORF
if (longORF == 0) {
  print("NoValidORFFound")
}else{
  finalpep <- translate(seq=realcds, frame=longf, sens=longs)
  cat(finalpep, sep="")
}

