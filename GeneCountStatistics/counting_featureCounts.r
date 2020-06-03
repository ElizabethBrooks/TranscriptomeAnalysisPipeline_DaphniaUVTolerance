#!/usr/bin/env Rscript
#Usage: Rscript counting_featureCounts.r readPath genomeFeaturesFile
#Usage Ex: Rscript counting_featureCounts.r /afs/crc.nd.edu/group/pfrenderlab/mendel/ebrooks/rnaseq/aligned_hisat2_run1 /afs/crc.nd.edu/group/pfrenderlab/mendel/ebrooks/rnaseq/PA42.3.0.annotation.18440.gff

#R script to perform analysis of aliged paired-end reads and generate gene count tables
#Install Rsubread using BiocManager, this should only need to be done once
#install.packages("BiocManager")
#BiocManager::install("Rsubread")
#Load the Rsubread library for featureCounts
library(Rsubread)
#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)
#Test if there is one input argument
if (length(args)!=2) {
  stop("One genome file path and an aligned/sorted read folder path must be supplied.n", 
  	call.=FALSE)
}
#Retrieve all sorted bam files in sub directories
bam.files <- list.files(path=args[1], pattern="accepted_hits.bam", full.names=TRUE, 
	recursive=TRUE)
#The mapped reads can be counted across genes by using the featureCounts function
# and input gtf/gff file, paired end reads, and 8 threads
fc <- featureCounts(bam.files, annot.ext=args[2], isGTFAnnotationFile=TRUE,
	GTF.attrType="ID", isPairedEnd=TRUE, nthreads=8)
#See what slots are stored in fc
names(fc)
#Take a look at the featurecounts stats
fc$stat
#Take a look at the dimensions to see the number of genes for samples
dim(fc$counts)
#Take a look at the first 6 lines
head(fc$counts)
#The annotation slot shows the annotation information that featureCounts used to summarise reads over genes
head(fc$annotation)