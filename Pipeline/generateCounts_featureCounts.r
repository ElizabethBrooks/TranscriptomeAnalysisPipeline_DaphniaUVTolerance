#!/usr/bin/env Rscript
#Usage: Rscript generateCounts_featureCounts.r
#Usage Ex: Rscript generateCounts_featureCounts.r
#R script to perform analysis of aliged paired-end reads and generate gene count tables
#Install Rsubread using BiocManager, this should only need to be done once
#install.packages("BiocManager")
#BiocManager::install("Rsubread")
#Load the Rsubread library for featureCounts
library(Rsubread)

#Retrieve all sorted bam files in directory and sub directories
bam.files <- list.files(path="/afs/crc.nd.edu/group/pfrenderlab/bateson/ebrooks/rnaseq/aligned_hisat2_run2", pattern="accepted_hits.sam", full.names=TRUE, recursive=TRUE)

#The mapped reads can be counted across genes by using the featureCounts function
# and input gtf/gff file, paired end reads, and 8 threads
fc <- featureCounts(bam.files, annot.ext="/afs/crc.nd.edu/group/pfrenderlab/bateson/ebrooks/rnaseq/PA42.3.0.annotation.18440.gff", isGTFAnnotationFile=TRUE, GTF.attrType="ID", isPairedEnd=TRUE, autosort=TRUE, nthreads=8)

# See what slots are stored in fc
names(fc)

## Take a look at the featurecounts stats
fc$stat

## Take a look at the dimensions to see the number of genes for samples
dim(fc$counts)
## Take a look at the first 6 lines
head(fc$counts)

#The annotation slot shows the annotation information that featureCounts used to summarise reads over genes
head(fc$annotation)