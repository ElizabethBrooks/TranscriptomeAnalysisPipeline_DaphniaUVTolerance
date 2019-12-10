#!/usr/bin/env Rscript
#Usage: Rscript generateCounts_featureCounts.r
#Usage Ex: Rscript generateCounts_featureCounts.r
#R script to perform analysis of aliged paired-end reads and generate gene count tables
#Install Rsubread using BiocManager, this should only need to be done once
#install.packages("BiocManager")
#BiocManager::install("Rsubread")
#Load the Rsubread library for featureCounts
library(Rsubread)

bam.files <- list.files(path = "./data", pattern = ".BAM$", full.names = TRUE)
bam.files

#The mapped reads can be counted across mouse genes by using the featureCounts function
fc <- featureCounts(bam.files, annot.inbuilt="mm10")

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