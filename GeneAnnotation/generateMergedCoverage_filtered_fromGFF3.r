#!/usr/bin/env Rscript
#R script to generate a gene model data track from a gff3 file
#Usage: Rscript generateGeneModel_fromGFF3.r gffPath scaffoldName geneName geneStart geneEnd bamFileList

#Install Rsubread using BiocManager, this should only need to be done once
#BiocManager::install("GenomicRanges")
#Load the Rsubread library for featureCounts
library(GenomicFeatures)
library(Gviz)
library(Rsamtools)

#Retrieve input arguments
args = commandArgs(trailingOnly=TRUE)

#Turn off UCSC chromosome names to switch to custom naming
options(ucscChromosomeNames=FALSE)

#Retrieve input gff file path
inputGFF = args[1]
#Store gff3 file as TxDb object
geneDB <- makeTxDbFromGFF(file=inputGFF, organism="Daphnia pulex")

#Retrieve output path
outPath = args[2]

#Set gene region name and coordinates
#Random ex: scaffold_1 gene1 2873-4893, scaffold_1 gene10 189307-194728, scaffold_1 gene100 1175334-1176811, scaffold_2 gene1319 95-3193
#High DE ex: scaffold_9 gene6568 2352829-2354488, scaffold_97 gene18147 72973-74469, scaffold_48 gene15097 51894-56155, scaffold_2 gene1861 3631041-3632552, 
# scaffold_6 gene4926 777761-779275, scaffold_8 gene5759 504618-506995
geneLoc = args[3]
geneName = args[4]
geneStart = as.integer(args[5])
geneEnd = as.integer(args[6])

#Retrieve input bam files
inputBamUV1 = args[7]
inputBamUV2 = args[8]
inputBamUV3 = args[9]
inputBamVIS1 = args[10]
inputBamVIS2 = args[11]
inputBamVIS3 = args[12]

#Index input bam files, if they have not already been indexed
#indexBam(inputBamUV1)
#indexBam(inputBamUV2)
#indexBam(inputBamUV3)
#indexBam(inputBamVIS1)
#indexBam(inputBamVIS2)
#indexBam(inputBamVIS3)

#Create the gene region track for specified chromosome or scaffold
geneTr <- GeneRegionTrack(geneDB, 
                          chromosome = geneLoc, 
                          start = geneStart,  
                          end = geneEnd, 
                          #collapseTranscripts = "longest",
                          col="blue",
                          fill="green",
                          size=8,
                          name = geneName)

#View features in gene model
#feature(geneTr)
#group(geneTr)
#head(gene(geneTr))
#head(transcript(geneTr))
#head(exon(geneTr))
#head(symbol(geneTr))

#Create the alignment tracks
#fill = "transparent"
alTrUV1 <- AlignmentsTrack(inputBamUV1, isPaired = TRUE, name = "Coverage", type = "coverage", col = "blue", fill = "blue")
alTrUV2 <- AlignmentsTrack(inputBamUV2, isPaired = TRUE, name = "Coverage", type = "coverage", col = "blue", fill = "blue")
alTrUV3 <- AlignmentsTrack(inputBamUV3, isPaired = TRUE, name = "Coverage", type = "coverage", col = "blue", fill = "blue")
alTrVIS1 <- AlignmentsTrack(inputBamVIS1, isPaired = TRUE, name = "Coverage", type = "coverage", col = "blue", fill = "blue")
alTrVIS2 <- AlignmentsTrack(inputBamVIS2, isPaired = TRUE, name = "Coverage", type = "coverage", col = "blue", fill = "blue")
alTrVIS3 <- AlignmentsTrack(inputBamVIS3, isPaired = TRUE, name = "Coverage", type = "coverage", col = "blue", fill = "blue")

#Generate track overlay
ot <- OverlayTrack(trackList = list(alTrUV1, alTrUV2, alTrUV3, alTrVIS1, alTrVIS2, alTrVIS3))

#Generate data track of coverages
#dt <- DataTrack(alTrUV1, alTrUV2, alTrUV3, alTrVIS1, alTrVIS2, alTrVIS3)

#Plot all data tracks
jpeg(file=paste(outPath,"/filteredCoverage_",geneName,".jpeg"), width=900, height=400)
plotTracks(c(geneTr, ot), 
           from = geneStart,
           to = geneEnd, 
           chromosome = geneLoc,
           #showId = TRUE,
           exonAnnotation = "feature",
           fontcolor.exon = 1,
           fontsize.exon = 8,
           thinBoxFeature = "UTR", 
           collapse = FALSE,
           extend.left = 500, 
           extend.right = 500,
           alpha.title = 1)
           #stackedBars = TRUE,
           #col.coverage = "green",
           #fill.coverage = "blue",
           #col.mates = "purple", 
           #col.gap = "orange"
dev.off()
