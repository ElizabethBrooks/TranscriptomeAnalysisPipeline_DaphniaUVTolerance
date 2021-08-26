#!/usr/bin/env Rscript
#R script to generate a gene model data track from a gff3 file
#Usage: Rscript generateGeneModel_fromGFF3.r gffPath

#Install Rsubread using BiocManager, this should only need to be done once
#BiocManager::install("GenomicRanges")
#Load the Rsubread library for featureCounts
library(GenomicFeatures)
library(Gviz)
library(Rsamtools)

#Turn off UCSC chromosome names to switch to custom naming
options(ucscChromosomeNames=FALSE)

#Retrieve input gff file path
inputGFF = "/home/mae/Documents/RNASeq_Workshop_ND/dp_gene15097.gff"
#Store gff3 file as TxDb object
geneDB <- makeTxDbFromGFF(file=inputGFF, organism="Daphnia melanica")

#Retrieve input bam files
inputBamUV1 = "/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/140327_I481_FCC3P1PACXX_L2_Pool_1_PA_UV/filteredMapQ.bam"
inputBamUV2 = "/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/140327_I481_FCC3P1PACXX_L3_Pool_2_PA_UV/filteredMapQ.bam"
inputBamUV3 = "/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/140327_I481_FCC3P1PACXX_L4_Pool_3_PA_UV/filteredMapQ.bam"
inputBamVIS1 = "/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/140327_I481_FCC3P1PACXX_L2_Pool_1_PA_VIS/filteredMapQ.bam"
inputBamVIS2 = "/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/140327_I481_FCC3P1PACXX_L3_Pool_2_PA_VIS/filteredMapQ.bam"
inputBamVIS3 = "/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/140327_I481_FCC3P1PACXX_L4_Pool_3_PA_VIS/filteredMapQ.bam"

#Index input bam files, if they have not already been indexed
#indexBam(inputBamUV1)
#indexBam(inputBamUV2)
#indexBam(inputBamUV3)
#indexBam(inputBamVIS1)
#indexBam(inputBamVIS2)
#indexBam(inputBamVIS3)

#Set gene region name and coordinates
#Random ex: scaffold_1 gene1 2873-4893, scaffold_1 gene10 189307-194728, scaffold_1 gene100 1175334-1176811, scaffold_2 gene1319 95-3193
#High DE ex: scaffold_9 gene6568 2352829-2354488, scaffold_97 gene18147 72973-74469, scaffold_48 gene15097 51894-56155, scaffold_2 gene1861 3631041-3632552, 
# scaffold_6 gene4926 777761-779275, scaffold_8 gene5759 504618-506995
geneLoc = "scaffold_48"
geneName = "gene15097"
geneStart = 51894
geneEnd = 56155

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
alTrUV1 <- AlignmentsTrack(inputBamUV1, isPaired = TRUE, name = "Coverage", type = "coverage", col = "red", fill = "red", alpha = 0.5)
alTrUV2 <- AlignmentsTrack(inputBamUV2, isPaired = TRUE, name = "Coverage", type = "coverage", col = "red", fill = "red", alpha = 0.5)
alTrUV3 <- AlignmentsTrack(inputBamUV3, isPaired = TRUE, name = "Coverage", type = "coverage", col = "red", fill = "red", alpha = 0.5)
alTrVIS1 <- AlignmentsTrack(inputBamVIS1, isPaired = TRUE, name = "Coverage", type = "coverage", col = "blue", fill = "blue", alpha = 0.5)
alTrVIS2 <- AlignmentsTrack(inputBamVIS2, isPaired = TRUE, name = "Coverage", type = "coverage", col = "blue", fill = "blue", alpha = 0.5)
alTrVIS3 <- AlignmentsTrack(inputBamVIS3, isPaired = TRUE, name = "Coverage", type = "coverage", col = "blue", fill = "blue", alpha = 0.5)

#Generate track overlay
ot <- OverlayTrack(trackList = list(alTrUV1, alTrUV2, alTrUV3, alTrVIS1, alTrVIS2, alTrVIS3))

#Generate data track of coverages
#dt <- DataTrack(alTrUV1, alTrUV2, alTrUV3, alTrVIS1, alTrVIS2, alTrVIS3)

#Plot all data tracks
jpeg(file=paste("/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/filteredIndCoverage_",geneName,".jpeg"), width=900, height=400)
plotTracks(c(geneTr, ot), 
           from = geneStart,
           to = geneEnd, 
           chromosome = geneLoc,
           showId = TRUE,
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
