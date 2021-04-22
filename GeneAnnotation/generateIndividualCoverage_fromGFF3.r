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
inputGFF = "/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/PA42.4.1.gff"
#Store gff3 file as TxDb object
geneDB <- makeTxDbFromGFF(file=inputGFF, organism="Daphnia melanica")

#Retrieve input bam files
inputBamUV1 = "/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/140327_I481_FCC3P1PACXX_L2_Pool_1_PA_UV/accepted_hits.bam"
inputBamUV2 = "/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/140327_I481_FCC3P1PACXX_L3_Pool_2_PA_UV/accepted_hits.bam"
inputBamUV3 = "/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/140327_I481_FCC3P1PACXX_L4_Pool_3_PA_UV/accepted_hits.bam"
inputBamVIS1 = "/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/140327_I481_FCC3P1PACXX_L2_Pool_1_PA_VIS/accepted_hits.bam"
inputBamVIS2 = "/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/140327_I481_FCC3P1PACXX_L3_Pool_2_PA_VIS/accepted_hits.bam"
inputBamVIS3 = "/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/140327_I481_FCC3P1PACXX_L4_Pool_3_PA_VIS/accepted_hits.bam"
#Output filtered bam file paths
destUV1 = "/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/filtered_hits_UV1.bam"
destUV2 = "/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/filtered_hits_UV2.bam"
destUV3 = "/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/filtered_hits_UV3.bam"
destVIS1 = "/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/filtered_hits_VIS1.bam"
destVIS2 = "/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/filtered_hits_VIS2.bam"
destVIS3 = "/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/filtered_hits_VIS3.bam"

#Index input bam files, if they have not already been indexed
#indexBam(inputBamUV1)
#indexBam(inputBamUV2)
#indexBam(inputBamUV3)
#indexBam(inputBamVIS1)
#indexBam(inputBamVIS2)
#indexBam(inputBamVIS3)

#Create bam file mapq filter
sbp <- ScanBamParam(what=scanBamWhat(), mapqFilter=40)

#Filter input bam files
scanBam(inputBamUV1, destUV1, param=sbp)
scanBam(inputBamUV2, destUV2, param=sbp)
scanBam(inputBamUV3, destUV3, param=sbp)
scanBam(inputBamVIS1, destVIS1, param=sbp)
scanBam(inputBamVIS2, destVIS2, param=sbp)
scanBam(inputBamVIS3, destVIS3, param=sbp)

#Index filtered bam files, if they have not already been indexed
indexBam(inputBamUV1)
indexBam(inputBamUV2)
indexBam(inputBamUV3)
indexBam(inputBamVIS1)
indexBam(inputBamVIS2)
indexBam(inputBamVIS3)

#Set gene region name and coordinates
#Random ex: scaffold_1 dp_gene1 2873-4893, scaffold_1 dp_gene10 189307-194728, scaffold_1 dp_gene100 1175334-1176811, scaffold_2 dp_gene1319 95-3193
#High DE ex: scaffold_9 dp_gene6568 2352829-2354488, scaffold_97 dp_gene18147 72973-74469, scaffold_48 dp_gene15097 51894-56155, scaffold_2 dp_gene1861 3631041-3632552, 
# scaffold_6 dp_gene4926 777761-779275, scaffold_8 dp_gene5759 504618-506995
geneLoc = "scaffold_8"
geneName = "dp_gene5759"
geneStart = 504618
geneEnd = 506995

#Create the gene region track for specified chromosome or scaffold
geneTr <- GeneRegionTrack(geneDB, 
                          chromosome = geneLoc, 
                          start = geneStart,  
                          end = geneEnd, 
                          col="blue",
                          fill="green",
                          size=10,
                          name = geneName)

#View features in gene model
#feature(geneTr)
#head(gene(geneTr))
#head(transcript(geneTr))
#head(exon(geneTr))
#head(symbol(geneTr))

#Create the alignment tracks
alTrUV1 <- AlignmentsTrack(gaUV1, isPaired = TRUE, name = "Alignments", type = "coverage", fill = "red")
alTrUV2 <- AlignmentsTrack(gaUV2, isPaired = TRUE, name = "Alignments", type = "coverage", fill = "yellow")
alTrUV3 <- AlignmentsTrack(gaUV3, isPaired = TRUE, name = "Alignments", type = "coverage", fill = "orange")
alTrVIS1 <- AlignmentsTrack(gaVIS1, isPaired = TRUE, name = "Alignments", type = "coverage", fill = "purple")
alTrVIS2 <- AlignmentsTrack(gaVIS2, isPaired = TRUE, name = "Alignments", type = "coverage", fill = "pink")
alTrVIS3 <- AlignmentsTrack(gaVIS3, isPaired = TRUE, name = "Alignments", type = "coverage", fill = "blue")

#Generate track overlay
ot <- OverlayTrack(trackList = list(alTrUV1, alTrUV2, alTrUV3, alTrVIS1, alTrVIS2, alTrVIS3))

#Plot all data tracks
jpeg(file=paste("/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/eachCoverage_",geneName,".jpeg"), width=900, height=400)
plotTracks(c(geneTr, ot), 
           from = geneStart,
           to = geneEnd, 
           chromosome = geneLoc, 
           exonAnnotation = "feature",
           fontcolor.exon = 1,
           fontsize.exon = 10,
           thinBoxFeature = "UTR", 
           collapse = FALSE,
           extend.left = 100, 
           extend.right = 100)
           #col.coverage = "green",
           #fill.coverage = "blue",
           #col.mates = "purple", 
           #col.gap = "orange"
dev.off()
