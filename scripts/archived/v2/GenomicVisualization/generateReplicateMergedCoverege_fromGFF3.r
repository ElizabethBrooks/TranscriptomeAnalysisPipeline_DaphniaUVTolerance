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
#Set output merged bam file paths
destUV = "/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/mergedGeneRegion_UV.bam"
destVIS = "/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/mergedGeneRegion_VIS.bam"

#Index input bam files, if they have not already been indexed
#indexBam(inputBamUV1)
#indexBam(inputBamUV2)
#indexBam(inputBamUV3)
#indexBam(inputBamVIS1)
#indexBam(inputBamVIS2)
#indexBam(inputBamVIS3)

#Set gene region name and coordinates
#Random ex: scaffold_1 dp_gene1 2873-4893, scaffold_1 dp_gene10 189307-194728, scaffold_1 dp_gene100 1175334-1176811, scaffold_2 dp_gene1319 95-3193
#High DE ex: scaffold_9 dp_gene6568 2352829-2354488
geneLoc = "scaffold_9"
geneName = "dp_gene6588"
extBy = 100
geneStart = 2352829 - extBy
geneEnd = 2354488 + extBy

#Retrieve gene features
transReg <- transcriptsBy(geneDB, by="gene")
transRange <- transReg[[geneName]]

#Merge UV bam files
mergeBam(files=c(inputBamUV1, inputBamUV2, inputBamUV3), region=transRange, destination=destUV, overwrite=TRUE)
indexBam(destUV)

#Merge VIS bam files
mergeBam(files=c(inputBamVIS1, inputBamVIS2, inputBamVIS3), region=transRange, destination=destVIS, overwrite=TRUE)
indexBam(destVIS)

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

#Create the read coverage tracks
#covTrUV <- DataTrack(destUV, isPaired = TRUE, type = "heatmap", name = "CoverageUV")
#covTrVIS <- DataTrack(destVIS, isPaired = TRUE, type = "heatmap", name = "CoverageVIS")
#hTrUV <- DataTrack(destUV, isPaired = TRUE, type = "h", name = "CoverageUV")
#hTrVIS <- DataTrack(destVIS, isPaired = TRUE, type = "h", name = "CoverageVIS")

#Create the alignment tracks
alTrUV <- AlignmentsTrack(destUV, isPaired = TRUE, name = "AlignmentsUV", type = "coverage")
alTrVIS <- AlignmentsTrack(destVIS, isPaired = TRUE, name = "AlignmentsVIS", type = "coverage")

#Generate track overlay
#ot <- OverlayTrack(trackList = list(alTrUV, alTrVIS))

#Plot all data tracks
plotTracks(c(geneTr, alTrUV, alTrVIS), 
           from = geneStart,
           to = geneEnd, 
           chromosome = geneLoc, 
           exonAnnotation = "feature",
           fontcolor.exon = 1,
           fontsize.exon = 10,
           thinBoxFeature = "UTR", 
           collapse = FALSE,
           col.coverage = "green",
           fill.coverage = "blue")
           #extend.left = 100, 
           #extend.right = 100,
           #col.mates = "purple", 
           #col.gap = "orange",