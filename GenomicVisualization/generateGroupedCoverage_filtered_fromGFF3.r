#!/usr/bin/env Rscript
#R script to generate a gene model data track from a gff3 file
#Usage: Rscript generateGeneModel_fromGFF3.r gffPath scaffoldName geneName geneStart geneEnd bamFile

#Install Rsubread using BiocManager, this should only need to be done once
#BiocManager::install("GenomicRanges")
#Load the Rsubread library for featureCounts
library(GenomicFeatures)
library(Gviz)
library(Rsamtools)
library(BSgenome)
library(BSgenome.Dpulex.KAP4)
library(rcartocolor)

#Turn off UCSC chromosome names to switch to custom naming
options(ucscChromosomeNames=FALSE)

#Retrieve input arguments
#args = commandArgs(trailingOnly=TRUE)

# set working directory
setwd("/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/ncbi_dataset/data/GCF_021134715.1/seqs")

# Plotting Palettes
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
# https://github.com/Nowosad/rcartocolor
plotColors <- carto_pal(12, "Safe")

# forge the BSgenome data package
#forgeBSgenomeDataPkg("/Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/BSgenome.Dpulex.KAP4-seed", replace=TRUE)

# set gene region name and coordinates
# photolyase (LOC124188748) NC_060018.1:363588-367123
geneLoc = "NC_060018.1"
geneName = "LOC124188748"
geneStart = 363588
geneEnd = 367123

# retrieve chromosome sequence
#sTrack <- SequenceTrack(Dpulex)
sTrack <- SequenceTrack(Dpulex, chromosome = geneLoc)

# view gene sequence
#plotTracks(sTrack, chromosome = geneLoc, from = geneStart, to = geneEnd)
#plotTracks(sTrack, from = geneStart, to = geneEnd)

#Retrieve input gff file path
inputGFF = "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Bioinformatics/variantsCalled_samtoolsBcftools/geneModels/LOC124188748.gff"
#Store gff3 file as TxDb object
geneDB <- makeTxDbFromGFF(file=inputGFF, format="gff3", organism="Daphnia pulex")

#Retrieve output path
outPath = "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Bioinformatics/variantsCalled_samtoolsBcftools/geneModels"

#Retrieve input bam files
inputBamE05 = "/Volumes/My Passport/OLYM_dMelUV/alignments/E05/filteredMapQ_readGroups.bam"
inputBamR2 = "/Volumes/My Passport/OLYM_dMelUV/alignments/R2/filteredMapQ_readGroups.bam"
inputBamY05 = "/Volumes/My Passport/OLYM_dMelUV/alignments/Y05/filteredMapQ_readGroups.bam"
inputBamY023 = "/Volumes/My Passport/OLYM_dMelUV/alignments/Y023/filteredMapQ_readGroups.bam"

#Index input bam files, if they have not already been indexed
#indexBam(inputBamE05)
#indexBam(inputBamR2)
#indexBam(inputBamY05)
#indexBam(inputBamY023)

#Create the gene region track for specified chromosome or scaffold
geneTr <- GeneRegionTrack(geneDB, 
                          chromosome = geneLoc, 
                          start = geneStart,  
                          end = geneEnd, 
                          #collapseTranscripts = "longest",
                          col=plotColors[4],
                          fill=plotColors[7],
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
alTrE05 <- AlignmentsTrack(inputBamE05, isPaired = TRUE, name = "E05", type = "coverage", col = plotColors[5], fill = plotColors[11])
alTrR2 <- AlignmentsTrack(inputBamR2, isPaired = TRUE, name = "R2", type = "coverage", col = plotColors[5], fill = plotColors[11])
alTrY05 <- AlignmentsTrack(inputBamY05, isPaired = TRUE, name = "Y05", type = "coverage", col = plotColors[5], fill = plotColors[11])
alTrY023 <- AlignmentsTrack(inputBamY023, isPaired = TRUE, name = "Y023", type = "coverage", col = plotColors[5], fill = plotColors[11])

#Generate track overlay
#ot <- OverlayTrack(trackList = list(alTrE05, alTrR2, alTrY05, alTrY023))

#Generate data track of coverages
#dt <- DataTrack(alTrUV1, alTrUV2, alTrUV3, alTrVIS1, alTrVIS2, alTrVIS3)

# plot all data tracks
jpeg(file=paste(outPath,"/filteredCoverage_geneModel_",geneName,".jpeg", sep = ""), width=900, height=400)
plotTracks(c(geneTr, sTrack, alTrE05, alTrR2, alTrY05, alTrY023), 
           #groups = c("E05", "R2", "Y05", "Y023"),
           from = geneStart,
           to = geneEnd, 
           chromosome = geneLoc,
           #showId = TRUE,
           exonAnnotation = "feature",
           fontcolor.exon = 1,
           fontsize.exon = 8,
           thinBoxFeature = "UTR", 
           collapse = FALSE,
           #extend.left = 500, 
           #extend.right = 500,
           alpha.title = 1,
           showIndels=TRUE)
           #stackedBars = TRUE,
           #col.coverage = "green",
           #fill.coverage = "blue",
           #col.mates = "purple", 
           #col.gap = "orange"
dev.off()

# set coordinates for cds 1
cdsStart = 364812
cdsEnd = 365085

# plot cds 1
plotTracks(c(geneTr, sTrack, alTrE05, alTrR2, alTrY05, alTrY023), 
           from = cdsStart,
           to = cdsEnd, 
           chromosome = geneLoc,
           exonAnnotation = "feature",
           fontcolor.exon = 1,
           fontsize.exon = 8,
           thinBoxFeature = "UTR", 
           collapse = FALSE,
           alpha.title = 1,
           showIndels=TRUE)

# set coordinates for cds 2
cdsStart = 365154
cdsEnd = 365279

# plot cds 2
plotTracks(c(geneTr, sTrack, alTrE05, alTrR2, alTrY05, alTrY023), 
           from = cdsStart,
           to = cdsEnd, 
           chromosome = geneLoc,
           exonAnnotation = "feature",
           fontcolor.exon = 1,
           fontsize.exon = 8,
           thinBoxFeature = "UTR", 
           collapse = FALSE,
           alpha.title = 1,
           showIndels=TRUE)

# set coordinates for cds 3
cdsStart = 365353
cdsEnd = 365555

# plot cds 3
plotTracks(c(geneTr, sTrack, alTrE05, alTrR2, alTrY05, alTrY023), 
           from = cdsStart,
           to = cdsEnd, 
           chromosome = geneLoc,
           exonAnnotation = "feature",
           fontcolor.exon = 1,
           fontsize.exon = 8,
           thinBoxFeature = "UTR", 
           collapse = FALSE,
           alpha.title = 1,
           showIndels=TRUE)

# set coordinates for cds 4
cdsStart = 365632
cdsEnd = 365686

# plot cds 4
plotTracks(c(geneTr, sTrack, alTrE05, alTrR2, alTrY05, alTrY023), 
           from = cdsStart,
           to = cdsEnd, 
           chromosome = geneLoc,
           exonAnnotation = "feature",
           fontcolor.exon = 1,
           fontsize.exon = 8,
           thinBoxFeature = "UTR", 
           collapse = FALSE,
           alpha.title = 1,
           showIndels=TRUE)

# set coordinates for cds 5
cdsStart = 365756
cdsEnd = 365859

# plot cds 5
plotTracks(c(geneTr, sTrack, alTrE05, alTrR2, alTrY05, alTrY023), 
           from = cdsStart,
           to = cdsEnd, 
           chromosome = geneLoc,
           exonAnnotation = "feature",
           fontcolor.exon = 1,
           fontsize.exon = 8,
           thinBoxFeature = "UTR", 
           collapse = FALSE,
           alpha.title = 1,
           showIndels=TRUE)

# set coordinates for cds 6
cdsStart = 365930
cdsEnd = 366128

# plot cds 6
plotTracks(c(geneTr, sTrack, alTrE05, alTrR2, alTrY05, alTrY023), 
           from = cdsStart,
           to = cdsEnd, 
           chromosome = geneLoc,
           exonAnnotation = "feature",
           fontcolor.exon = 1,
           fontsize.exon = 8,
           thinBoxFeature = "UTR", 
           collapse = FALSE,
           alpha.title = 1,
           showIndels=TRUE)

# set coordinates for cds 7
cdsStart = 366201
cdsEnd = 366353

# plot cds 7
plotTracks(c(geneTr, sTrack, alTrE05, alTrR2, alTrY05, alTrY023), 
           from = cdsStart,
           to = cdsEnd, 
           chromosome = geneLoc,
           exonAnnotation = "feature",
           fontcolor.exon = 1,
           fontsize.exon = 8,
           thinBoxFeature = "UTR", 
           collapse = FALSE,
           alpha.title = 1,
           showIndels=TRUE)

# set coordinates for cds 8
cdsStart = 366418
cdsEnd = 366457

# plot cds 8
plotTracks(c(geneTr, sTrack, alTrE05, alTrR2, alTrY05, alTrY023), 
           from = cdsStart,
           to = cdsEnd, 
           chromosome = geneLoc,
           exonAnnotation = "feature",
           fontcolor.exon = 1,
           fontsize.exon = 8,
           thinBoxFeature = "UTR", 
           collapse = FALSE,
           alpha.title = 1,
           showIndels=TRUE)

# set coordinates for cds 9
cdsStart = 366520
cdsEnd = 366642

# plot cds 9
plotTracks(c(geneTr, sTrack, alTrE05, alTrR2, alTrY05, alTrY023), 
           from = cdsStart,
           to = cdsEnd, 
           chromosome = geneLoc,
           exonAnnotation = "feature",
           fontcolor.exon = 1,
           fontsize.exon = 8,
           thinBoxFeature = "UTR", 
           collapse = FALSE,
           alpha.title = 1,
           showIndels=TRUE)

# set coordinates for cds 10
cdsStart = 366711
cdsEnd = 367002

# plot cds 10
plotTracks(c(geneTr, sTrack, alTrE05, alTrR2, alTrY05, alTrY023), 
           from = cdsStart,
           to = cdsEnd, 
           chromosome = geneLoc,
           exonAnnotation = "feature",
           fontcolor.exon = 1,
           fontsize.exon = 8,
           thinBoxFeature = "UTR", 
           collapse = FALSE,
           alpha.title = 1,
           showIndels=TRUE)
