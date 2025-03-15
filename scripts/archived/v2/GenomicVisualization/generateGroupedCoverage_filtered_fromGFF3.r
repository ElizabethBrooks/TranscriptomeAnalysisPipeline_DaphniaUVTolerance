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

# set output path
outPath = "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Bioinformatics/variantsCalled_samtoolsBcftools/geneModels"

# Plotting Palettes
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
# https://github.com/Nowosad/rcartocolor
plotColors <- carto_pal(12, "Safe")

# forge the BSgenome data package
#forgeBSgenomeDataPkg("/Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/BSgenome.Dpulex.KAP4-seed", replace=TRUE)


# retrieve chromosome sequence
#sTrack <- SequenceTrack(Dpulex)
sTrack <- SequenceTrack(Dpulex, chromosome = geneLoc)

# view gene sequence
#plotTracks(sTrack, chromosome = geneLoc, from = geneStart, to = geneEnd)
#plotTracks(sTrack, from = geneStart, to = geneEnd)

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

#Create the alignment tracks
#fill = "transparent"
alTrE05 <- AlignmentsTrack(inputBamE05, isPaired = TRUE, name = "E05", type = "coverage", col = plotColors[5], fill = plotColors[11])
alTrR2 <- AlignmentsTrack(inputBamR2, isPaired = TRUE, name = "R2", type = "coverage", col = plotColors[5], fill = plotColors[11])
alTrY05 <- AlignmentsTrack(inputBamY05, isPaired = TRUE, name = "Y05", type = "coverage", col = plotColors[5], fill = plotColors[11])
alTrY023 <- AlignmentsTrack(inputBamY023, isPaired = TRUE, name = "Y023", type = "coverage", col = plotColors[5], fill = plotColors[11])

# set gene region name and coordinates
# photolyase (LOC124188748) NC_060018.1:363588-367123
geneLoc = "NC_060018.1"
geneName = "LOC124188748"
geneStart = 363588
geneEnd = 367123
# NFKB (LOC124190117) NC_060019.1:3336006-3343487
#geneLoc = "NC_060019.1"
#geneName = "LOC124190117"
#geneStart = 3336006
#geneEnd = 3343487
# adrenodoxin (LOC124193197) NC_060020.1:3488210-3489206
#geneLoc = "NC_060020.1"
#geneName = "LOC124193197"
#geneStart = 3488210
#geneEnd = 3489206
# mucin2 (LOC124202665) NC_060025.1:8795802-8799903
#geneLoc = "NC_060025.1"
#geneName = "LOC124202665"
#geneStart = 8795802
#geneEnd = 8799903

## Important
## geneName="XXXX"
## cat /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/ncbi_dataset/data/GCF_021134715.1/genomic.gff | grep -w "$geneName" > /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Bioinformatics/variantsCalled_samtoolsBcftools/geneModels/$geneName.gff
#Retrieve input gene gff file path
inputGFF = paste("/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Bioinformatics/variantsCalled_samtoolsBcftools/geneModels/",geneName,".gff", sep = "")
#Store gff3 file as TxDb object
geneDB <- makeTxDbFromGFF(file=inputGFF, format="gff3", organism="Daphnia pulex")

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

#Generate track overlay
#ot <- OverlayTrack(trackList = list(alTrE05, alTrR2, alTrY05, alTrY023))

#Generate data track of coverages
#dt <- DataTrack(alTrUV1, alTrUV2, alTrUV3, alTrVIS1, alTrVIS2, alTrVIS3)

# plot all data tracks
jpeg(file=paste(outPath,"/filteredCoverage_geneModel_",geneName,".jpeg", sep = ""), width=900, height=400)
plotTracks(list(sTrack, geneTr, alTrE05, alTrR2, alTrY05, alTrY023), 
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
#cdsStart = 364812
#cdsEnd = 365085

# plot cds 1
#plotTracks(c(geneTr, sTrack, alTrE05, alTrR2, alTrY05, alTrY023), 
           #from = cdsStart,
           #to = cdsEnd, 
           #chromosome = geneLoc,
           #exonAnnotation = "feature",
           #fontcolor.exon = 1,
           #fontsize.exon = 8,
           #thinBoxFeature = "UTR", 
           #collapse = FALSE,
           #alpha.title = 1,
           #showIndels=TRUE)
