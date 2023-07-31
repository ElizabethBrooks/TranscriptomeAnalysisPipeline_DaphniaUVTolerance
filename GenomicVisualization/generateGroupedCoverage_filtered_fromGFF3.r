#!/usr/bin/env Rscript
#R script to generate a gene model data track from a gff3 file
#Usage: Rscript generateGeneModel_fromGFF3.r gffPath scaffoldName geneName geneStart geneEnd bamFile

#Install Rsubread using BiocManager, this should only need to be done once
#BiocManager::install("GenomicRanges")
#Load the Rsubread library for featureCounts
library(GenomicFeatures)
library(Gviz)
library(Rsamtools)

#Retrieve input arguments
#args = commandArgs(trailingOnly=TRUE)

#Turn off UCSC chromosome names to switch to custom naming
options(ucscChromosomeNames=FALSE)

#Retrieve input gff file path
inputGFF = "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Bioinformatics/variantsCalled_samtoolsBcftools/geneModels/LOC124188748.gff"
#Store gff3 file as TxDb object
geneDB <- makeTxDbFromGFF(file=inputGFF, format="gff3", organism="Daphnia pulex")

#Retrieve output path
outPath = "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Bioinformatics/variantsCalled_samtoolsBcftools/geneModels"

#Set gene region name and coordinates
#Random ex: scaffold_1 gene1 2873-4893, scaffold_1 gene10 189307-194728, scaffold_1 gene100 1175334-1176811, scaffold_2 gene1319 95-3193
#High DE ex: scaffold_9 gene6568 2352829-2354488, scaffold_97 gene18147 72973-74469, scaffold_48 gene15097 51894-56155, scaffold_2 gene1861 3631041-3632552, 
# scaffold_6 gene4926 777761-779275, scaffold_8 gene5759 504618-506995
geneLoc = "NC_060018.1"
geneName = "LOC124188748"
geneStart = 363588
geneEnd = 367123

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
alTrE05 <- AlignmentsTrack(inputBamE05, isPaired = TRUE, name = "E05", type = "coverage", col = "blue", fill = "blue")
alTrR2 <- AlignmentsTrack(inputBamR2, isPaired = TRUE, name = "R2", type = "coverage", col = "purple", fill = "purple")
alTrY05 <- AlignmentsTrack(inputBamY05, isPaired = TRUE, name = "Y05", type = "coverage", col = "orange", fill = "orange")
alTrY023 <- AlignmentsTrack(inputBamY023, isPaired = TRUE, name = "Y023", type = "coverage", col = "red", fill = "red")

#Generate track overlay
#ot <- OverlayTrack(trackList = list(alTrE05, alTrR2, alTrY05, alTrY023))

#Generate data track of coverages
#dt <- DataTrack(alTrUV1, alTrUV2, alTrUV3, alTrVIS1, alTrVIS2, alTrVIS3)

#Plot all data tracks
#jpeg(file=paste(outPath,"/filteredCoverage_",geneName,".jpeg", sep = ""), width=900, height=400)
plotTracks(c(geneTr, alTrE05, alTrR2, alTrY05, alTrY023), 
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
           alpha.title = 1)
           #stackedBars = TRUE,
           #col.coverage = "green",
           #fill.coverage = "blue",
           #col.mates = "purple", 
           #col.gap = "orange"
#dev.off()
