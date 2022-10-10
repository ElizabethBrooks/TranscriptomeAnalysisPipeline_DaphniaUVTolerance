#!/usr/bin/env Rscript

# script to create a network for a subset of samples using WGNCA
# usage: Rscript dataInput_subset_WGCNA.R workingDir countsFile startCounts endCounts startSubset endSubset subsetTag traitsFile
# usage ex: Rscript dataInput_subset_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/DEGenotypes/glmQLF_normalizedCounts.csv 1 24 1 12 tol /Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_treatment_WGCNA_Olympics.csv
# usage ex: Rscript dataInput_subset_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/DEGenotypes/glmQLF_normalizedCounts.csv 1 24 13 24 nTol /Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_treatment_WGCNA_Olympics.csv
# usage ex: Rscript dataInput_subset_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/DEGenotypes/glmQLF_normalizedCounts.csv 1 24 1 6 Y05 /Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_treatment_WGCNA_Olympics.csv
# usage ex: Rscript dataInput_subset_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/DEGenotypes/glmQLF_normalizedCounts.csv 1 24 7 12 Y023 /Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_treatment_WGCNA_Olympics.csv
# usage ex: Rscript dataInput_subset_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/DEGenotypes/glmQLF_normalizedCounts.csv 1 24 13 18 E05 /Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_treatment_WGCNA_Olympics.csv
# usage ex: Rscript dataInput_subset_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/DEGenotypes/glmQLF_normalizedCounts.csv 1 24 19 24 R2 /Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_treatment_WGCNA_Olympics.csv
# usage ex: Rscript dataInput_subset_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_genotype_WGCNA /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/DEGenotypes/glmQLF_normalizedCounts.csv 1 24 1 12 tol /Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_treatment_WGCNA_Olympics.csv
# usage ex: Rscript dataInput_subset_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_genotype_WGCNA /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/DEGenotypes/glmQLF_normalizedCounts.csv 1 24 13 24 nTol /Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_treatment_WGCNA_Olympics.csv
# usage ex: Rscript dataInput_subset_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_genotype_WGCNA /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/DEGenotypes/glmQLF_normalizedCounts.csv 1 24 1 6 Y05 /Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_treatment_WGCNA_Olympics.csv
# usage ex: Rscript dataInput_subset_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_genotype_WGCNA /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/DEGenotypes/glmQLF_normalizedCounts.csv 1 24 7 12 Y023 /Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_treatment_WGCNA_Olympics.csv
# usage ex: Rscript dataInput_subset_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_genotype_WGCNA /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/DEGenotypes/glmQLF_normalizedCounts.csv 1 24 13 18 E05 /Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_treatment_WGCNA_Olympics.csv
# usage ex: Rscript dataInput_subset_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_genotype_WGCNA /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/DEGenotypes/glmQLF_normalizedCounts.csv 1 24 19 24 R2 /Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_treatment_WGCNA_Olympics.csv

#Load the WGCNA and edgeR packages
library(WGCNA)

#The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)

#Set working directory
workingDir = args[1];
#workingDir="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_genotype_WGCNA"
setwd(workingDir)

#Import normalized gene count data
inputTable <- read.csv(file=args[2], row.names="gene")[ ,args[3]:args[4]]
#inputTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/DEGenotypes/glmQLF_normalizedCounts.csv", row.names="gene", header=TRUE)[ ,1:24]

#Subset input counts by genotype or factor
inputTable_subset <- inputTable[,args[5]:args[6]]
#inputTable_subset <- inputTable[,1:6]

# retrieve subsetTag tag
tag <- args[7]
#tag <- "Y05"

# load in the trait data
allTraits = read.csv(args[8])
#allTraits = read.csv("/Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_consensus_WGCNA_Olympics.csv")
dim(allTraits)
names(allTraits)

# transpose each subset
datExpr0 = data = as.data.frame(t(inputTable_subset))
names(datExpr0) = rownames(inputTable_subset)
rownames(datExpr0) = names(inputTable_subset)

#Check the genes across all samples
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

# remove the offending genes and samples from the data
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# cluster the samples to see if there are any obvious outliers
sampleTree = hclust(dist(datExpr0), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
exportFile <- paste(tag, "sampleClustering.png", sep="_")
png(file = exportFile, width = 12, height = 9, units="in", res=150)
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
#abline(h = 15, col = "red")
dev.off()

# Determine cluster under the line
#clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
#table(clust)
# clust 1 contains the samples we want to keep.
#keepSamples = (clust==1)

# filter out samples
datExpr <- datExpr0
#datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# Form a data frame analogous to expression data that will hold the clinical traits
samples = rownames(datExpr)
traitRows = match(samples, allTraits$sample)
datTraits = allTraits[traitRows, -1]
rownames(datTraits) = allTraits[traitRows, 1]
collectGarbage()

# Re-cluster samples
exportFile <- paste(tag, "sampleDendrogram_traitHeatmap.png", sep="_")
png(file = exportFile, units="in", res=150)
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

# save the relevant expression and trait data for use in the next steps
exportFile <- paste(tag, "dataInput.RData", sep="-")
save(datExpr, datTraits, file = exportFile)
