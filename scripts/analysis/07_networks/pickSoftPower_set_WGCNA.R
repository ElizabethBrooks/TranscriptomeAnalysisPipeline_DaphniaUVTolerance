#!/usr/bin/env Rscript

# script to create a network for a subset of samples using WGNCA
# usage: Rscript pickSoftPower_set_WGCNA.R workingDir subsetTag
# usage ex: Rscript pickSoftPower_set_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/ensembl/GCA_021134715.1/biostatistics/NetworkAnalysis/WGCN_tolerance_WGCNA OLYM
# usage ex: Rscript pickSoftPower_set_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/ensembl/GCA_021134715.1/biostatistics/NetworkAnalysis/WGCN_tolerance_WGCNA tol
# usage ex: Rscript pickSoftPower_set_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/ensembl/GCA_021134715.1/biostatistics/NetworkAnalysis/WGCN_tolerance_WGCNA nTol
# usage ex: Rscript pickSoftPower_set_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/ensembl/GCA_021134715.1/biostatistics/NetworkAnalysis/WGCN_tolerance_WGCNA Y05
# usage ex: Rscript pickSoftPower_set_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/ensembl/GCA_021134715.1/biostatistics/NetworkAnalysis/WGCN_tolerance_WGCNA Y023
# usage ex: Rscript pickSoftPower_set_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/ensembl/GCA_021134715.1/biostatistics/NetworkAnalysis/WGCN_tolerance_WGCNA E05
# usage ex: Rscript pickSoftPower_set_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/ensembl/GCA_021134715.1/biostatistics/NetworkAnalysis/WGCN_tolerance_WGCNA R2

#Load the WGCNA and edgeR packages
library(WGCNA)

#The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

# Allow multi-threading within WGCNA. At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads()

#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)

#Set working directory
workingDir = args[1]
#workingDir="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_WGCNA"
setwd(workingDir)

# retrieve genotype tag
genotype <- args[2]
#genotype <- "Y05"

#Import normalized gene count data
# Load the data saved in the first part
importFile <- paste(genotype, "dataInput.RData", sep="-")
lnames = load(file = importFile)

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=36, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results
cex1 = 0.9
exportFile <- paste(genotype, "SoftPowers.png", sep="_")
png(file = exportFile, wi = 9, he = 5, units="in", res=150)
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
abline(h=0.90,col="blue")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
