#!/usr/bin/env Rscript

# script to identify set-specific and consensos modules
# usage: Rscript eigengeneNetworks_consensus_tolerance_WGCNA.R workingDir minModuleSize
# usage ex: Rscript eigengeneNetworks_consensus_tolerance_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/ensembl/GCA_021134715.1/biostatistics/NetworkAnalysis/WGCN_tolerance_WGCNA 60

#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)

#Set working directory
workingDir = args[1];
#workingDir="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/ensembl/GCA_021134715.1/biostatistics/NetworkAnalysis/WGCN_tolerance_WGCNA"
setwd(workingDir)

# load the WGCNA package
library(WGCNA)

# the following setting is important, do not omit.
options(stringsAsFactors = FALSE)

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = as.numeric(args[2])
#minModuleSize = 60

# retrieve the set tag
setTag <- "tolerance"

# number of subsets
nSets <- 2

# For easier labeling of plots, create a vector holding descriptive names of the sets
setLabels = c("Tolerant", "Not Tolerant")
shortLabels = c("tol", "nTol")

# Load the data saved in the first part
importFile <- paste("Consensus-dataInput", setTag, sep="-")
importFile <- paste(importFile, "RData", sep=".")
lnames = load(file = importFile);
#The variable lnames contains the names of loaded variables.
lnames

# Load the results of network analysis
importFile <- paste("Consensus-NetworkConstruction-man", minModuleSize, sep="-")
importFile <- paste(importFile, "RData", sep=".")
lnames = load(file = importFile)
lnames

# Create a variable treatment that will hold just the treatment trait in both sets
treatment = vector(mode = "list", length = nSets)
for (set in 1:nSets){
treatment[[set]] = list(data = as.data.frame(Traits[[set]]$data$treatment))
names(treatment[[set]]$data) = "treatment"
}
# Recalculate consMEs to give them color names
consMEsC = multiSetMEs(multiExpr, universalColors = moduleColors)
# We add the treatment trait to the eigengenes and order them by consesus hierarchical clustering:
MET = consensusOrderMEs(addTraitToMEs(consMEsC, treatment))

# now call the function plotEigengeneNetworks that performs the differential analysis
sizeGrWindow(8,10)
exportFile <- paste("EigengeneNetworks", minModuleSize, sep="_")
exportFile <- paste(exportFile, "png", sep=".")
png(file = exportFile, width= 8, height = 10, units="in", res=150)
par(cex = 0.9)
plotEigengeneNetworks(MET, setLabels, marDendro = c(0,2,2,1), marHeatmap = c(3,3,2,1),
zlimPreservation = c(0.5, 1), xLabelsAngle = 90)
dev.off()
