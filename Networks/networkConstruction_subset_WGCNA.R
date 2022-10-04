#!/usr/bin/env Rscript

# script to create a network for a subset of samples using WGNCA
# usage: Rscript networkConstruction_subset_WGCNA.R workingDir countsFile startCounts endCounts startSubset endSubset
# usage ex: Rscript networkConstruction_subset_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_WGCNA Y05 20
# usage ex: Rscript networkConstruction_subset_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_WGCNA Y023 12
# usage ex: Rscript networkConstruction_subset_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_WGCNA E05 14
# usage ex: Rscript networkConstruction_subset_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_WGCNA R2 9

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
workingDir = args[1];
#workingDir="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_WGCNA"
setwd(workingDir)

# retrieve genotype tag
genotype <- args[2]
#genotype <- "Y05"

#Import normalized gene count data
# Load the data saved in the first part
importFile <- paste(genotype, "dataInput.RData", sep="-")
lnames = load(file = importFile)

# set the soft thresholding power
softPower = as.numeric(args[3])
# Y05
#softPower = 20
# Y023
#softPower = 12
# E05
#softPower = 14
# R2
#softPower = 9

# determine adjacency
adjacency = adjacency(datExpr, power = softPower)

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")
# Plot the resulting clustering tree (dendrogram)
exportFile <- paste(genotype, "geneClustering.pdf", sep="_")
pdf(file = exportFile)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 100
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
exportFile <- paste(genotype, "dynamicTreeCut.pdf", sep="_")
pdf(file = exportFile)
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
exportFile <- paste(genotype, "dynamicTreeCut.pdf", sep="_")
pdf(file = exportFile)
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
# choose a height cut of 0.25, corresponding to correlation of 0.75, to merge
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()

# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

# plot the gene dendrogram again, with the 
# original and merged module colors underneath
exportFile <- paste(genotype, "geneDendro-3.pdf", sep="_")
pdf(file = exportFile, wi = 9, he = 6)
sizeGrWindow(12, 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
exportFile <- paste(genotype, "networkConstruction-stepByStep.RData", sep="-")
save(MEs, moduleLabels, moduleColors, geneTree, file = exportFile)
