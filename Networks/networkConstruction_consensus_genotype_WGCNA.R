#!/usr/bin/env Rscript

# script to construct a network for a set of samples using WGNCA
# usage: Rscript networkConstruction_consensus_WGCNA.R workingDir softPower minModuleSize
# usage ex: Rscript networkConstruction_consensus_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_genotype_WGCNA 30 100
# usage ex: Rscript networkConstruction_consensus_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_genotype_WGCNA 14 30
# usage ex: Rscript networkConstruction_consensus_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_genotype_WGCNA 14 60

#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)

#Set working directory
workingDir = args[1];
#workingDir="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_genotype_WGCNA"
setwd(workingDir)

# load the WGCNA package
library(WGCNA)

# the following setting is important, do not omit.
options(stringsAsFactors = FALSE)

# Allow multi-threading within WGCNA. At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads()

# Load the data saved in the first part
setTag <- "genotype"
importFile <- paste("Consensus-dataInput", setTag, sep="-")
importFile <- paste(importFile, "RData", sep=".")
lnames = load(file = importFile)
#The variable lnames contains the names of loaded variables.
lnames
# Get the number of sets in the multiExpr structure.
nSets = checkSets(multiExpr)$nSets

# set soft thresholding power
softPower = as.numeric(args[2])
#softPower = 30
#softPower = 14

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = as.numeric(args[3])
#minModuleSize = 60

# Initialize an appropriate array to hold the adjacencies
adjacencies = array(0, dim = c(nSets, nGenes, nGenes))
# Calculate adjacencies in each individual data set
for (set in 1:nSets){
  adjacencies[set, , ] = abs(cor(multiExpr[[set]]$data, use = "p"))^softPower
}

# Initialize an appropriate array to hold the TOMs
TOM = array(0, dim = c(nSets, nGenes, nGenes))
# Calculate TOMs in each individual data set
for (set in 1:nSets){
  TOM[set, , ] = TOMsimilarity(adjacencies[set, , ])
}

# Define the reference percentile
scaleP = 0.95
# Set RNG seed for reproducibility of sampling
set.seed(12345)
# Sample sufficiently large number of TOM entries
nSamples = as.integer(1/(1-scaleP) * 1000)
# Choose the sampled TOM entries
scaleSample = sample(nGenes*(nGenes-1)/2, size = nSamples)
TOMScalingSamples = list()
# These are TOM values at reference percentile
scaleQuant = rep(1, nSets)
# Scaling powers to equalize reference TOM values
scalePowers = rep(1, nSets)
# Loop over sets
for (set in 1:nSets){
  # Select the sampled TOM entries
  TOMScalingSamples[[set]] = as.dist(TOM[set, , ])[scaleSample]
  # Calculate the 95th percentile
  scaleQuant[set] = quantile(TOMScalingSamples[[set]],
                             probs = scaleP, type = 8)
  # Scale the male TOM
  if (set>1){
    scalePowers[set] = log(scaleQuant[1])/log(scaleQuant[set])
    TOM[set, ,] = TOM[set, ,]^scalePowers[set]
  }
}

# For plotting, also scale the sampled TOM entries
scaledTOMSamples = list()
for (set in 1:nSets){
  scaledTOMSamples[[set]] = TOMScalingSamples[[set]]^scalePowers[set]
}

# genotype set
# Y05 vs Y023
# Open a suitably sized graphics window
png(file = "ConsensusTOMScaling-QQPlot_Y05_Y023.png", wi = 6, he = 6, units="in", res=150)
sizeGrWindow(6,6)
# qq plot of the unscaled samples
qqUnscaled = qqplot(TOMScalingSamples[[1]], TOMScalingSamples[[2]], plot.it = TRUE, cex = 0.6,
                    xlab = paste("TOM in", setLabels[1]), ylab = paste("TOM in", setLabels[2]),                    
                    main = "Q-Q plot of TOM", pch = 20)
# qq plot of the scaled samples
qqScaled = qqplot(scaledTOMSamples[[1]], scaledTOMSamples[[2]], plot.it = FALSE)
points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20)
abline(a=0, b=1, col = "blue")
legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch = 20, col = c("black", "red"))
dev.off()
# E05 vs R2
# Open a suitably sized graphics window
png(file = "ConsensusTOMScaling-QQPlot_E05_R2.png", wi = 6, he = 6, units="in", res=150)
sizeGrWindow(6,6)
# qq plot of the unscaled samples
qqUnscaled = qqplot(TOMScalingSamples[[3]], TOMScalingSamples[[4]], plot.it = TRUE, cex = 0.6,
                    xlab = paste("TOM in", setLabels[3]), ylab = paste("TOM in", setLabels[4]),                    
                    main = "Q-Q plot of TOM", pch = 20)
# qq plot of the scaled samples
qqScaled = qqplot(scaledTOMSamples[[3]], scaledTOMSamples[[4]], plot.it = FALSE)
points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20)
abline(a=0, b=1, col = "blue")
legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch = 20, col = c("black", "red"))
dev.off()
# Y05 vs E05
# Open a suitably sized graphics window
png(file = "ConsensusTOMScaling-QQPlot_Y05_E05.png", wi = 6, he = 6, units="in", res=150)
sizeGrWindow(6,6)
# qq plot of the unscaled samples
qqUnscaled = qqplot(TOMScalingSamples[[1]], TOMScalingSamples[[3]], plot.it = TRUE, cex = 0.6,
                    xlab = paste("TOM in", setLabels[1]), ylab = paste("TOM in", setLabels[3]),                   
                    main = "Q-Q plot of TOM", pch = 20)
# qq plot of the scaled samples
qqScaled = qqplot(scaledTOMSamples[[1]], scaledTOMSamples[[3]], plot.it = FALSE)
points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20)
abline(a=0, b=1, col = "blue")
legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch = 20, col = c("black", "red"))
dev.off()
# Y023 vs R2
# Open a suitably sized graphics window
png(file = "ConsensusTOMScaling-QQPlot_Y023_R2.png", wi = 6, he = 6, units="in", res=150)
sizeGrWindow(6,6)
# qq plot of the unscaled samples
qqUnscaled = qqplot(TOMScalingSamples[[2]], TOMScalingSamples[[4]], plot.it = TRUE, cex = 0.6,
                    xlab = paste("TOM in", setLabels[2]), ylab = paste("TOM in", setLabels[4]),                    
                    main = "Q-Q plot of TOM", pch = 20)
# qq plot of the scaled samples
qqScaled = qqplot(scaledTOMSamples[[2]], scaledTOMSamples[[4]], plot.it = FALSE)
points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20)
abline(a=0, b=1, col = "blue")
legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch = 20, col = c("black", "red"))
dev.off()

# calculate the consensus Topological Overlap by taking the component-wise (“parallel”) 
# minimum of the TOMs in individual sets
consensusTOM = pmin(TOM[1, , ], TOM[2, , ])

# Clustering
consTree = hclust(as.dist(1-consensusTOM), method = "average")

# Module identification using dynamic tree cut:
unmergedLabels = cutreeDynamic(dendro = consTree, distM = 1-consensusTOM,
                               deepSplit = 2, cutHeight = 0.995,
                               minClusterSize = minModuleSize,
                               pamRespectsDendro = FALSE )
unmergedColors = labels2colors(unmergedLabels)

# quick summary of the module detection
table(unmergedLabels)

# plot the consensus gene dendrogram together with the preliminary module colors
exportFile <- paste("ConsensusDendrogram", minModuleSize, sep="_")
exportFile <- paste(exportFile, "png", sep=".")
png(file = exportFile, wi = 8, he = 6, units="in", res=150)
sizeGrWindow(8,6)
plotDendroAndColors(consTree, unmergedColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Calculate module eigengenes
unmergedMEs = multiSetMEs(multiExpr, colors = NULL, universalColors = unmergedColors)
# Calculate consensus dissimilarity of consensus module eigengenes
consMEDiss = consensusMEDissimilarity(unmergedMEs);
# Cluster consensus modules
consMETree = hclust(as.dist(consMEDiss), method = "average");
# Plot the result
exportFile <- paste("ConsensusModuleClustering", minModuleSize, sep="_")
exportFile <- paste(exportFile, "png", sep=".")
png(file = exportFile, wi = 7, he = 6, units="in", res=150)
sizeGrWindow(7,6)
par(mfrow = c(1,1))
plot(consMETree, main = "Consensus clustering of consensus module eigengenes",
     xlab = "", sub = "")
abline(h=0.25, col = "red")
dev.off()

# automatically merge module bellow the merging threshold
merge = mergeCloseModules(multiExpr, unmergedLabels, cutHeight = 0.25, verbose = 3)

# Numeric module labels
moduleLabels = merge$colors
# Convert labels to colors
moduleColors = labels2colors(moduleLabels)
# Eigengenes of the new merged modules:
consMEs = merge$newMEs

#  plot the gene dendrogram again, this time with both the 
# unmerged and the merged module colors
exportFile <- paste("ConsensusDendrogram_mergedModules", minModuleSize, sep="_")
exportFile <- paste(exportFile, "png", sep=".")
png(file = exportFile, wi = 9, he = 6, units="in", res=150)
sizeGrWindow(9,6)
plotDendroAndColors(consTree, cbind(unmergedColors, moduleColors),
                    c("Unmerged", "Merged"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# save the information necessary in the subsequent steps
exportFile <- paste("Consensus-NetworkConstruction-man", minModuleSize, sep="-")
exportFile <- paste(exportFile, "RData", sep=".")
save(consMEs, moduleColors, moduleLabels, consTree, file=exportFile) 
