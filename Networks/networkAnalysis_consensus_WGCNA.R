#!/usr/bin/env Rscript

# script to identify set-specific and consensos modules
# usage: Rscript networkAnalysis_consensus_WGCNA.R workingDir countsFile subsetTag minModuleSize
# usage ex: Rscript networkAnalysis_consensus_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA tol 60 tolerance
# usage ex: Rscript networkAnalysis_consensus_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA nTol 60 tolerance
# usage ex: Rscript networkAnalysis_consensus_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA Y05 100 tolerance
# usage ex: Rscript networkAnalysis_consensus_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA Y023 100 tolerance
# usage ex: Rscript networkAnalysis_consensus_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA E05 100 tolerance
# usage ex: Rscript networkAnalysis_consensus_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA R2 100 tolerance

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

# retrieve subset tag
subsetTag <- args[2]
#subsetTag <- "Y05"

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = as.numeric(args[3])
#minModuleSize = 60

# retrieve the set tag
setTag <- args[4]

# load the subset data and rename them so that 
# they do not conflict with the consensus data
importFile <- paste(subsetTag, minModuleSize, sep="_")
importFile <- paste(importFile, "networkConstruction-stepByStep.RData", sep="-")
lnames = load(file = importFile)
lnames

# Rename variables to avoid conflicts
subsetLabels = moduleLabels
subsetColors = moduleColors
subsetTree = geneTree
subsetMEs = orderMEs(MEs, greyName = "ME0")

# load the results of the consensus module identification:
importFile <- paste("Consensus-NetworkConstruction-man", minModuleSize, sep="-")
importFile <- paste(importFile, "RData", sep=".")
lnames = load(importFile)
lnames

# restrict modules to genes that occur in both sets
importFile <- paste("Consensus-dataInput", setTag, sep="-")
importFile <- paste(importFile, "RData", sep=".")
load(file = importFile)
importFile <- paste(subsetTag, "dataInput.RData", sep="-")
load(file = importFile)
subsetGenes = colnames(datExpr)
consGenes = mtd.colnames(multiExpr)
common = intersect(subsetGenes, consGenes)
subsetColors.common = subsetColors[match(common, subsetGenes)];
moduleColors.common = moduleColors[match(common, consGenes)];

# Isolate the module labels in the order they appear in ordered module eigengenes
#subsetModuleLabels = substring(names(subsetMEs), 3)
consModuleLabels = substring(names(consMEs[[1]]$data), 3)
subsetModules = substring(names(subsetMEs), 3)
#consModules = substring(names(consMEs[[1]]$data), 3)
# Convert the numeric module labels to color labels
#subsetModules = labels2colors(as.numeric(subsetModuleLabels))
consModules = labels2colors(as.numeric(consModuleLabels))
# Numbers of subset and consensus modules
nSubsetMods = length(subsetModules)
nConsMods = length(consModules)
# Initialize tables of p-values and of the corresponding counts
pTable = matrix(0, nrow = nSubsetMods, ncol = nConsMods)
CountTbl = matrix(0, nrow = nSubsetMods, ncol = nConsMods)
# Execute all pairwaise comparisons
for (smod in 1:nSubsetMods){
  for (cmod in 1:nConsMods){
    subsetMembers = (subsetColors.common == subsetModules[smod])
    consMembers = (moduleColors.common == consModules[cmod])
    pTable[smod, cmod] = -log10(fisher.test(subsetMembers, consMembers, alternative = "greater")$p.value)
    CountTbl[smod, cmod] = sum(subsetColors.common == subsetModules[smod] & moduleColors.common ==
                                 consModules[cmod])
  }
}

# Truncate p values smaller than 10^{-50} to 10^{-50}
pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)])
pTable[pTable>50 ] = 50
# Marginal counts (really module sizes)
subsetModTotals = apply(CountTbl, 1, sum)
consModTotals = apply(CountTbl, 2, sum)
# setup plot title
plotTitle <- paste("Correspondence of", subsetTag, "Set-Specific and Consensus Modules", sep=" ")
# Actual plotting
exportFile <- paste(subsetTag, "ConsensusVsSubsetModules", sep="_")
exportFile <- paste(exportFile, minModuleSize, sep="_")
exportFile <- paste(exportFile, "png", sep=".")
png(file = exportFile, wi = 10, he = 7, units="in", res=150)
sizeGrWindow(10,7 )
par(mfrow=c(1,1))
par(cex = 1.0)
par(mar=c(8, 10.4, 2.7, 1)+0.3)
# Use function labeledHeatmap to produce the color-coded table with all the trimmings
labeledHeatmap(Matrix = pTable,
               xLabels = paste(" ", consModules),
               yLabels = paste(" ", subsetModules),
               colorLabels = TRUE,
               xSymbols = paste("Cons ", consModules, ": ", consModTotals, sep=""),
               ySymbols = paste("Subset ", subsetModules, ": ", subsetModTotals, sep=""),
               textMatrix = CountTbl,
               colors = blueWhiteRed(100)[50:100],
               main = plotTitle,
               cex.text = 1.0, cex.lab = 1.0, setStdMargins = FALSE)
dev.off()
