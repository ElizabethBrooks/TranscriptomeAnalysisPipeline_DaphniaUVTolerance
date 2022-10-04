# script to identify set-specific and consensos modules
# usage: Rscript networkAnalysis_consensus_WGCNA.R workingDir countsFile startCounts endCounts startSubset endSubset
# usage ex: Rscript networkAnalysis_consensus_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_WGCNA Y05
# usage ex: Rscript networkAnalysis_consensus_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_WGCNA Y023
# usage ex: Rscript networkAnalysis_consensus_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_WGCNA E05
# usage ex: Rscript networkAnalysis_consensus_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_WGCNA R2

#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)

#Set working directory
workingDir = args[1];
#workingDir="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_WGCNA"
setwd(workingDir)

# load the WGCNA package
library(WGCNA)

# the following setting is important, do not omit.
options(stringsAsFactors = FALSE)

# retrieve genotype tag
genotype <- args[2]
#genotype <- "Y05"

# load the subset data and rename them so that 
# they do not conflict with the consensus data
importFile <- paste(genotype, "networkConstruction-stepByStep.RData", sep="-")
lnames = load(file = importFile)
lnames

# Rename variables to avoid conflicts
subsetLabels = moduleLabels
subsetColors = moduleColors
subsetTree = geneTree
subsetMEs = orderMEs(MEs, greyName = "ME0")

# load the results of the consensus module identification:
lnames = load("Consensus-NetworkConstruction-man.RData")
lnames

# restrict modules to genes that occur in both sets
load(file = "Consensus-dataInput.RData")
importFile <- paste(genotype, "dataInput.RData", sep="-")
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
pTable = matrix(0, nrow = nSubsetMods, ncol = nConsMods);
CountTbl = matrix(0, nrow = nSubsetMods, ncol = nConsMods);
# Execute all pairwaise comparisons
for (smod in 1:nSubsetMods){
  for (cmod in 1:nConsMods){
    subsetMembers = (subsetColors.common == subsetModules[smod]);
    consMembers = (moduleColors.common == consModules[cmod]);
    pTable[smod, cmod] = -log10(fisher.test(subsetMembers, consMembers, alternative = "greater")$p.value);
    CountTbl[smod, cmod] = sum(subsetColors.common == subsetModules[smod] & moduleColors.common ==
                                 consModules[cmod])
  }
}

# Truncate p values smaller than 10^{-50} to 10^{-50}
pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
pTable[pTable>50 ] = 50 ;
# Marginal counts (really module sizes)
subsetModTotals = apply(CountTbl, 1, sum)
consModTotals = apply(CountTbl, 2, sum)
# setup plot title
plotTitle <- paste("Correspondence of", genotype, "Set-Specific and Consensus Modules", sep=" ")
# Actual plotting
exportFile <- paste(genotype, "ConsensusVsSubsetModules.pdf", sep="_")
pdf(file = exportFile, wi = 10, he = 7);
sizeGrWindow(10,7 );
par(mfrow=c(1,1));
par(cex = 1.0);
par(mar=c(8, 10.4, 2.7, 1)+0.3);
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
               cex.text = 1.0, cex.lab = 1.0, setStdMargins = FALSE);
dev.off()
