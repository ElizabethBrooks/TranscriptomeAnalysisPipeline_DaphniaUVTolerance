# script to identify set-specific and consensos modules
# usage: Rscript relateModsToTraits_consensus_WGCNA.R workingDir minModuleSize setTag
# usage ex: Rscript relateModsToTraits_consensus_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/ensembl/GCA_021134715.1/biostatistics/NetworkAnalysis/WGCN_tolerance_WGCNA 60 tolerance

#Retrieve input file name of gene counts
#args = commandArgs(trailingOnly=TRUE)

#Set working directory
#workingDir = args[1];
workingDir="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/ensembl/GCA_021134715.1/biostatistics/NetworkAnalysis/WGCN_tolerance_WGCNA"
setwd(workingDir)

# load the WGCNA package
library(WGCNA)

# the following setting is important, do not omit.
options(stringsAsFactors = FALSE)

# We like large modules, so we set the minimum module size relatively high:
#minModuleSize = as.numeric(args[2])
minModuleSize = 60

# retrieve set tag
#setTag <- args[3]
setTag <- "tolerance"

# load the input consensus data
importFile <- paste("Consensus-dataInput", setTag, sep="-")
importFile <- paste(importFile, "RData", sep=".")
lnames = load(file = importFile)
lnames

# Also load results of network analysis
importFile <- paste("Consensus-NetworkConstruction-man", minModuleSize, sep="-")
importFile <- paste(importFile, "RData", sep=".")
lnames = load(file = importFile);
lnames
exprSize = checkSets(multiExpr)
nSets = exprSize$nSets

# Set up variables to contain the module-trait correlations
moduleTraitCor = list()
moduleTraitPvalue = list()
# Calculate the correlations
for (set in 1:nSets){
  moduleTraitCor[[set]] = cor(consMEs[[set]]$data, Traits[[set]]$data, use = "p");
  moduleTraitPvalue[[set]] = corPvalueFisher(moduleTraitCor[[set]], exprSize$nSamples[set]);
}

# Convert numerical lables to colors for labeling of modules in the plot
MEColors = labels2colors(as.numeric(substring(names(consMEs[[1]]$data), 3)));
MEColorNames = paste("ME", MEColors, sep="");
# Open a suitably sized window (the user should change the window size if necessary)
sizeGrWindow(10,7)
#png(file = "Plots/ModuleTraitRelationships-female.png", wi = 10, he = 7, units="in", res=150);
# Plot the module-trait relationship table for set number 1
set = 1
textMatrix = paste(signif(moduleTraitCor[[set]], 2), "\n(",
                   signif(moduleTraitPvalue[[set]], 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[[set]])
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = moduleTraitCor[[set]],
               xLabels = names(Traits[[set]]$data),
               yLabels = MEColorNames,
               ySymbols = MEColorNames,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module--trait relationships in", setLabels[set]))
dev.off();
# Plot the module-trait relationship table for set number 2
set = 2
textMatrix = paste(signif(moduleTraitCor[[set]], 2), "\n(",
                   signif(moduleTraitPvalue[[set]], 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[[set]])
sizeGrWindow(10,7)
#png(file = "Plots/ModuleTraitRelationships-male.png", wi = 10, he = 7, units="in", res=150);
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = moduleTraitCor[[set]],
               xLabels = names(Traits[[set]]$data),
               yLabels = MEColorNames,
               ySymbols = MEColorNames,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module--trait relationships in", setLabels[set]))
dev.off();

# Initialize matrices to hold the consensus correlation and p-value
consensusCor = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]));
consensusPvalue = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]));
# Find consensus negative correlations
negative = moduleTraitCor[[1]] < 0 & moduleTraitCor[[2]] < 0;
consensusCor[negative] = pmax(moduleTraitCor[[1]][negative], moduleTraitCor[[2]][negative]);
consensusPvalue[negative] = pmax(moduleTraitPvalue[[1]][negative], moduleTraitPvalue[[2]][negative]);
# Find consensus positive correlations
positive = moduleTraitCor[[1]] > 0 & moduleTraitCor[[2]] > 0;
consensusCor[positive] = pmin(moduleTraitCor[[1]][positive], moduleTraitCor[[2]][positive]);
consensusPvalue[positive] = pmax(moduleTraitPvalue[[1]][positive], moduleTraitPvalue[[2]][positive]);

# display the consensus module–trait relationships again using a color-coded table
textMatrix = paste(signif(consensusCor, 2), "\n(",
                   signif(consensusPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[[set]])
sizeGrWindow(10,7)
#png(file = "Plots/ModuleTraitRelationships-consensus.png", wi = 10, he = 7, units="in", res=150);
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = consensusCor,
               xLabels = names(Traits[[set]]$data),
               yLabels = MEColorNames,
               ySymbols = MEColorNames,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Consensus module--trait relationships across\n",
                            paste(setLabels, collapse = " and ")))
dev.off();
