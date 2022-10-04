#Retrieve input file name of gene counts
#args = commandArgs(trailingOnly=TRUE)

#Set working directory
#workingDir = args[1];
workingDir="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_WGCNA"
setwd(workingDir)

# load the WGCNA package
library(WGCNA)

# the following setting is important, do not omit.
options(stringsAsFactors = FALSE)

# Allow multi-threading within WGCNA. At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
#enableWGCNAThreads()

# retrieve genotype tag
#genotype <- args[2]
genotype <- "Y05"

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
subsetModuleLabels = substring(names(subsetMEs), 3)
consModuleLabels = substring(names(consMEs[[1]]$data), 3)
# Convert the numeric module labels to color labels
subsetModules = labels2colors(as.numeric(subsetModuleLabels))
consModules = labels2colors(as.numeric(subsetModuleLabels))
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

