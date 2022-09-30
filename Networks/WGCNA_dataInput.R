#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("edgeR")

#Load the WGCNA and edgeR packages
library(WGCNA)
library(edgeR)

#The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)

#Set working directory
#workingDir = args[1];
workingDir="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4"
setwd(workingDir)

#Import normalized gene count data
#inputTable <- read.csv(file=args[2], row.names="gene")[ ,args[3]:args[4]]
inputTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/DEGenotypes/glmQLF_normalizedCounts.csv", row.names="gene", header=TRUE)[ ,1:24]

#Subset input counts by genotype
inputTable_Y05 <- inputTable[,1:6]
inputTable_Y023 <- inputTable[,7:12]
inputTable_E05 <- inputTable[,13:18]
inputTable_R2 <- inputTable[,19:24]

# We work with four sets
nSets <- 4
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("Y05 Tolerant", "Y023 Tolerant", "E05 Not Tolerant", "R2 Not Tolerant")
shortLabels = c("Y05", "Y023", "E05", "R2")


# form multi-set expression data
multiExpr = vector(mode = "list", length = nSets)

# transpose and add each subset
multiExpr[[1]] = list(data = as.data.frame(t(inputTable_Y05[])))
names(multiExpr[[1]]$data) = rownames(inputTable_Y05)
rownames(multiExpr[[1]]$data) = names(inputTable_Y05)
multiExpr[[2]] = list(data = as.data.frame(t(inputTable_Y023[])))
names(multiExpr[[2]]$data) = rownames(inputTable_Y023)
rownames(multiExpr[[2]]$data) = names(inputTable_Y023)
multiExpr[[3]] = list(data = as.data.frame(t(inputTable_R2[])))
names(multiExpr[[3]]$data) = rownames(inputTable_R2)
rownames(multiExpr[[3]]$data) = names(inputTable_R2)
multiExpr[[4]] = list(data = as.data.frame(t(inputTable_E05[])))
names(multiExpr[[4]]$data) = rownames(inputTable_E05)
rownames(multiExpr[[4]]$data) = names(inputTable_E05)

# Check that the data has the correct format for many functions operating on multiple sets:
exprSize = checkSets(multiExpr)
exprSize

#Check the genes across all samples
gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
gsg$allOK

if (!gsg$allOK){
  # Print information about the removed genes
  if (sum(!gsg$goodGenes) > 0){
    printFlush(paste("Removing genes:", 
                     paste(names(multiExpr[[1]]$data)[!gsg$goodGenes],
                           collapse = ", ")))
  }
  for (set in 1:exprSize$nSets){
    if (sum(!gsg$goodSamples[[set]])){
      printFlush(paste("In set", 
                       setLabels[set], 
                       "removing samples",
                       paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], 
                             collapse = ", ")))
    # Remove the offending genes and samples
    multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], 
                                                  gsg$goodGenes]
    }
  }
  # Update exprSize
  exprSize = checkSets(multiExpr)
}

#  cluster the samples on their Euclidean distance, separately in each set
sampleTrees = list()
for (set in 1:nSets){
  sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}

# plot dendrogrems
pdf(file = "SampleClustering.pdf", width = 12, height = 12)
par(mfrow=c(4,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets){
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7)
}
dev.off()

# load in the trait data
allTraits = read.csv("/Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_WGCNA_Olympics.csv");
dim(allTraits)
names(allTraits)

# Form a multi-set structure that will hold the clinical traits.
Traits = vector(mode="list", length = nSets)
for (set in 1:nSets){
  setSamples = rownames(multiExpr[[set]]$data)
  traitRows = match(setSamples, allTraits[,1])
  Traits[[set]] = list(data = allTraits[traitRows, -1])
  rownames(Traits[[set]]$data) = allTraits[traitRows, 1]
}
collectGarbage()

# define data set dimensions
nGenes = exprSize$nGenes
nSamples = exprSize$nSamples

#  save the relevant data for use in the subsequent analysis
save(multiExpr, 
     Traits, 
     nGenes, 
     nSamples, 
     setLabels, 
     shortLabels, 
     exprSize,
     file = "Consensus-dataInput.RData")
