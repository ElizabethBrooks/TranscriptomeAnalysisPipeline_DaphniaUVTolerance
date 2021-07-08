#Set working directory
#workingDir = args[1];
workingDir="/home/mae/Documents/RNASeq_Workshop_ND/WGCNA_PA42_v4.1/effectSubsets"
setwd(workingDir); 

# Load the WGCNA package
library(WGCNA)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments. 
# See note above.
#enableWGCNAThreads()


# Load the data saved in the first part
lnamesInter = load(file = "PA42_v4.1_dataInputInter.RData");

powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Choose a set of interaction soft-thresholding powers
# Call the network topology analysis function
sft = pickSoftThreshold(datExprInter, powerVector = powers, verbose = 5)
# Plot the results:
jpeg("thresholdingPowersInter.jpg", width = 960, height = 480)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


#Construct the network in blocks of the specified size
netInter = blockwiseModules(datExprInter, power = 8,
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "PA42TOMInter_threshold8_signed", 
                       verbose = 3, maxBlockSize = 15000)


# open a graphics window for interaction
jpeg("clusterDendrogramInter_threshold8_signed.jpg", width = 960, height = 960)
# Convert labels to colors for plotting
mergedColors = labels2colors(netInter$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(netInter$dendrograms[[1]], mergedColors[netInter$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


#Add interaction module colors
moduleLabels = netInter$colors
moduleColors = labels2colors(netInter$colors)
MEs = netInter$MEs;
geneTree = netInter$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "PA42_v4.1_networkConstructionInter_auto_threshold8_signed.RData")
