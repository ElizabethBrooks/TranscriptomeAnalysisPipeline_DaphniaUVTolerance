#Set working directory
workingDir = args[1];
#workingDir="/home/mae/Documents/RNASeq_Workshop_ND/WGCNA_PA42_v4.1/effectSubsets"
setwd(workingDir); 

# Load the WGCNA package
library(WGCNA)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

# Load the expression and trait data saved in the first part
lnames1 = load(file = "PA42_v4.1_dataInputInter.RData");
#The variable lnames contains the names of loaded variables.
#lnames

# Load network data saved in the second part.
lnames2 = load(file = "PA42_v4.1_networkConstructionInter_auto_threshold8_signed.RData");


nGenes = ncol(datExprInter)
nSamples = nrow(datExprInter)


# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(datExprInter, power = 8);


#Generate a network heatmap for a subset of genes
nSelect = 400
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window
#sizeGrWindow(9,9)
jpeg("networkHeatmapInter_subset400.jpg", width = 960, height = 960)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing 
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^10;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
dev.off()


# Recalculate module eigengenes
MEs = moduleEigengenes(datExprInter, moduleColors)$eigengenes
# Isolate treatment from the clinical traits
treatment = as.data.frame(datTraitsInter$treatment);
names(treatment) = "treatment"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, treatment))
# Plot the relationships among the eigengenes and the trait
#sizeGrWindow(5,7.5);
jpeg("moduleEigengenesInter_treatment.jpg", width = 960, height = 960)
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)
dev.off()

# Plot the dendrogram
#sizeGrWindow(6,6);
jpeg("eigengenesDendrogramInter_treatment.jpg", width = 960, height = 960)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
dev.off()

# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
jpeg("eigengenesHeatmapInter_treatment.jpg", width = 960, height = 960)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()


# Recalculate module eigengenes
MEs = moduleEigengenes(datExprInter, moduleColors)$eigengenes
# Isolate tolerance from the clinical traits
tolerance = as.data.frame(datTraitsInter$tolerance);
names(tolerance) = "tolerance"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, tolerance))
# Plot the relationships among the eigengenes and the trait
#sizeGrWindow(5,7.5);
jpeg("moduleEigengenesInter_tolerance.jpg", width = 960, height = 960)
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)
dev.off()

# Plot the dendrogram
#sizeGrWindow(6,6);
jpeg("eigengenesDendrogramInter_tolerance.jpg", width = 960, height = 960)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
dev.off()

# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
jpeg("eigengenesHeatmapInter_tolerance.jpg", width = 960, height = 960)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()


# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
#sizeGrWindow(9,9)
jpeg("networkHeatmapInter.jpg", width = 960, height = 960)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
dev.off()
