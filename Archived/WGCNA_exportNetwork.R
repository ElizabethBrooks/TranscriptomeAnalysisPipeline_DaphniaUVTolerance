#Set working directory
#workingDir = args[1];
workingDir="/home/mae/Documents/RNASeq_Workshop_ND/WGCNA_PA42_v4.1/allGenes"
setwd(workingDir); 

# Load the WGCNA package
library(WGCNA)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

# Load the expression and trait data saved in the first part
lnames = load(file = "PA42_v4.1_dataInput.RData");
#The variable lnames contains the names of loaded variables.
#lnames

# Load network data saved in the second part.
lnames = load(file = "PA42_v4.1_networkConstruction_auto_threshold8_signed.RData");
#lnames


# Recalculate topological overlap
#TOM = TOMsimilarityFromExpr(datExpr, power = 7);
# Read in the annotation file
#annot = read.csv(file = "GeneAnnotation.csv");
# Select module
#module = "brown";
# Select module probes
#probes = names(datExpr)
#inModule = (moduleColors==module);
#modProbes = probes[inModule];
# Select the corresponding Topological Overlap
#modTOM = TOM[inModule, inModule];
#dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into an edge list file VisANT can read
#vis = exportNetworkToVisANT(modTOM,
#                            file = paste("VisANTInput-", module, ".txt", sep=""),
#                            weighted = TRUE,
#                            threshold = 0,
#                            probeToGene = data.frame(annot$substanceBXH, annot$gene_symbol) )


#Prepare data for VisANT
#nTop = 30;
#IMConn = softConnectivity(datExpr[, modProbes]);
#top = (rank(-IMConn) <= nTop)
#vis = exportNetworkToVisANT(modTOM[top, top],
#                            file = paste("VisANTInput-", module, "-top30.txt", sep=""),
#                            weighted = TRUE,
#                            threshold = 0,
#                            probeToGene = data.frame(annot$substanceBXH, annot$gene_symbol) )


#Prepare data for cytoscape
# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power = 8);
# Read in the annotation file
annot = read.csv(file = "geneAnnotations.csv");


# Select modules
modules = c("pink");
# Select module probes
probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modGenes = annot$entrezID[match(modProbes, annot$gene)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
#                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.1,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])


# Select modules
modules = c("brown");
# Select module probes
probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modGenes = annot$entrezID[match(modProbes, annot$gene)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
#                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.1,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])


# Read in the annotation file
annot = read.csv(file = "geneAnnotations_entrezSubset.csv");
# Select modules
modules = c("lightgreen");
# Select module probes
probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modGenes = annot$entrezID[match(modProbes, annot$gene)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
#                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.2,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])


# Select modules
modules = c("lightcyan");
# Select module probes
probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modGenes = annot$entrezID[match(modProbes, annot$gene)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
#                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.2,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])
