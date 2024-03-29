#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)
#Test if there is one input argument

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

# Load network data saved in the second part.
lnames2 = load(file = "PA42_v4.1_networkConstructionInter_auto_threshold8_signed.RData");


# Recalculate topological overlap
#TOM = TOMsimilarityFromExpr(datExprInter, power = 7);
# Read in the annotation file
#annot = read.csv(file = "GeneAnnotation.csv");
# Select module
#module = "brown";
# Select module probes
#probes = names(datExprInter)
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
#IMConn = softConnectivity(datExprInter[, modProbes]);
#top = (rank(-IMConn) <= nTop)
#vis = exportNetworkToVisANT(modTOM[top, top],
#                            file = paste("VisANTInput-", module, "-top30.txt", sep=""),
#                            weighted = TRUE,
#                            threshold = 0,
#                            probeToGene = data.frame(annot$substanceBXH, annot$gene_symbol) )


#Prepare data for cytoscape
# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExprInter, power = 8);
# Read in the annotation file
annot = read.csv(file = "geneAnnotations.csv");


# Select modules
modules = c("turquoise");
# Select module probes
probes = names(datExprInter)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modGenes = annot$entrezID[match(modProbes, annot$gene)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputInter-edges-", paste(modules, collapse="-"), ".txt", sep=""),
#                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.1,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])


# Select modules
modules = c("green");
# Select module probes
probes = names(datExprInter)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modGenes = annot$entrezID[match(modProbes, annot$gene)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputInter-edges-", paste(modules, collapse="-"), ".txt", sep=""),
#                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.1,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])


# Read in the annotation file
annot = read.csv(file = "geneAnnotations_entrezSubset.csv");
# Select modules
modules = c("blue");
# Select module probes
probes = names(datExprInter)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modGenes = annot$entrezID[match(modProbes, annot$gene)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputInter-edges-", paste(modules, collapse="-"), ".txt", sep=""),
#                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.2,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])


# Select modules
modules = c("pink");
# Select module probes
probes = names(datExprInter)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modGenes = annot$entrezID[match(modProbes, annot$gene)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInputInter-edges-", paste(modules, collapse="-"), ".txt", sep=""),
#                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.2,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])
