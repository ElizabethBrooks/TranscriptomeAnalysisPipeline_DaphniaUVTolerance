#!/usr/bin/env Rscript

# R script to create euler diagrams

#install.packages("eulerr")

# turn off scientific notation
options(scipen = 999)

# load librarys
library(topGO)
library(eulerr)
library(rcartocolor)
library(stringr)

# Plotting Palettes
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
# https://github.com/Nowosad/rcartocolor
plotColors <- carto_pal(12, "Safe")
plotColorSubset <- c(plotColors[3], plotColors[2], plotColors[1])

# read in BP GO term enrichment results
treatmentGO <- read.csv("/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Genotypes/GOAnalysis/treatment_BP_GO_terms.csv")
toleranceGO <- read.csv("/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Genotypes/GOAnalysis/tolerance_BP_GO_terms.csv")
interactionGO <- read.csv("/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Genotypes/GOAnalysis/interaction_BP_GO_terms.csv")

# keep only sig
treatmentGOSig <- na.omit(toleranceGO[toleranceGO$weightFisher < 0.05,])
toleranceGOSig <- na.omit(interactionGO[interactionGO$weightFisher < 0.05,])
interactionGOSig <- na.omit(treatmentGO[treatmentGO$weightFisher < 0.05,])

# convert Annotated to numeric
treatmentGOSig$Size <- as.numeric(treatmentGOSig$Annotated)
toleranceGOSig$Size <- as.numeric(toleranceGOSig$Annotated)
interactionGOSig$Size <- as.numeric(interactionGOSig$Annotated)

# keep only small
treatmentGOSmall <- treatmentGOSig[treatmentGOSig$Size < 100,]
toleranceGOSmall <- toleranceGOSig[toleranceGOSig$Size < 100,]
interactionGOSmall <- interactionGOSig[interactionGOSig$Size < 100,]

# import DEGs
treatmentTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Genotypes/glmQLF_2WayANOVA_treatment_topTags_LFC1.2.csv", row.names="gene")
toleranceTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Genotypes/glmQLF_2WayANOVA_tolerance_topTags_LFC1.2.csv", row.names="gene")
interactionTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Genotypes/glmQLF_2WayANOVA_interaction_topTags_LFC1.2.csv", row.names="gene")

# keep only sig
treatmentSig <- treatmentTable[treatmentTable$FDR < 0.05,]
toleranceSig <- toleranceTable[toleranceTable$FDR < 0.05,]
interactionSig <- interactionTable[interactionTable$FDR < 0.05,]

# setup gene sets
geneSet_treatment <- rownames(treatmentSig)
geneSet_tolerance <- rownames(toleranceSig)
geneSet_interaction <- rownames(interactionSig)

# DE list
glm_list_venn <- list(treatment = geneSet_treatment, 
                      tolerance = geneSet_tolerance,
                      interaction = geneSet_interaction)

# DE euler diagram
euler_plot <- euler(glm_list_venn)
plot(euler_plot, fills = plotColorSubset)

# GO annotations

## GO enrichment
# read in GO mappings
GOmaps <- readMappings(file = "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/email/geneToGO_tagged_map.txt")
# create function to return list of interesting DE genes (0 == not significant, 1 == significant)
get_interesting_DE_genes <- function(geneUniverse){
  interesting_DE_genes <- rep(0, length(geneUniverse))
  for(i in 1:length(geneUniverse)){
    if(geneUniverse[i] < 0.05){
      interesting_DE_genes[i] = 1
    }
  }
  interesting_DE_genes <- setNames(interesting_DE_genes, names(geneUniverse))
  return(interesting_DE_genes)
}

# create named list of all genes (gene universe) and p-values
# the gene universe is set to be the list of all genes contained in the gene2GO list of annotated genes
list_genes_interaction <- as.numeric(interactionTable$FDR)
list_genes_interaction <- setNames(list_genes_interaction, rownames(interactionTable))
list_genes_interaction_filtered <- list_genes_interaction[names(list_genes_interaction) %in% names(GOmaps)]
list_genes_treatment <- as.numeric(interactionTable$FDR)
list_genes_treatment <- setNames(list_genes_treatment, rownames(interactionTable))
list_genes_treatment_filtered <- list_genes_treatment[names(list_genes_treatment) %in% names(GOmaps)]
list_genes_tolerance <- as.numeric(interactionTable$FDR)
list_genes_tolerance <- setNames(list_genes_tolerance, rownames(interactionTable))
list_genes_tolerance_filtered <- list_genes_tolerance[names(list_genes_tolerance) %in% names(GOmaps)]

# create topGOdata objects for enrichment analysis
BP_GO_data_interaction <- new('topGOdata', ontology = 'BP', allGenes = list_genes_interaction_filtered, 
                              geneSel = get_interesting_DE_genes, nodeSize = 10, annot = annFUN.gene2GO, 
                              gene2GO = GOmaps)
BP_GO_data_treatment <- new('topGOdata', ontology = 'BP', allGenes = list_genes_treatment_filtered, 
                            geneSel = get_interesting_DE_genes, nodeSize = 10, annot = annFUN.gene2GO, 
                            gene2GO = GOmaps)
BP_GO_data_tolerance <- new('topGOdata', ontology = 'BP', allGenes = list_genes_tolerance_filtered, 
                            geneSel = get_interesting_DE_genes, nodeSize = 10, annot = annFUN.gene2GO, 
                            gene2GO = GOmaps)

# retrieve geneIDs associated with all GO terms
# https://support.bioconductor.org/p/29775/
allGO_interaction = genesInTerm(BP_GO_data_interaction)
allGO_treatment = genesInTerm(BP_GO_data_treatment)
allGO_tolerance = genesInTerm(BP_GO_data_tolerance)

# combine GO lists
allGO <- modifyList(allGO_interaction, allGO_treatment)
allGO <- modifyList(allGO, allGO_tolerance)

# loop over each GO term ID and add to the list
# then create euler diagrams
# treatment
treatmentTermList <- list(treatment = geneSet_treatment, 
                          tolerance = geneSet_tolerance,
                          interaction = geneSet_interaction)
for(i in 1:length(treatmentGOSmall)){
  index <- i+3
  treatmentTermList[[index]] <- unlist(allGO[treatmentGOSmall$GO.ID[i]])
  names(treatmentTermList)[[index]] <- treatmentGOSmall$GO.ID[i]
}
euler_plot_treatmentGO <- euler(treatmentTermList)
plot(euler_plot_treatmentGO)
# tolerance
toleranceTermList <- list(treatment = geneSet_treatment, 
                          tolerance = geneSet_tolerance,
                          interaction = geneSet_interaction)
for(i in 1:length(toleranceGOSmall)){
  index <- i+length(treatmentGOSmall)+3
  toleranceTermList[[index]] <- unlist(allGO[toleranceGOSmall$GO.ID[i]])
  names(toleranceTermList)[[index]] <- toleranceGOSmall$GO.ID[i]
}
euler_plot_toleranceGO <- euler(toleranceTermList)
plot(euler_plot_toleranceGO)
# interaction
interactionTermList <- list(treatment = geneSet_treatment, 
                          tolerance = geneSet_tolerance,
                          interaction = geneSet_interaction)
for(i in 1:length(interactionGOSmall)){
  index <- i+length(treatmentGOSmall)+length(toleranceGOSmall)+3
  interactionTermList[[index]] <- unlist(allGO[interactionGOSmall$GO.ID[i]])
  names(interactionTermList)[[index]] <- interactionGOSmall$GO.ID[i]
}
euler_plot_interactionGO <- euler(interactionTermList)
plot(euler_plot_interactionGO)
# combined
termList <- modifyList(treatmentTermList, toleranceTermList)
termList <- modifyList(interactionTermList, termList)
euler_plot_combinedGO <- euler(termList)
plot(euler_plot_combinedGO)
