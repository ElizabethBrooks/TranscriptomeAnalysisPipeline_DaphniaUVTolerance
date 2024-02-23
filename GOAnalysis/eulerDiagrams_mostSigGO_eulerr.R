#!/usr/bin/env Rscript

# R script to create euler diagrams

#install.packages("eulerr")

# turn off scientific notation
options(scipen = 999)

# set the working directory
workingDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Genotypes"
setwd(workingDir)

# load librarys
library(topGO)
library(eulerr)
library(rcartocolor)
library(stringr)

# euler settings
eulerr_options(labels = list (fontsize = 30), quantities = list (fontsize = 30))

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
treatmentGOSig <- na.omit(treatmentGO[treatmentGO$weightFisher < 0.001,])
toleranceGOSig <- na.omit(toleranceGO[toleranceGO$weightFisher < 0.001,])
interactionGOSig <- na.omit(interactionGO[interactionGO$weightFisher < 0.001,])

# convert Annotated to numeric
treatmentGOSig$Size <- as.numeric(treatmentGOSig$Annotated)
toleranceGOSig$Size <- as.numeric(toleranceGOSig$Annotated)
interactionGOSig$Size <- as.numeric(interactionGOSig$Annotated)

# keep only small
treatmentGOSmall <- treatmentGOSig#[treatmentGOSig$Size < 75,]
toleranceGOSmall <- toleranceGOSig#[toleranceGOSig$Size < 75,]
interactionGOSmall <- interactionGOSig#[interactionGOSig$Size < 75,]

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
glm_list_venn <- list(Treatment = geneSet_treatment, 
                      Tolerance = geneSet_tolerance,
                      Interaction = geneSet_interaction)

# DE euler diagram
euler_plot <- euler(glm_list_venn)
jpeg("sigDEGS_euler.jpg")
plot(euler_plot, fills = plotColorSubset, quantities = list(type = "counts"))
dev.off()


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
list_genes_treatment <- as.numeric(treatmentTable$FDR)
list_genes_treatment <- setNames(list_genes_treatment, rownames(treatmentTable))
list_genes_treatment_filtered <- list_genes_treatment[names(list_genes_treatment) %in% names(GOmaps)]
list_genes_tolerance <- as.numeric(toleranceTable$FDR)
list_genes_tolerance <- setNames(list_genes_tolerance, rownames(toleranceTable))
list_genes_tolerance_filtered <- list_genes_tolerance[names(list_genes_tolerance) %in% names(GOmaps)]
list_genes_interaction <- as.numeric(interactionTable$FDR)
list_genes_interaction <- setNames(list_genes_interaction, rownames(interactionTable))
list_genes_interaction_filtered <- list_genes_interaction[names(list_genes_interaction) %in% names(GOmaps)]

# create topGOdata objects for enrichment analysis
BP_GO_data_treatment <- new('topGOdata', ontology = 'BP', allGenes = list_genes_treatment_filtered, 
                            geneSel = get_interesting_DE_genes, nodeSize = 10, annot = annFUN.gene2GO, 
                            gene2GO = GOmaps)
BP_GO_data_tolerance <- new('topGOdata', ontology = 'BP', allGenes = list_genes_tolerance_filtered, 
                            geneSel = get_interesting_DE_genes, nodeSize = 10, annot = annFUN.gene2GO, 
                            gene2GO = GOmaps)
BP_GO_data_interaction <- new('topGOdata', ontology = 'BP', allGenes = list_genes_interaction_filtered, 
                              geneSel = get_interesting_DE_genes, nodeSize = 10, annot = annFUN.gene2GO, 
                              gene2GO = GOmaps)

# retrieve geneIDs associated with all GO terms
# https://support.bioconductor.org/p/29775/
allGO_treatment = genesInTerm(BP_GO_data_treatment)
allGO_tolerance = genesInTerm(BP_GO_data_tolerance)
allGO_interaction = genesInTerm(BP_GO_data_interaction)

# combine GO lists
allGO <- modifyList(allGO_interaction, allGO_treatment)
allGO <- modifyList(allGO, allGO_tolerance)

# loop over each GO term ID and add to the list
# then create euler diagrams
# treatment
treatmentTermList <- list()
for(i in 1:length(treatmentGOSmall)){
  treatmentTermList[[i]] <- unlist(allGO[treatmentGOSmall$GO.ID[i]])
}
names(treatmentTermList) <- treatmentGOSmall$GO.ID
treatmentTermList <- modifyList(glm_list_venn, treatmentTermList)
euler_plot_treatmentGO <- euler(treatmentTermList)
jpeg("GOAnalysis/sigBP_treatment_euler.jpg")
plot(euler_plot_treatmentGO)
dev.off()
# tolerance
toleranceTermList <- list()
for(i in 1:length(toleranceGOSmall)){
  toleranceTermList[[i]] <- unlist(allGO[toleranceGOSmall$GO.ID[i]])
}
names(toleranceTermList) <- toleranceGOSmall$GO.ID
toleranceTermList <- modifyList(glm_list_venn, toleranceTermList)
euler_plot_toleranceGO <- euler(toleranceTermList)
jpeg("GOAnalysis/sigBP_tolerance_euler.jpg")
plot(euler_plot_toleranceGO)
dev.off()
# interaction
interactionTermList <- list()
for(i in 1:length(interactionGOSmall)){
  interactionTermList[[i]] <- unlist(allGO[interactionGOSmall$GO.ID[i]])
}
names(interactionTermList) <- interactionGOSmall$GO.ID
interactionTermList <- modifyList(glm_list_venn, interactionTermList)
euler_plot_interactionGO <- euler(interactionTermList)
jpeg("GOAnalysis/sigBP_interaction_euler.jpg")
plot(euler_plot_interactionGO)
dev.off()
# combined
termList <- modifyList(glm_list_venn, toleranceTermList)
termList <- modifyList(treatmentTermList, termList)
termList <- modifyList(interactionTermList, termList)
euler_plot_combinedGO <- euler(termList)
jpeg("GOAnalysis/sigBP_combined_euler.jpg")
plot(euler_plot_combinedGO)
dev.off()

# geneIDs with repair or other interesting GO terms
# The repair pathway GO terms we specifically explored in the analysis were DNA repair (GO:0006281), mismatch repair (MMR; GO:0006298), base excision repair (BER; GO:0006284), homologous recombination (HR; GO:0035825), nucleotide excision repair (NER; GO:0006289), intrastrand crosslink repair (ICL repair; GO:0036297), double strand break repair (DSBR; GO:0006302), single strand break repair (SSBR; GO:0000012).
DNAR <- unlist(allGO["GO:0006281"]) # not sig
MMR <- unlist(allGO["GO:0006298"]) # not sig
BER <- unlist(allGO["GO:0006284"]) # not sig
HR <- unlist(allGO["GO:0035825"]) # not sig
NER <- unlist(allGO["GO:0006289"]) # not sig
ICLR <- unlist(allGO["GO:0036297"]) # not sig
DSBR <- unlist(allGO["GO:0006302"]) # not sig
SSBR <- unlist(allGO["GO:0000012"]) # not sig
# Interestingly, the cellular response to oxidative stress (GO:0034599) was significantly enriched in the treatment effect. The heme biosynthetic process (GO:0006783) and alpha-amino acid biosynthetic process (GO:1901607) were enriched for both the treatment and interaction effects (S. Figure 2).
#oxidativeStress <- unlist(allGO["GO:0034599"]) # treatment sig
#hemeProcess <- unlist(allGO["GO:0006783"]) # interaction and treatment sig
#alphaProcess <- unlist(allGO["GO:1901607"]) # interaction and treatment sig
# The radiation response terms included response to radiation	(GO:0009314), response to ionizing radiation (GO:0010212), response to gamma radiation (GO:0010332), cellular response to radiation (GO:0071478), cellular response to gamma radiation (GO:0071480), cellular response to ionizing radiation (GO:0071479), regulation of response to gamma radiation (GO:2001228), regulation of cellular response to gamma radiation (GO:1905843), positive regulation of response to gamma radiation (GO:2001230), negative regulation of response to gamma radiation (GO:2001229), 
# response to X-ray (GO:0010165), positive regulation of cellular response to gamma radiation	(GO:1905845), negative regulation of cellular response to gamma radiation (GO:1905844), cellular response to X-ray (GO:0071481), regulation of cellular response to X-ray (GO:2000683), negative regulation of cellular response to X-ray (GO:2000684), positive regulation of cellular response to X-ray (GO:2000685), 
# phototransduction UV (GO:0007604), response to UV (GO:0009411), response to UV-A (GO:0070141), detection of UV (GO:0009589), cellular response to UV (GO:0034644), cellular response to UV-A (GO:0071492), response to UV-B	4 (GO:0010224), cellular response to UV-B (GO:0071493), response to UV-C (GO:0010225), 
# regulation of mRNA stability involved in cellular response to UV (GO:1902629), cellular response to UV-C (GO:0071494), regulation of translation involved in cellular response to UV (GO:1904803).
#radiation <- unlist(allGO["GO:0009314"]) # not sig
#cell <- unlist(allGO["GO:0071478"]) # not sig
# The stress terms included response to stress (GO:0006950), response to oxidative stress (GO:0006979), cellular response to oxidative stress (GO:0034599), 
# positive regulation of translation in response to oxidative stress (GO:0032939), cellular response to reactive oxygen (GO:0034614), negative regulation of cellular response to oxidative stress (GO:1900408), positive regulation of cellular response to oxidative stress (GO:1900409), 
# negative regulation of translation in response to oxidative stress (GO:0032938), regulation of translation in response to oxidative stress (GO:0043556), Mcs4 RR-MAPKKK complex (GO:1990315), cellular response to hydroperoxide (GO:0071447), 
# oxidative stress-induced premature senescence (GO:0090403), symbiont defense to host-produced reactive oxygen species (GO:0052164), regulation of cellular response to oxidative stress (GO:1900407).
stress <- unlist(allGO["GO:0006950"]) # interaction sig
oxidative <- unlist(allGO["GO:0006979"]) # interaction sig
cellOxidative <- unlist(allGO["GO:0034599"]) # treatment sig
#cellReactiveOxy <- unlist(allGO["GO:0034614"]) # not sig

# euler diagram with interesting repair GO terms
glm_list_venn_repairGO <-list(tolerance = geneSet_tolerance,
                              treatment = geneSet_treatment, 
                              interaction = geneSet_interaction,
                              DNAR = DNAR,
                              MMR = MMR,
                              BER = BER,
                              HR = HR,
                              NER = NER,
                              ICLR = ICLR,
                              DSBR = DSBR)
#SSBR = SSBR)
euler_plot_repairGO <- euler(glm_list_venn_repairGO)
jpeg("GOAnalysis/notSigGO_DNARepair_euler.jpg")
plot(euler_plot_repairGO)
dev.off()

# euler diagram with significant stress GO terms
glm_list_venn_stressGO <-list(tolerance = geneSet_tolerance,
                              treatment = geneSet_treatment, 
                              interaction = geneSet_interaction,
                              "GO:0006950" = stress,
                              "GO:0006979" = oxidative,
                              "GO:0034599" = cellOxidative)
#cellReactiveOxy = cellReactiveOxy)
euler_plot_stressGO <- euler(glm_list_venn_stressGO)#, shape = "ellipse")
jpeg("GOAnalysis/sigBP_stress_euler.jpg")
plot(euler_plot_stressGO, quantities = list(type = c("counts")))#, fills = plotColors[1:6])
dev.off()
