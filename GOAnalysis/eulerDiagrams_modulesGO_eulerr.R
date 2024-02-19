#!/usr/bin/env Rscript

# R script to create euler diagrams

#install.packages("eulerr")

# turn off scientific notation
options(scipen = 999)

# set the working directory
workingDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Genotypes"
setwd(workingDir)

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
#jpeg("sigDEGS_euler.jpg")
plot(euler_plot, fills = plotColorSubset, quantities = list(type = c("percent", "counts")))
#dev.off()

# import modules GO data
modulesNames <- list.files(path="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Genotypes/GOAnalysis_OLYM_30/", pattern = "_BP_sigGO_terms.csv$")
modulesNames <- str_remove(modulesNames, "_BP_sigGO_terms.csv")
modulesFiles <- list.files(path="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Genotypes/GOAnalysis_OLYM_30/", pattern = "_BP_sigGO_terms.csv$", full.names = TRUE)
modulesGO <- lapply(modulesFiles, read.csv)
names(modulesGO) <- modulesNames

# retrieve subset tag
set <- "OLYM"

# set the minimum module size
minModSize <- "30"

# retrieve WGCNA directory
inDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Genotypes"

# set the full subset tag name
tag <- paste(set, minModSize, sep="_")

# Load the expression and trait data saved in the first part
importFile <- paste(set, "dataInput.RData", sep="-")
importFile <- paste(inDir, importFile, sep="/")
lnames1 = load(file = importFile)

# Load network data saved in the second part
importFile <- paste(tag, "networkConstruction-stepByStep.RData", sep="-")
importFile <- paste(inDir, importFile, sep="/")
lnames2 = load(file = importFile)

# create list of module colors mapped to numbers
numMods <- length(unique(moduleColors))
colorTable <- data.frame(
  color = unique(moduleColors),
  number = seq(from = 1, to = numMods, by = 1)
)

# initialize module data frame
resultsTable <- data.frame(
  gene = character(),
  color = character(),
  number = numeric()
)

# match gene IDs with module colors
for(i in 1:numMods){
  gene <- names(datExpr)[moduleColors==colorTable[i,1]]
  color <- rep(colorTable[i,1], length(gene))
  number <- rep(colorTable[i,2], length(gene))
  moduleData <- cbind(gene, color, number)
  resultsTable <- rbind(resultsTable, moduleData)
}

# setup data frame of module sizes
modSizes <- data.frame(
  number = numeric(),
  color = character(),
  size = numeric()
)

# retrieve module sizes
for(i in 1:numMods){
  number <- head(resultsTable[resultsTable$number == i,3],1)
  color <- head(resultsTable[resultsTable$number == i,2],1)
  size <- nrow(resultsTable[resultsTable$number == i,])
  moduleData <- cbind(number, color, size)
  modSizes <-rbind(modSizes, moduleData)
}

# read in GO mappings
GOmaps <- readMappings(file = "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/email/geneToGO_tagged_map.txt")

# create named list of all genes (gene universe) and p-values. The gene universe is set to be
# the list of all genes contained in the gene2GO list of annotated genes.
list_genes <- as.numeric(resultsTable$number)
list_genes <- setNames(list_genes, resultsTable$gene)
list_genes_filtered <- list_genes[names(list_genes) %in% names(GOmaps)]

# retrieve genes in modules with significant DNA repair BP terms
#steelblueGenes <- resultsTable[resultsTable$color == "steelblue",]$gene
#floralwhiteGenes <- resultsTable[resultsTable$color == "floralwhite",]$gene
#salmon4Genes <- resultsTable[resultsTable$color == "salmon4",]$gene

# GO annotations

# function to remove null lists
# https://stackoverflow.com/questions/26539441/remove-null-elements-from-list-of-lists
## A helper function that tests whether an object is either NULL _or_ 
## a list of NULLs
is.NullOb <- function(x) is.null(x) | all(sapply(x, is.null))
## Recursively step down into list, removing all such objects 
rmNullObs <- function(x) {
  x <- Filter(Negate(is.NullOb), x)
  lapply(x, function(x) if (is.list(x)) rmNullObs(x) else x)
}

# create list for lists of GO terms
modulesGOList <- split(colorTable, seq(nrow(colorTable)))
names(modulesGOList) <- colorTable$color

## GO enrichment
#Loop through each module
for(j in 1:numMods){
  # status message
  print(paste("Processing", colorTable[j,1]))
  
  #Create function to return list of interesting DE genes (0 == not significant, 1 == significant)
  get_interesting_DE_genes <- function(geneUniverse){
    interesting_DE_genes <- rep(0, length(geneUniverse))
    for(i in 1:length(geneUniverse)){
      if(geneUniverse[i] == colorTable[j,2]){
        interesting_DE_genes[i] = 1
      }
    }
    interesting_DE_genes <- setNames(interesting_DE_genes, names(geneUniverse))
    return(interesting_DE_genes)
  }
  
  #Create topGOdata objects for enrichment analysis (1 for each ontology)
  BP_GO_data <- new('topGOdata', ontology = 'BP', allGenes = list_genes_filtered, 
                    geneSel = get_interesting_DE_genes, nodeSize = 10, annot = annFUN.gene2GO, 
                    gene2GO = GOmaps)

  # retrieve geneIDs associated with all GO terms
  # https://support.bioconductor.org/p/29775/
  allGO = genesInTerm(BP_GO_data)
  
  # geneIDs with repair or other interesting GO terms
  # The repair pathway GO terms we specifically explored in the analysis were DNA repair (GO:0006281), mismatch repair (MMR; GO:0006298), base excision repair (BER; GO:0006284), homologous recombination (HR; GO:0035825), nucleotide excision repair (NER; GO:0006289), intrastrand crosslink repair (ICL repair; GO:0036297), double strand break repair (DSBR; GO:0006302), single strand break repair (SSBR; GO:0000012).
  repairTerms <- list(DNAR = unlist(allGO["GO:0006281"]),
                        MMR = unlist(allGO["GO:0006298"]),
                        BER = unlist(allGO["GO:0006284"]),
                        HR = unlist(allGO["GO:0035825"]),
                        NER = unlist(allGO["GO:0006289"]),
                        ICLR = unlist(allGO["GO:0036297"]),
                        DSBR = unlist(allGO["GO:0006302"]),
                        SSBR = unlist(allGO["GO:0000012"])
  )
  
  # The radiation response terms included response to radiation (GO:0009314), cellular response to radiation (GO:0071478), phototransduction UV (GO:0007604),
  # response to UV (GO:0009411), response to UV-A (GO:0070141), detection of UV (GO:0009589), cellular response to UV (GO:0034644), cellular response to UV-A (GO:0071492), 
  # response to UV-B4 (GO:0010224), cellular response to UV-B (GO:0071493), response to UV-C (GO:0010225), regulation of mRNA stability involved in cellular response to UV (GO:1902629), 
  # cellular response to UV-C (GO:0071494), regulation of translation involved in cellular response to UV (GO:1904803).
  radiationTerms <- list(radiation = unlist(allGO["GO:0009314"]),
                          ionizing = unlist(allGO["GO:0010212"]),
                          gamma = unlist(allGO["GO:0010332"]),
                          cell = unlist(allGO["GO:0071478"]),
                          cellGamma = unlist(allGO["GO:0071480"]),
                          cellIonizing = unlist(allGO["GO:0071479"]),
                          regGamma = unlist(allGO["GO:2001228"]),
                          regCellGamma = unlist(allGO["GO:1905843"]),
                          xRay = unlist(allGO["GO:0010165"]),
                          cellXRay = unlist(allGO["GO:0071481"]),
                          regCellXRay = unlist(allGO["GO:2000683"]),
                          photoUV = unlist(allGO["GO:0007604"]),
                          UV = unlist(allGO["GO:0009411"]),
                          UVA = unlist(allGO["GO:0070141"]),
                          detectUV = unlist(allGO["GO:0009589"]),
                          cellUV = unlist(allGO["GO:0034644"]),
                          cellUVA = unlist(allGO["GO:0071492"]),
                          UVB4 = unlist(allGO["GO:0010224"]),
                          cellUVB = unlist(allGO["GO:0071493"]),
                          UVC = unlist(allGO["GO:0010225"]),
                          mRNACellUV = unlist(allGO["GO:1902629"]),
                          cellUVC = unlist(allGO["GO:0071494"]),
                          transCellUV = unlist(allGO["GO:1904803"])
  )

  # The stress terms included response to stress (GO:0006950), response to oxidative stress (GO:0006979), cellular response to oxidative stress (GO:0034599), 
  # cellular response to reactive oxygen (GO:0034614), regulation of translation in response to oxidative stress (GO:0043556), regulation of cellular response to oxidative stress (GO:1900407).
  stressTerms <- list(stress = unlist(allGO["GO:0006950"]),
                        oxidative = unlist(allGO["GO:0006979"]),
                        cellOxidative = unlist(allGO["GO:0034599"]),
                        cellReactiveOxy = unlist(allGO["GO:0034614"]),
                        regTransOxidative = unlist(allGO["GO:0043556"]),
                        MAPKKK = unlist(allGO["GO:1990315"]),
                        hydroperoxide = unlist(allGO["GO:0071447"]),
                        senescence = unlist(allGO["GO:0090403"]),
                        symbiont = unlist(allGO["GO:0052164"]),
                        regCellOxidative = unlist(allGO["GO:1900407"])
  )
  
  # add list of lists of GO terms
  modulesGOList[[colorTable[j,1]]] <- list(repairTerms = repairTerms,
                                           radiationTerms = radiationTerms,
                                           stressTerms = stressTerms)
}

# remove nulls
modulesGOList_cleaned <- rmNullObs(modulesGOList)

#Loop through each module
for(j in 1:1){
  # status message
  print(paste("Plotting", colorTable[j,1]))
  
  # create plot list
  repairModulesGO <- list()
  repairModulesGO <-modulesGOList_cleaned[[colorTable[j,1]]][["repairTerms"]]
  repairModulesGO[[length(repairModulesGO)+1]] <- c(resultsTable[resultsTable$color == colorTable[j,1],]$gene)
  names(repairModulesGO[[length(repairModulesGO)]]) <- rep(colorTable[j,1], length(repairModulesGO[[length(repairModulesGO)]]))
  names(repairModulesGO) <- c(names(repairModulesGO)[1:length(repairModulesGO)-1], colorTable[j,1])
  
  # euler diagram with interesting repair modules GO terms
  euler_plot_repairModulesGO <- euler(repairModulesGO)
  fileOut <- paste(colorTable[j,1], "repair_euler.jpg", sep = "_")
  fileOut <- paste("GOAnalysis_OLYM_30", fileOut, sep = "/")
  jpeg(filename = fileOut)
  plot(euler_plot_repairModulesGO)
  dev.off()
  
  # create plot list
  radiationModulesGO <- list()
  radiationModulesGO <-modulesGOList_cleaned[[colorTable[j,1]]][["radiationTerms"]]
  radiationModulesGO[[length(radiationModulesGO)+1]] <- c(resultsTable[resultsTable$color == colorTable[j,1],]$gene)
  names(radiationModulesGO[[length(radiationModulesGO)]]) <- rep(colorTable[j,1], length(radiationModulesGO[[length(radiationModulesGO)]]))
  names(radiationModulesGO) <- c(names(radiationModulesGO)[1:length(radiationModulesGO)-1], colorTable[j,1])
  
  # euler diagram with interesting repair modules GO terms
  euler_plot_radiationModulesGO <- euler(radiationModulesGO)
  fileOut <- paste(colorTable[j,1], "radiation_euler.jpg", sep = "_")
  fileOut <- paste("GOAnalysis_OLYM_30", fileOut, sep = "/")
  jpeg(filename = fileOut)
  plot(euler_plot_radiationModulesGO)
  dev.off()
  
  # create plot list
  stressModulesGO <- list()
  stressModulesGO <-modulesGOList_cleaned[[colorTable[j,1]]][["stressTerms"]]
  stressModulesGO[[length(stressModulesGO)+1]] <- c(resultsTable[resultsTable$color == colorTable[j,1],]$gene)
  names(stressModulesGO[[length(stressModulesGO)]]) <- rep(colorTable[j,1], length(stressModulesGO[[length(stressModulesGO)]]))
  names(stressModulesGO) <- c(names(stressModulesGO)[1:length(stressModulesGO)-1], colorTable[j,1])
  
  # euler diagram with interesting repair modules GO terms
  euler_plot_stressModulesGO <- euler(stressModulesGO)
  fileOut <- paste(colorTable[j,1], "stress_euler.jpg", sep = "_")
  fileOut <- paste("GOAnalysis_OLYM_30", fileOut, sep = "/")
  jpeg(filename = fileOut)
  plot(euler_plot_stressModulesGO)
  dev.off()
}
