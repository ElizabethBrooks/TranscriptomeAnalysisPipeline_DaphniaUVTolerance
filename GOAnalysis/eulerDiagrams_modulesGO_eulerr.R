#!/usr/bin/env Rscript

# R script to create euler diagrams

#install.packages("eulerr")

# turn off scientific notation
options(scipen = 999)

# set the working directory
workingDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Genotypes/GOAnalysis_OLYM_30"
setwd(workingDir)

# load librarys
library(topGO)
library(eulerr)
library(rcartocolor)
library(stringr)
library(tidyr)

# euler settings
eulerr_options(labels = list(fontsize = 30), quantities = list(fontsize = 30), fills = list (alpha = 0.75), legend = list(fontsize = 30))

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

## PUB
# DE euler diagram
euler_plot <- euler(glm_list_venn)
#jpeg("sigDEGS_euler.jpg")
plot(euler_plot, fills = plotColorSubset, quantities = list(type = "counts"), labels = NULL, legend = list(labels = c("Treatment", "Tolerance", "Interaction")))
#dev.off()

# import modules GO data
modulesNames <- list.files(path="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Genotypes/GOAnalysis_OLYM_30/", pattern = "_BP_sigGO_terms.csv.flt$")
modulesNames <- str_remove(modulesNames, "_BP_sigGO_terms.csv.flt")
modulesFiles <- list.files(path="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Genotypes/GOAnalysis_OLYM_30/", pattern = "_BP_sigGO_terms.csv.flt$", full.names = TRUE)
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

# retrieve dN dS values
#dNdSTable <- read.csv(file = args[6], row.names="geneID")
dNdSTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/SelectionTests/Pulex_Olympics_kaksResults.csv", row.names="geneID")

# remove NAs
dNdSSubset <- na.omit(dNdSTable)

# merging data frames
# https://sparkbyexamples.com/r-programming/r-join-data-frames-with-examples/#full-outer-join

# add geneID column
geneID <- row.names(dNdSSubset)
dNdSSubset <- cbind(geneID,dNdSSubset)
colnames(resultsTable)[1] ="geneID"

# full outer join data frames
resultsTable <- merge(x = resultsTable, y = dNdSSubset, 
                      by = "geneID", all=TRUE)

# set color and numner tags for NAs
resultsTable$color <- resultsTable$color %>% replace_na('None')
resultsTable$number <- resultsTable$number %>% replace_na('0')

# remove NAs
resultsTable <- na.omit(resultsTable)

# update the number of modules
numMods <- numMods + 1

# update the table of module colors
colorTable[nrow(colorTable) + 1,] <- c("None","0")

# remove Nones
resultsTable_subset <- resultsTable[!grepl("None", resultsTable$color),]


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

# combine all tables
effectTable <- rbind(interactionSig, treatmentSig, toleranceSig)

# add geneID column
effectTable$geneID <- row.names(effectTable)

# subset dN dS values to remove outliers
dNdSSubset <- dNdSSubset[dNdSSubset$dNdS < 99,]

# subset dN dS values for positively selected genes
dNdSSubset <- dNdSSubset[dNdSSubset$dNdS > 1,]

# setup gene sets
geneSet_positive <- row.names(dNdSSubset)
geneSet_DEGs <- effectTable$geneID
geneSet_Modules <- resultsTable_subset$gene

## PUB
# euler diagram of modules, DEGs, and positively selected genes
positiveSet_list <- list(DEGs = geneSet_DEGs,
                         Positive = geneSet_positive,
                         Modules = geneSet_Modules)
euler_plot_modulesDEGsPositive <- euler(positiveSet_list)
#jpeg(filename = "positive_DEGs_modules_euler.jpg")
plot(euler_plot_modulesDEGsPositive, fills = list(fill = c(plotColors)), quantities = list(type = "counts"), labels = NULL, legend = list(labels = c("DEGs", "Positive", "Modules")))
#dev.off()

# euler diagram of modules, separate DEGs, and positively selected genes
separateSet_list <- list(Tolerance = geneSet_tolerance,
                         Treatment = geneSet_treatment, 
                         Interaction = geneSet_interaction,
                         Positive = geneSet_positive,
                         Modules = geneSet_Modules)
euler_plot_modulesSeparateDEGsPositive <- euler(separateSet_list)
#jpeg(filename = "positive_separateDEGs_modules_euler.jpg")
plot(euler_plot_modulesSeparateDEGsPositive, fills = plotColors, quantities = list(type = "counts"))
#dev.off()

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
  #jpeg(filename = fileOut)
  plot(euler_plot_repairModulesGO)
  #dev.off()
  
  # create plot list
  radiationModulesGO <- list()
  radiationModulesGO <-modulesGOList_cleaned[[colorTable[j,1]]][["radiationTerms"]]
  radiationModulesGO[[length(radiationModulesGO)+1]] <- c(resultsTable[resultsTable$color == colorTable[j,1],]$gene)
  names(radiationModulesGO[[length(radiationModulesGO)]]) <- rep(colorTable[j,1], length(radiationModulesGO[[length(radiationModulesGO)]]))
  names(radiationModulesGO) <- c(names(radiationModulesGO)[1:length(radiationModulesGO)-1], colorTable[j,1])
  
  # euler diagram with interesting repair modules GO terms
  euler_plot_radiationModulesGO <- euler(radiationModulesGO)
  fileOut <- paste(colorTable[j,1], "radiation_euler.jpg", sep = "_")
  #jpeg(filename = fileOut)
  plot(euler_plot_radiationModulesGO)
  #dev.off()
  
  # create plot list
  stressModulesGO <- list()
  stressModulesGO <-modulesGOList_cleaned[[colorTable[j,1]]][["stressTerms"]]
  stressModulesGO[[length(stressModulesGO)+1]] <- c(resultsTable[resultsTable$color == colorTable[j,1],]$gene)
  names(stressModulesGO[[length(stressModulesGO)]]) <- rep(colorTable[j,1], length(stressModulesGO[[length(stressModulesGO)]]))
  names(stressModulesGO) <- c(names(stressModulesGO)[1:length(stressModulesGO)-1], colorTable[j,1])
  
  # euler diagram with interesting repair modules GO terms
  euler_plot_stressModulesGO <- euler(stressModulesGO)
  fileOut <- paste(colorTable[j,1], "stress_euler.jpg", sep = "_")
  #jpeg(filename = fileOut)
  plot(euler_plot_stressModulesGO)
  #dev.off()
}

# significant repair modules diagrams
DDR <- unlist(allGO["GO:0006974"]) # DNA damage response
DNAR <- unlist(allGO["GO:0006281"])
MMR <- unlist(allGO["GO:0006298"])
ICLR <- unlist(allGO["GO:0036297"])
BER <- unlist(allGO["GO:0006284"])
NER <- unlist(allGO["GO:0006289"])
# "None" "GO:0006281" 
#noneModulesGO <- c(resultsTable[resultsTable$color == "None",]$gene)
#names(noneModulesGO) <- rep("None", length(noneModulesGO))
# euler diagram with interesting repair GO terms
#noneModule_repair_euler <-list(#tolerance = geneSet_tolerance,
                            #treatment = geneSet_treatment, 
                            #interaction = geneSet_interaction,
                            #none = noneModulesGO,
                            #"GO:0006281" = DNAR)
#noneModule_repair_euler_plot <- euler(noneModule_repair_euler)
##jpeg("sigModulesGO_none_repair_euler.jpg")
#plot(noneModule_repair_euler_plot, fills = c(plotColors[1:4], plotColors[6]))
##dev.off()
# "salmon4" "GO:0006298" 
salmon4ModulesGO <- c(resultsTable[resultsTable$color == "salmon4",]$gene)
names(salmon4ModulesGO) <- rep("salmon4", length(salmon4ModulesGO))
# "steelblue" "GO:0036297"
steelblueModulesGO <- c(resultsTable[resultsTable$color == "steelblue",]$gene)
names(steelblueModulesGO) <- rep("steelblue", length(steelblueModulesGO))
# "floralwhite" "GO:0006284"
# "floralwhite" "GO:0006289"
floralwhiteModulesGO <- c(resultsTable[resultsTable$color == "floralwhite",]$gene)
names(floralwhiteModulesGO) <- rep("floralwhite", length(floralwhiteModulesGO))
# euler diagram with interesting repair GO terms
modules_repair_euler <-list(#Tolerance = geneSet_tolerance,
                            #Treatment = geneSet_treatment, 
                            #Interaction = geneSet_interaction,
                            salmon4 = salmon4ModulesGO,
                            steelblue = steelblueModulesGO,
                            floralwhite = floralwhiteModulesGO,
                            #"GO:0006281" = DNAR,
                            "GO:0006298" = MMR,
                            "GO:0036297" = ICLR,
                            "GO:0006284" = BER, 
                            "GO:0006289" = NER)
modules_repair_euler_plot <- euler(modules_repair_euler)
#jpeg("sigModules_GO_salmon4_steelblue_repair_euler.jpg")
plot(modules_repair_euler_plot, fills = plotColors)
#dev.off()

# euler diagram with interesting repair GO terms for salmon4
modules_repair_salmon4_euler <-list(#tolerance = geneSet_tolerance,
                            treatment = geneSet_treatment, 
                            #interaction = geneSet_interaction,
                            salmon4 = salmon4ModulesGO,
                            #"GO:0006281" = DNAR,
                            "GO:0006298" = MMR)
modules_repair_salmon4_euler_plot <- euler(modules_repair_salmon4_euler)
#jpeg("sigModules_GO_salmon4_repair_euler.jpg")
plot(modules_repair_salmon4_euler_plot, fills = plotColors, quantities = list(type = "counts"))
#dev.off()

# euler diagram with interesting repair GO terms for salmon4
modules_repair_DEGs_salmon4_euler <-list(tolerance = geneSet_tolerance,
                                    treatment = geneSet_treatment, 
                                    interaction = geneSet_interaction,
                                    salmon4 = salmon4ModulesGO,
                                    "GO:0006281" = DNAR,
                                    "GO:0006298" = MMR)
modules_repair_DEGs_salmon4_euler_plot <- euler(modules_repair_DEGs_salmon4_euler)
#jpeg("sigModules_GO_DEGs_salmon4_repair_euler.jpg")
plot(modules_repair_DEGs_salmon4_euler_plot, fills = plotColors, quantities = list(type = "counts"))
#dev.off()

# significant radiation module diagram
radiation <- unlist(allGO["GO:0009314"])
cell <- unlist(allGO["GO:0071478"])
# "purple" "GO:0071478"
purpleModulesGO <- c(resultsTable[resultsTable$color == "steelblue",]$gene)
names(purpleModulesGO) <- rep("steelblue", length(purpleModulesGO))
# euler diagram with interesting repair GO terms
modules_radiation_euler <-list(purple = purpleModulesGO,
                            "GO:0071478" = cell,
                            "GO:0009314" = radiation)
modules_radiation_euler_plot <- euler(modules_radiation_euler)
#jpeg("sigModules_GO_purple_radiation_euler.jpg")
plot(modules_radiation_euler_plot, fills = plotColors[1:2])
#dev.off()

# significant stress modules diagram
stress <- unlist(allGO["GO:0006950"]) # interaction sig
oxidative <- unlist(allGO["GO:0006979"]) # interaction sig
cellOxidative <- unlist(allGO["GO:0034599"]) # treatment sig
cellReactiveOxy <- unlist(allGO["GO:0034614"])
# "darkolivegreen" "GO:0034614" 
darkolivegreenblueModulesGO <- c(resultsTable[resultsTable$color == "darkolivegreen",]$gene)
names(darkolivegreenblueModulesGO) <- rep("darkolivegreen", length(darkolivegreenblueModulesGO))
# "salmon4" "GO:0034599"
# euler diagram with interesting repair GO terms
modules_stress_euler <-list(#tolerance = geneSet_tolerance,
                            #treatment = geneSet_treatment, 
                            #interaction = geneSet_interaction,
                            salmon4 = salmon4ModulesGO,
                            darkolivegreen = darkolivegreenblueModulesGO,
                            "GO:0034614" = cellReactiveOxy,
                            "GO:0034599" = cellOxidative,
                            "GO:0006950" = stress,
                            "GO:0006979" = oxidative)
modules_stress_euler_plot <- euler(modules_stress_euler)
#jpeg("sigModules_GO_salmon4_darkolivegreen_stress_euler.jpg")
plot(modules_stress_euler_plot, fills = plotColors)
#dev.off()

# PUB
# salmon4 and DEGs
salmon4_euler <-list(Tolerance = geneSet_tolerance,
                     salmon4 = salmon4ModulesGO,
                     Treatment = geneSet_treatment, 
                     Interaction = geneSet_interaction)
                     #Positive = geneSet_positive)
salmon4_euler_plot <- euler(salmon4_euler)
#jpeg("sigModules_DEGs_salmon4_euler.jpg")
plot(salmon4_euler_plot, fills = plotColors, quantities = list(type = "counts"), labels = NULL, legend = list(labels = c("Tolerance", "salmon4", "Treatment", "Interaction")))
#dev.off()

# PUB
# salmon4 and DEGs and GO
salmon4_GO_euler <-list(Treatment = geneSet_treatment, 
                        salmon4 = salmon4ModulesGO,
                        "GO:0006298" = MMR,
                        "GO:0034599" = cellOxidative)
                        #Positive = geneSet_positive)
salmon4_GO_euler_plot <- euler(salmon4_GO_euler)
#jpeg("sigModules_GO_DEGs_GO_positive_salmon4_euler.jpg")
plot(salmon4_GO_euler_plot, fills = c(plotColors[1:4]), quantities = list(type = "counts"), labels = NULL, legend = list(labels = c("Treatment", "salmon4", "GO:0006298", "GO:0034599")))
#dev.off()

## PUB?
# potentially locally adapted and UVR responsive modules with DEGs
sienna3ModulesGO <- c(resultsTable[resultsTable$color == "sienna3",]$gene)
names(sienna3ModulesGO) <- rep("sienna3", length(sienna3ModulesGO))
skyblueModulesGO <- c(resultsTable[resultsTable$color == "skyblue",]$gene)
names(skyblueModulesGO) <- rep("skyblue", length(skyblueModulesGO))
salmon4_skyblue_sienna3_DEGs_euler <-list(skyblue = skyblueModulesGO,
                             salmon4 = salmon4ModulesGO,
                             sienna3 = sienna3ModulesGO,
                             Treatment = geneSet_treatment)
                             #Positive = geneSet_positive)
salmon4_skyblue_sienna3_DEGs_euler_plot <- euler(salmon4_skyblue_sienna3_DEGs_euler)
#jpeg("sigModules_treatment_salmon4_skyblue_sienna3_euler.jpg")
plot(salmon4_skyblue_sienna3_DEGs_euler_plot, fills = c(plotColors[1:3], plotColors[6]), quantities = list(type = "counts"))
#dev.off()

# potentially locally adapted and UVR responsive modules with DEGs
salmon4_skyblue_sienna3_euler <-list(Tolerance = geneSet_tolerance,
                             Treatment = geneSet_treatment, 
                             Interaction = geneSet_interaction,
                             salmon4 = salmon4ModulesGO,
                             skyblue = skyblueModulesGO,
                             sienna3 = sienna3ModulesGO)
#Positive = geneSet_positive)
salmon4_skyblue_sienna3_euler_plot <- euler(salmon4_skyblue_sienna3_euler)
#jpeg("sigModules_DEGs_salmon4_skyblue_sienna3_euler.jpg")
plot(salmon4_skyblue_sienna3_euler_plot, fills = plotColors, quantities = list(type = "counts"))
#dev.off()

# salmon4 and skyblue and DEGs
salmon4_skyblue_GO_euler <-list(Treatment = geneSet_treatment, 
                                salmon4 = salmon4ModulesGO,
                                skyblue = skyblueModulesGO,
                                "GO:0006298" = MMR,
                                "GO:0034599" = cellOxidative)
                                #Positive = geneSet_positive)
salmon4_skyblue_GO_euler_plot <- euler(salmon4_skyblue_GO_euler)
#jpeg("sigModules_GO_DEGs_salmon4_skyblue_euler.jpg")
plot(salmon4_skyblue_GO_euler_plot, fills = plotColors, quantities = list(type = "counts"))
#dev.off()

# PUB
# lightyellow and DEGs
PRR <- unlist(data.frame(unlist(allGO["GO:0006301"])), use.names = FALSE)
lightyellowModulesGO <- c(resultsTable[resultsTable$color == "lightyellow",]$gene)
names(lightyellowModulesGO) <- rep("lightyellow", length(lightyellowModulesGO))
lightyellow_GO_euler <-list(Tolerance = geneSet_tolerance,
                                    Treatment = geneSet_treatment, 
                                    Interaction = geneSet_interaction,
                                    #salmon4 = salmon4ModulesGO,
                                    lightyellow = lightyellowModulesGO)
                                    #"GO:0006298" = MMR,
                                    #"GO:0034599" = cellOxidative)
                                    #"GO:0006301" = PRR)
                                    #Positive = geneSet_positive)
lightyellow_GO_euler_plot <- euler(lightyellow_GO_euler)
#jpeg("sigModules_GO_DEGs_lightyellow_euler.jpg")
plot(lightyellow_GO_euler_plot, fills = plotColors, quantities = list(type = "counts"), labels = NULL, legend = list(labels = c("Tolerance", "Treatment", "Interaction", "lightyellow")))
#dev.off()
