#!/usr/bin/env Rscript

# install packages
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("edgeR")
#install.packages("ggrepel")
#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/ggpubr")
#install.packages("ggpubr")
#install.packages("eulerr")

# turn off scientific notation
options(scipen = 999)

# load librarys
library(topGO)
library(ggplot2)
library(ggrepel)
library(ggVennDiagram)
library(rcartocolor)
library(ggpubr)
library(eulerr)

# Plotting Palettes
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
# https://github.com/Nowosad/rcartocolor
plotColors <- carto_pal(12, "Safe")
plotColorSubset <- c(plotColors[3], plotColors[2], plotColors[1])

# turn off scientific notation
options(scipen = 999)

# set the working directory
workingDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/AnnotationAnalysis/Genotypes"
setwd(workingDir)

# import module analysis results
anovaTable <- read.csv("/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Genotypes/DEGsANOVA_OLYM_30/aov_summary_pValues.csv")
positiveTable <- read.csv("/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/selectionTests/Genotypes/fisherTest_positiveSelection_modules.csv")

# update anova module and column names
anovaTable$module <- gsub("ME", "", anovaTable$module)
colnames(anovaTable) <- c("color","treatment","tolerance","toleranceGenotype","treatmentTolerance","treatmentToleranceGenotype")

# import GO annotations
goTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/email/geneToGO_map.txt", sep = "\t", header = FALSE)

# update column names
colnames(goTable) <- c("geneID", "GO")

# add gene tag
goTable$geneID <- paste("gene", goTable$geneID, sep="-")

# retrieve dN dS values
dNdSTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/SelectionTests/Pulex_Olympics_kaksResults.csv", row.names="geneID")

# filter to keep dNdS > 1 & dNdS < 99
dNdSSubset <- dNdSTable[dNdSTable$dNdS > 1 & dNdSTable$dNdS < 99,]

# import DEGs
interactionTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Genotypes/glmQLF_2WayANOVA_interaction_topTags_LFC1.2.csv", row.names="gene")
treatmentTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Genotypes/glmQLF_2WayANOVA_treatment_topTags_LFC1.2.csv", row.names="gene")
toleranceTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Genotypes/glmQLF_2WayANOVA_tolerance_topTags_LFC1.2.csv", row.names="gene")

# keep only sig
interactionSig <- interactionTable[interactionTable$FDR < 0.05,]
treatmentSig <- treatmentTable[treatmentTable$FDR < 0.05,]
toleranceSig <- toleranceTable[toleranceTable$FDR < 0.05,]

# combine all tables
effectTable <- rbind(interactionSig, treatmentSig, toleranceSig)

# add geneID column
effectTable$geneID <- row.names(effectTable)

# import modules
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

# full outer join data frames
moduleCombinedTable <- merge(x = modSizes, y = anovaTable, 
                      by = "color", all=TRUE)
moduleCombinedTable <- merge(x = moduleCombinedTable, y = positiveTable, 
                             by = "color", all=TRUE)

# export module info
outDir="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Genotypes"
exportFile <- paste(tag, "moduleInfo_size_anova_positive.csv", sep="-")
exportFile <- paste(outDir, exportFile, sep="/")
write.csv(file=exportFile, moduleCombinedTable, row.names=FALSE)


# DEGs

# setup gene sets
geneSet_treatment <- rownames(treatmentSig)
geneSet_tolerance <- rownames(toleranceSig)
geneSet_interaction <- rownames(interactionSig)
glm_list_venn <- list(treatment = geneSet_treatment, 
                      tolerance = geneSet_tolerance,
                      interaction = geneSet_interaction)

# euler diagram
euler_plot <- euler(glm_list_venn)
plot(euler_plot, fills = plotColorSubset)

# create venn diagram
jpeg("DEGs_venn.jpg")
ggVennDiagram(glm_list_venn, label_alpha=0.25, category.names = c("Treatment","Tolerance","Interaction")) +
  scale_colour_discrete(type = plotColorSubset)
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
oxidativeStress <- unlist(allGO["GO:0034599"]) # treatment sig
hemeProcess <- unlist(allGO["GO:0006783"]) # interaction and treatment sig
alphaProcess <- unlist(allGO["GO:1901607"]) # interaction and treatment sig
# The radiation response terms included response to radiation	(GO:0009314), response to ionizing radiation (GO:0010212), response to gamma radiation (GO:0010332), cellular response to radiation (GO:0071478), cellular response to gamma radiation (GO:0071480), cellular response to ionizing radiation (GO:0071479), regulation of response to gamma radiation (GO:2001228), regulation of cellular response to gamma radiation (GO:1905843), positive regulation of response to gamma radiation (GO:2001230), negative regulation of response to gamma radiation (GO:2001229), 
# response to X-ray (GO:0010165), positive regulation of cellular response to gamma radiation	(GO:1905845), negative regulation of cellular response to gamma radiation (GO:1905844), cellular response to X-ray (GO:0071481), regulation of cellular response to X-ray (GO:2000683), negative regulation of cellular response to X-ray (GO:2000684), positive regulation of cellular response to X-ray (GO:2000685), 
# phototransduction UV (GO:0007604), response to UV (GO:0009411), response to UV-A (GO:0070141), detection of UV (GO:0009589), cellular response to UV (GO:0034644), cellular response to UV-A (GO:0071492), response to UV-B	4 (GO:0010224), cellular response to UV-B (GO:0071493), response to UV-C (GO:0010225), 
# regulation of mRNA stability involved in cellular response to UV (GO:1902629), cellular response to UV-C (GO:0071494), regulation of translation involved in cellular response to UV (GO:1904803).
radiation <- unlist(allGO["GO:0009314"]) # not sig
#ionizing <- unlist(allGO["GO:0010212"])
#gamma <- unlist(allGO["GO:0010332"])
cell <- unlist(allGO["GO:0071478"]) # not sig
#cellGamma <- unlist(allGO["GO:0071480"])
#cellIonizing <- unlist(allGO["GO:0071479"])
#regGamma <- unlist(allGO["GO:2001228"])
#regCellGamma <- unlist(allGO["GO:1905843"])
#posRegGamma <- unlist(allGO["GO:2001230"])
#negRegGamma <- unlist(allGO["GO:2001229"])
#xRay <- unlist(allGO["GO:0010165"])
#posRegCellGamma <- unlist(allGO["GO:1905845"])
#negRegCellGamma <- unlist(allGO["GO:1905844"])
#cellXRay <- unlist(allGO["GO:0071481"])
#regCellXRay <- unlist(allGO["GO:2000683"])
#negRegCellXRay <- unlist(allGO["GO:2000684"])
#posRegCellXRay <- unlist(allGO["GO:2000685"])
#photoUV <- unlist(allGO["GO:0007604"])
#UV <- unlist(allGO["GO:0009411"])
#UVA <- unlist(allGO["GO:0070141"])
#detectUV <- unlist(allGO["GO:0009589"])
#cellUV <- unlist(allGO["GO:0034644"])
#cellUVA <- unlist(allGO["GO:0071492"])
#UVB4 <- unlist(allGO["GO:0010224"])
#cellUVB <- unlist(allGO["GO:0071493"])
#UVC <- unlist(allGO["GO:0010225"])
#mRNACellUV <- unlist(allGO["GO:1902629"])
#cellUVC <- unlist(allGO["GO:0071494"])
#transCellUV <- unlist(allGO["GO:1904803"])
# The stress terms included response to stress (GO:0006950), response to oxidative stress (GO:0006979), cellular response to oxidative stress (GO:0034599), 
# positive regulation of translation in response to oxidative stress (GO:0032939), cellular response to reactive oxygen (GO:0034614), negative regulation of cellular response to oxidative stress (GO:1900408), positive regulation of cellular response to oxidative stress (GO:1900409), 
# negative regulation of translation in response to oxidative stress (GO:0032938), regulation of translation in response to oxidative stress (GO:0043556), Mcs4 RR-MAPKKK complex (GO:1990315), cellular response to hydroperoxide (GO:0071447), 
# oxidative stress-induced premature senescence (GO:0090403), symbiont defense to host-produced reactive oxygen species (GO:0052164), regulation of cellular response to oxidative stress (GO:1900407).
stress <- unlist(allGO["GO:0006950"]) # interaction sig
oxidative <- unlist(allGO["GO:0006979"]) # interaction sig
cellOxidative <- unlist(allGO["GO:0034599"]) # treatment sig
#posRegTransOxidative <- unlist(allGO["GO:0032939"])
cellReactiveOxy <- unlist(allGO["GO:0034614"]) # not sig
#negRegCellOxidative <- unlist(allGO["GO:1900408"])
#posRegCellOxidative <- unlist(allGO["GO:1900409"])
#negRegTransOxidative <- unlist(allGO["GO:0032938"])
#regTransOxidative <- unlist(allGO["GO:0043556"])
#MAPKKK <- unlist(allGO["GO:1990315"])
#hydroperoxide <- unlist(allGO["GO:0071447"])
#senescence <- unlist(allGO["GO:0090403"])
#symbiont <- unlist(allGO["GO:0052164"])
#regCellOxidative <- unlist(allGO["GO:1900407"])

# euler diagram with interesting radiation GO terms
glm_list_venn_radGO <-list(treatment = geneSet_treatment, 
                        tolerance = geneSet_tolerance,
                        interaction = geneSet_interaction,
                        radiation = radiation,
                        cell = cell)
euler_plot_radGO <- euler(glm_list_venn_radGO)
plot(euler_plot_radGO)

# euler diagram with interesting repair GO terms
glm_list_venn_repairGO <-list(treatment = geneSet_treatment, 
                        tolerance = geneSet_tolerance,
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
plot(euler_plot_repairGO)

# euler diagram with significant stress GO terms
glm_list_venn_stressGO <-list(treatment = geneSet_treatment, 
                        tolerance = geneSet_tolerance,
                        interaction = geneSet_interaction,
                        stress = stress,
                        oxidative = oxidative,
                        cell = cellOxidative)
                        #cellReactiveOxy = cellReactiveOxy)
euler_plot_stressGO <- euler(glm_list_venn_stressGO)
plot(euler_plot_stressGO)

# euler diagram with significant GO terms
glm_list_venn_sigGO <-list(treatment = geneSet_treatment, 
                        tolerance = geneSet_tolerance,
                        interaction = geneSet_interaction,
                        oxidative = oxidative,
                        stress = stress)
                        #heme = hemeProcess,
                        #alpha = alphaProcess)
euler_plot_sigGO <- euler(glm_list_venn_sigGO)
plot(euler_plot_sigGO)

# euler diagram with interesting all GO terms
#glm_list_venn_allGO <-list(treatment = geneSet_treatment, 
#                              tolerance = geneSet_tolerance,
#                              interaction = geneSet_interaction,
#                              DNAR = DNAR,
#                              MMR = MMR,
#                              BER = BER,
#                              HR = HR,
#                              NER = NER,
#                              ICLR = ICLR,
#                              DSBR = DSBR,
#                              stress = stress,
#                              oxidative = oxidative,
#                              cellOxidative = cellOxidative,
#                              cellReactiveOxy = cellReactiveOxy,
#                              radiation = radiation,
#                              cell = cell)
#SSBR = SSBR)
#euler_plot_allGO <- euler(glm_list_venn_allGO)
#plot(euler_plot_allGO)
  
# setup gene sets
geneSet_KAP4 <- goTable$geneID
geneSet_DEGs <- effectTable$geneID
geneSet_Modules <- resultsTable$gene

# create combined list of gene sets
annotationSet_list <- list(KAP4 = geneSet_KAP4, 
                    DEGs = geneSet_DEGs,
                    Modules = geneSet_Modules)

# create venn diagram
jpeg("annotatedGO_DEGs_modules_venn.jpg")
ggVennDiagram(annotationSet_list, label_alpha=0.25, category.names = c("KAP4","DEGs","Modules")) +
  scale_colour_discrete(type = plotColorSubset)
dev.off()


# positively selected

# setup gene sets
geneSet_positive <- row.names(dNdSSubset)

# create combined list of gene sets
positiveSet_list <- list(Positive = geneSet_positive, 
                     DEGs = geneSet_DEGs,
                     Modules = geneSet_Modules)

# create venn diagram
jpeg("positive_DEGs_modules_venn.jpg")
ggVennDiagram(positiveSet_list, label_alpha=0.25, category.names = c("Positive","DEGs","Modules")) +
  scale_colour_discrete(type = plotColorSubset)
dev.off()


# combined

# create combined list of gene sets
combinedSet_list <- list(KAP4 = geneSet_KAP4,
                        Positive = geneSet_positive, 
                        DEGs = geneSet_DEGs,
                        Modules = geneSet_Modules)

# create venn diagram
jpeg("annotatedGO_positive_DEGs_modules_venn.jpg")
ggVennDiagram(combinedSet_list, label_alpha=0.25, category.names = c("GO","Positive","DE","Modules")) +
  scale_colour_discrete(type = c(plotColorSubset, plotColors[10]))
dev.off()


# combine venn diagrams
# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/

# store each venn
effectVenn <- ggVennDiagram(glm_list_venn, label_alpha=0.25, category.names = c("Treatment","Tolerance","Interaction")) +
  scale_colour_discrete(type = plotColorSubset)
annotationVenn <- ggVennDiagram(annotationSet_list, label_alpha=0.25, category.names = c("KAP4","DEGs","Modules")) +
  scale_colour_discrete(type = plotColorSubset)
positiveVenn <- ggVennDiagram(positiveSet_list, label_alpha=0.25, category.names = c("Positive","DEGs","Modules")) +
  scale_colour_discrete(type = plotColorSubset)
combinedVenn <- ggVennDiagram(combinedSet_list, label_alpha=0.25, category.names = c("GO","Positive","DE","Modules")) +
  scale_colour_discrete(type = c(plotColorSubset, plotColors[10]))

# one figure in row 1 and two figures in row 2
jpeg("DEGs_combined_venn.jpg")
ggarrange(effectVenn, combinedVenn, ncol = 2, labels = c("A", "B"))
dev.off()

# one figure in row 1 and two figures in row 2
#ggarrange(effectVenn,
#          ggarrange(annotationVenn, positiveVenn, ncol = 2, labels = c("B", "C")),
#          nrow = 2, 
#          labels = "A"
#) 

