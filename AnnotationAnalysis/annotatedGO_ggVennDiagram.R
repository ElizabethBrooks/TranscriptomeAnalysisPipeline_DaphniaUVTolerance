#!/usr/bin/env Rscript

# install packages
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("edgeR")
#install.packages("ggrepel")
#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/ggpubr")
#install.packages("ggpubr")

# turn off scientific notation
options(scipen = 999)

# load librarys
library(ggplot2)
library(ggrepel)
library(ggVennDiagram)
library(rcartocolor)
library(ggpubr)

# Plotting Palettes
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
# https://github.com/Nowosad/rcartocolor
plotColors <- carto_pal(12, "Safe")
plotColorSubset <- c(plotColors[4], plotColors[5], plotColors[6])

# turn off scientific notation
options(scipen = 999)

# set the working directory
workingDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/AnnotationAnalysis"
setwd(workingDir)

# import module analysis results
anovaTable <- read.csv("/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Tolerance/ANOVA_OLYM_30/aov_summary_pValues.csv")
positiveTable <- read.csv("/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/selectionTests/fisherTest_positiveSelection_modules.csv")

# update anova module and column names
anovaTable$module <- gsub("ME", "", anovaTable$module)
colnames(anovaTable) <- c("color","treatment","tolerance","interaction")

# import GO annotations
goTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/email/geneToGO_map.txt", sep = "\t", header = FALSE)

# update column names
colnames(goTable) <- c("geneID", "GO")

# add gene tag
goTable$geneID <- paste("gene", goTable$geneID, sep="-")

# retrieve dN dS values
dNdSTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/SelectionTests/Pulex_Olympics_kaksResults.csv", row.names="geneID")

# import DEGs
interactionTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Tolerance/glmQLF_2WayANOVA_interaction_topTags_LFC1.2.csv", row.names="gene")
treatmentTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Tolerance/glmQLF_2WayANOVA_treatment_topTags_LFC1.2.csv", row.names="gene")
toleranceTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Tolerance/glmQLF_2WayANOVA_tolerance_topTags_LFC1.2.csv", row.names="gene")

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
inDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Tolerance"

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
  color = character()
)

# match gene IDs with module colors
for(i in 1:numMods){
  gene <- names(datExpr)[moduleColors==colorTable[i,1]]
  color <- rep(colorTable[i,1], length(gene))
  moduleData <- cbind(gene, color)
  resultsTable <- rbind(resultsTable, moduleData)
}

# setup data frame of module sizes
modSizes <- data.frame(
  color = character(),
  size = numeric()
)

# retrieve module sizes
for(i in 1:numMods){
  color <- head(resultsTable[resultsTable$number == i,2],1)
  size <- nrow(resultsTable[resultsTable$number == i,])
  moduleData <- cbind(color, size)
  modSizes <-rbind(modSizes, moduleData)
}

# full outer join data frames
moduleCombinedTable <- merge(x = modSizes, y = anovaTable, 
                      by = "color", all=TRUE)
moduleCombinedTable <- merge(x = moduleCombinedTable, y = positiveTable, 
                             by = "color", all=TRUE)

# export module info
outDir="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Tolerance"
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
# create venn diagram
jpeg("DEGs_venn.jpg")
ggVennDiagram(glm_list_venn, label_alpha=0.25, category.names = c("Treatment","Tolerance","Interaction")) +
  scale_colour_discrete(type = plotColorSubset)
dev.off()


# GO annotations

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
geneSet_positive <- row.names(dNdSTable)

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
  scale_colour_discrete(type = c(plotColorSubset,"#661100"))
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

# one figure in row 1 and two figures in row 2
ggarrange(effectVenn,                                                 # First row with scatter plot
          ggarrange(annotationVenn, positiveVenn, ncol = 2, labels = c("B", "C")), # Second row with box and dot plots
          nrow = 2, 
          labels = "A"                                        # Labels of the scatter plot
) 
