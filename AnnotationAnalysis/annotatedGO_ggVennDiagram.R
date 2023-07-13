#!/usr/bin/env Rscript

#Install edgeR and statmod, this should only need to be done once
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("edgeR")
#install.packages("ggrepel")

# turn off scientific notation
options(scipen = 999)

# load librarys
library(ggplot2)
library(ggrepel)
library(ggVennDiagram)
library(rcartocolor)

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

# import GO annotations
goTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/email/geneToGO_map.txt", sep = "\t")


# import DEGs
interactionTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Tolerance/glmQLF_2WayANOVA_interaction_topTags_LFC1.2.csv", row.names="gene")
treatmentTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Tolerance/glmQLF_2WayANOVA_treatment_topTags_LFC1.2.csv", row.names="gene")
toleranceTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Tolerance/glmQLF_2WayANOVA_tolerance_topTags_LFC1.2.csv", row.names="gene")

# keep only sig
interactionSig <- interactionTable[interactionTable$FDR < 0.05,]
treatmentSig <- treatmentTable[treatmentTable$FDR < 0.05,]
toleranceSig <- toleranceTable[toleranceTable$FDR < 0.05,]

# import modules




# GO annotations

# setup gene sets
geneSet_KAP4 <- rownames()
geneSet_DEGs <- rownames()
geneSet_Modules <- rownames()

# create combined list of gene sets
geneSet_list <- list(KAP4 = geneSet_KAP4, 
                    DEGs = geneSet_DEGs,
                    Modules = geneSet_Modules)

# create venn diagram
jpeg("annotatedGO_venn.jpg")
ggVennDiagram(glm_list_venn, label_alpha=0.25, category.names = c("KAP4","DEGs","Modules")) +
  scale_colour_discrete(type = plotColorSubset)
dev.off()

