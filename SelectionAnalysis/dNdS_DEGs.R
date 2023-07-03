#!/usr/bin/env Rscript

# load libraries
library(ggplot2)
library(rcartocolor)
library(tidyr)

# turn off scientific notation
options(scipen = 999)

# set the working directory
workingDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/SelectionTests"
setwd(workingDir)

# Plotting Palettes
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
# https://github.com/Nowosad/rcartocolor
plotColors <- carto_pal(12, "Safe")
plotColorSubset <- c(plotColors[4], plotColors[11], plotColors[5], plotColors[9])

# retrieve dN dS values
dNdSTable <- read.csv(file="Pulex_Olympics_kaksResults.csv", row.names="geneID")

# remove NAs
dNdSSubset <- na.omit(dNdSTable)

# view outliers
dNdSSubset[dNdSSubset$dS > 10,]

# subset dN dS values to remove outliers
dNdSSubset <- dNdSSubset[dNdSSubset$dS < 10,]

# retrieve DEGs
interactionTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Tolerance/glmQLF_2WayANOVA_interaction_topTags_LFC1.2.csv")
treatmentTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Tolerance/glmQLF_2WayANOVA_treatment_topTags_LFC1.2.csv")
toleranceTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Tolerance/glmQLF_2WayANOVA_tolerance_topTags_LFC1.2.csv")

# keep only sig
interactionSig <- subset(interactionTable, interactionTable$FDR < 0.05)
treatmentSig <- subset(treatmentTable, treatmentTable$FDR < 0.05)
toleranceSig <- subset(toleranceTable, toleranceTable$FDR < 0.05)

# add effect tags
interactionSig$Effect <- "Interaction"
treatmentSig$Effect <- "Treatment"
toleranceSig$Effect <- "Tolerance"

# add geneID column
geneID <- row.names(dNdSSubset)
dNdSSubset <- cbind(geneID,dNdSSubset)
geneID <- row.names(interactionSig)
interactionSig <- cbind(geneID,interactionSig)
geneID <- row.names(treatmentSig)
treatmentSig <- cbind(geneID,treatmentSig)
geneID <- row.names(toleranceSig)
toleranceSig <- cbind(geneID,toleranceSig)

# keep necessary columns
dNdSSubset <- dNdSSubset[,c("geneID","dN","dS","dNdS")]
interactionSubset <- interactionSig[,c("geneID","Effect")]
treatmentSubset <- treatmentSig[,c("geneID","Effect")]
toleranceSubset <- toleranceSig[,c("geneID","Effect")]

# merging data frames
# https://sparkbyexamples.com/r-programming/r-join-data-frames-with-examples/#full-outer-join

# combine all tables
plotTableSubset <- rbind(interactionSubset, treatmentSubset, toleranceSubset)

# full outer join
plotTable <- merge(x = dNdSSubset, y = plotTableSubset, 
                   by = "geneID", all=TRUE)

# set effect tag for NAs
plotTable$Effect <- plotTable$Effect %>% replace_na('None')

# add selection type
plotTable$Selection[plotTable$dNdS == 99] <- "Error"
plotTable$Selection[plotTable$dNdS > 1 & plotTable$dNdS < 99] <- "Positive"
plotTable$Selection[plotTable$dNdS < 1] <- "Negative"

# remove NAs
plotTable <- na.omit(plotTable)

# set row order for plotting
plotTable$Order <- ifelse(plotTable$Effect=="None", 1, 2)
plotTable <- plotTable[order(plotTable$Order),]


# plotting

# box plot colored by selection
jpeg("dNdS_DEGs_boxPlot.jpg")
ggplot(plotTable, aes(x=Effect, y=dNdS)) +
  geom_boxplot() +
  theme_minimal()
dev.off()

# box plot colored by selection
jpeg("dNdS_DEGs_selection_boxPlot.jpg")
ggplot(plotTable, aes(x=Effect, y=dNdS, fill=Selection)) +
  geom_boxplot() +
  theme_minimal() +
  scale_colour_discrete(type = plotColorSubset)
dev.off()

# box plot colored by selection with labels
#jpeg("dNdS_DEGs_selection_boxPlot_labeled.jpg")
#ggplot(plotTable, aes(x=Effect, y=dNdS, fill=Selection)) +
#  geom_boxplot() +
#  theme_minimal() +
#  ggrepel::geom_text_repel(data = labelSetAll, aes(label = labelSetAll$geneID), max.overlaps=100) +
#  scale_colour_discrete(type = plotColorSubset)
#dev.off()

# scatter plots
# http://www.sthda.com/english/wiki/ggplot2-scatter-plots-quick-start-guide-r-software-and-data-visualization
# https://stackoverflow.com/questions/15015356/how-to-do-selective-labeling-with-ggplot-geom-point
# http://www.sthda.com/english/wiki/ggplot2-point-shapes
# https://www.statology.org/ggplot-multiple-data-frames/

# scatter plot colored by selection
# shaped by DEG effect set
jpeg("dNdS_DEGs_selection_scatterPlot.jpg")
ggplot(plotTable, aes(x=dS, y=dN, shape=Selection, color=Effect)) +
  theme_minimal() +
  geom_point() +
  scale_shape_manual(values=c(2, 1, 3)) +
  scale_colour_discrete(type = plotColorSubset)
dev.off()
