#!/usr/bin/env Rscript

# load libraries
library(ggplot2)
library(ghibli)

# turn off scientific notation
options(scipen = 999)

# set the working directory
workingDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/SelectionTests"
setwd(workingDir)

# Plotting Palettes
# retrieve the vector of colors associated with PonyoMedium
ghibli_colors <- ghibli_palette("PonyoMedium", type = "discrete")
# vector with a subset of colors associated with PonyoMedium
ghibli_subset <- c(ghibli_colors[3], ghibli_colors[6], ghibli_colors[4])

# retrieve dN dS values
inputTable <- read.csv(file="Pulex_Olympics_kaksResults.csv", row.names="geneID")


# merging data frames
# https://sparkbyexamples.com/r-programming/r-join-data-frames-with-examples/#full-outer-join

# remove genes with dN/dS > 3
subsetTable <- inputTable[inputTable$dNdS < 3,]

# retrieve DEGs
interactionTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Tolerance/glmQLF_2WayANOVA_interaction_topTags_LFC1.2.csv")
treatmentTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Tolerance/glmQLF_2WayANOVA_treatment_topTags_LFC1.2.csv")
toleranceTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Tolerance/glmQLF_2WayANOVA_tolerance_topTags_LFC1.2.csv")

# keep only sig
interactionTable <- subset(interactionTable, interactionTable$FDR < 0.05)
treatmentTable <- subset(treatmentTable, treatmentTable$FDR < 0.05)
toleranceTable <- subset(toleranceTable, toleranceTable$FDR < 0.05)

# add effect tags
interactionTable$Effect <- "Interaction"
treatmentTable$Effect <- "Treatment"
toleranceTable$Effect <- "Tolerance"

# add geneID column
geneID <- row.names(subsetTable)
subsetTable <- cbind(geneID,subsetTable)
geneID <- row.names(interactionTable)
interactionTable <- cbind(geneID,interactionTable)
geneID <- row.names(treatmentTable)
treatmentTable <- cbind(geneID,treatmentTable)
geneID <- row.names(toleranceTable)
toleranceTable <- cbind(geneID,toleranceTable)

# full outer join data frames
interactionTable <- merge(x = interactionTable, y = subsetTable, 
                          by = "geneID", all=TRUE)
treatmentTable <- merge(x = treatmentTable, y = subsetTable, 
                        by = "geneID", all=TRUE)
toleranceTable <- merge(x = toleranceTable, y = subsetTable, 
                          by = "geneID", all=TRUE)

# remove rows with NAs
interactionTable <- na.omit(interactionTable)
treatmentTable <- na.omit(treatmentTable)
toleranceTable <- na.omit(toleranceTable)

# combine all tables
all_plotTable <- rbind(interactionTable, treatmentTable, toleranceTable)

# add column for identifying mode of selection
all_plotTable$Selection <- "NA"
all_plotTable$Selection[all_plotTable$dNdS > 1] <- "Positive"
all_plotTable$Selection[all_plotTable$dNdS < 1] <- "Negative"

# create label set
labelSetAll <- all_plotTable[all_plotTable$dNdS > 1,]

# box plot colored by selection
jpeg("dNdS_DEGs_boxPlot.jpg")
ggplot(all_plotTable, aes(x=Effect, y=dNdS)) +
  geom_boxplot() +
  theme_minimal()
dev.off()

# box plot colored by selection
jpeg("dNdS_DEGs_selection_boxPlot.jpg")
ggplot(all_plotTable, aes(x=Effect, y=dNdS, fill=Selection)) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_discrete(type = ghibli_subset, breaks = c("Positive", "Negative"))
dev.off()

# box plot colored by selection with labels
jpeg("dNdS_DEGs_selection_boxPlot_labeled.jpg")
ggplot(all_plotTable, aes(x=Effect, y=dNdS, fill=Selection)) +
  geom_boxplot() +
  theme_minimal() +
  ggrepel::geom_text_repel(data = labelSetAll, aes(label = labelSetAll$geneID), max.overlaps=100) +
  scale_fill_discrete(type = ghibli_subset, breaks = c("Positive", "Negative"))
dev.off()


# violin plots
# http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization

# violin plot color by selection
jpeg("dNdS_DEGs_violinPlot.jpg")
ggplot(all_plotTable, aes(x=Effect, y=dNdS)) +
  theme_minimal() +
  geom_violin(trim=FALSE) 
dev.off()

# violin plot color by selection
jpeg("dNdS_DEGs_selection_violinPlot.jpg")
ggplot(all_plotTable, aes(x=Effect, y=dNdS, fill=Selection)) +
  theme_minimal() +
  geom_violin(trim=FALSE) +
  scale_fill_discrete(type = ghibli_subset, breaks = c("Positive", "Negative"))
dev.off()

# violin plot color by selection with labels
jpeg("dNdS_modules_selection_violinPlot_labeled.jpg")
ggplot(all_plotTable, aes(x=Effect, y=dNdS, fill=Selection)) +
  theme_minimal() +
  ggrepel::geom_text_repel(data = labelSetAll, aes(label = labelSetAll$geneID), max.overlaps=100) +
  geom_violin(trim=FALSE) +
  scale_fill_discrete(type = ghibli_subset, breaks = c("Positive", "Negative")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
