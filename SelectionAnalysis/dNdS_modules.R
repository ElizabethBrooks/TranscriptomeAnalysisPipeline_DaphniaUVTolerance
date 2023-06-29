#!/usr/bin/env Rscript

# load libraries
library(ggplot2)
library(ghibli)

# turn off scientific notation
options(scipen = 999)

# set the working directory
workingDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/SelectionTests"
setwd(workingDir)

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

#GO enrichment
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


# Plotting Palettes
# retrieve the vector of colors associated with PonyoMedium
ghibli_colors <- ghibli_palette("PonyoMedium", type = "discrete")
# vector with a subset of colors associated with PonyoMedium
ghibli_subset <- c(ghibli_colors[3], ghibli_colors[6], ghibli_colors[4])

# retrieve dN dS values
inputTable <- read.csv(file="Pulex_Olympics_kaksResults.csv", row.names="geneID")


# merging data frames
# https://sparkbyexamples.com/r-programming/r-join-data-frames-with-examples/#full-outer-join

# remove genes with dN/dS > 7
subsetTable <- inputTable[inputTable$dNdS < 7,]
#subsetTable <- inputTable

# add geneID column
geneID <- row.names(subsetTable)
subsetTable <- cbind(geneID,subsetTable)
colnames(resultsTable)[1] ="geneID"

# full outer join data frames
resultsTable <- merge(x = resultsTable, y = subsetTable, 
                          by = "geneID", all=TRUE)

# remove rows with NAs
resultsTable <- na.omit(resultsTable)

# add column for identifying mode of selection
resultsTable$Selection <- "NA"
resultsTable$Selection[resultsTable$dNdS > 1] <- "Positive"
resultsTable$Selection[resultsTable$dNdS < 1] <- "Negative"


# box plot colored by selection
jpeg("dNdS_modules_boxPlot.jpg")
ggplot(resultsTable, aes(x=color, y=dNdS)) +
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# box plot colored by selection
jpeg("dNdS_modules_selection_boxPlot.jpg")
ggplot(resultsTable, aes(x=color, y=dNdS, fill=Selection)) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_discrete(type = ghibli_subset, breaks = c("Positive", "Negative")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# violin plots
# http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization

# violin plot color by selection
jpeg("dNdS_modules_violinPlot.jpg")
ggplot(resultsTable, aes(x=color, y=dNdS)) +
  theme_minimal() +
  geom_violin(trim=FALSE) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# violin plot color by selection
jpeg("dNdS_modules_selection_violinPlot.jpg")
ggplot(resultsTable, aes(x=color, y=dNdS, fill=Selection)) +
  theme_minimal() +
  geom_violin(trim=FALSE) +
  scale_fill_discrete(type = ghibli_subset, breaks = c("Positive", "Negative")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

