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

# remove genes with dN/dS = 99
# https://ocw.mit.edu/courses/6-877j-computational-evolutionary-biology-fall-2005/9a6d5e515fb1e7608eb3919855b01880_pamlfaqs.pdf
#subsetTable <- inputTable[inputTable$dNdS < 99,]
subsetTable <- inputTable

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
resultsTable$Selection[resultsTable$dNdS == 99] <- "Error"
resultsTable$Selection[resultsTable$dNdS > 1 & resultsTable$dNdS < 99] <- "Positive"
resultsTable$Selection[resultsTable$dNdS < 1] <- "Negative"

# create label set
labelSetNeg <- resultsTable[resultsTable$Selection == "Negative",]
labelSetPos <- resultsTable[resultsTable$Selection == "Positive",]
labelSetError <- resultsTable[resultsTable$Selection == "Error",]


# scatter plots
# http://www.sthda.com/english/wiki/ggplot2-scatter-plots-quick-start-guide-r-software-and-data-visualization
# https://stackoverflow.com/questions/15015356/how-to-do-selective-labeling-with-ggplot-geom-point
# http://www.sthda.com/english/wiki/ggplot2-point-shapes

# scatter plot colored by selection
# shaped by DEG effect set
jpeg("dNdS_modules_selection_scatterPlot_all.jpg")
ggplot(resultsTable, aes(x=dS, y=dN, shape=Selection, color=color)) +
  theme_minimal() +
  geom_point() +
  scale_shape_manual(values=c(2, 1, 3))
dev.off()

# scatter plot colored by selection
# shaped by DEG effect set
# with positvely selected genes labeled by geneID
jpeg("dNdS_modules_selection_scatterPlot_labeled_all.jpg")
ggplot(resultsTable, aes(x=dS, y=dN, shape=Selection, color=color)) +
  theme_minimal() +
  geom_point() +
  ggrepel::geom_text_repel(data = labelSetPos, aes(label = labelSetPos$geneID), max.overlaps=100, color="Black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_shape_manual(values=c(2, 1, 3))
dev.off()

# scatter plot colored by selection
# shaped by DEG effect set
# with error genes labeled by geneID
jpeg("dNdS_modules_selection_scatterPlot_labeled_error.jpg")
ggplot(resultsTable, aes(x=dS, y=dN, shape=Selection, color=color)) +
  theme_minimal() +
  geom_point() +
  ggrepel::geom_text_repel(data = labelSetNeg, aes(label = labelSetNeg$geneID), max.overlaps=100, color="Black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_shape_manual(values=c(2, 1, 3))
dev.off()

# scatter plot colored by selection
# shaped by DEG effect set
# with error genes labeled by geneID
jpeg("dNdS_modules_selection_scatterPlot_labeled_error.jpg")
ggplot(resultsTable, aes(x=dS, y=dN, shape=Selection, color=color)) +
  theme_minimal() +
  geom_point() +
  ggrepel::geom_text_repel(data = labelSetError, aes(label = labelSetError$geneID), max.overlaps=100, color="Black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_shape_manual(values=c(2, 1, 3))
dev.off()



