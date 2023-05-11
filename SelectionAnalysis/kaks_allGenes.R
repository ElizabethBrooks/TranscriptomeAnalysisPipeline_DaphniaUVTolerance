#!/usr/bin/env Rscript

# load libraries
library(ggplot2)
library(ghibli)

# turn off scientific notation
options(scipen = 999)

# set the working directory
setwd("/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Bioinformatics/variantsCalled_samtoolsBcftools/selectionTests")

# retrieve dN dS values
inputTable <- read.csv(file="Pulex_Olympics_kaksResults.csv", row.names="geneID")

# remove genes with dN/dS > 10
subsetTable <- inputTable[inputTable$dNdS < 10,]

# scatter dS vs dN
ggplot(subsetTable, aes(x=dS, y=dN)) + 
  geom_point()

# scatter dN/dS vs dS
ggplot(subsetTable, aes(x=dS, y=dNdS)) + 
  geom_point()

# genes under positive selection
posTable <- subsetTable[subsetTable$dNdS > 1,]

# scatter dS vs dN
ggplot(posTable, aes(x=dS, y=dN)) + 
  geom_point()

# scatter dN/dS vs dS
ggplot(posTable, aes(x=dS, y=dNdS)) + 
  geom_point()

# genes under purifying selection
negTable <- subsetTable[subsetTable$dNdS < 1,]

# scatter dS vs dN
ggplot(negTable, aes(x=dS, y=dN)) + 
  geom_point()

# scatter dN/dS vs dS
ggplot(negTable, aes(x=dS, y=dNdS)) + 
  geom_point()


# density plots
# http://www.sthda.com/english/articles/32-r-graphics-essentials/133-plot-one-variable-frequency-graph-density-distribution-and-more/

# create a basic plot
subsetPlot <- ggplot(subsetTable, aes(x = dNdS))

# histogram
subsetPlot + geom_histogram(bins = 30, color = "black", fill = "gray") 

