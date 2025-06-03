#!/usr/bin/env Rscript

#Set working directory
#workingDir = args[1];
workingDir = "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Genotypes"
setwd(workingDir)

# Load the expression and trait data saved in the first part
#importFile <- paste(set, "dataInput.RData", sep="-")
importFile <- "OLYM-dataInput.RData"
lnames1 = load(file = importFile)

# Load network data saved in the second part
#importFile <- paste(tag, "networkConstruction-stepByStep.RData", sep="-")
importFile <- "OLYM_30-networkConstruction-stepByStep.RData"
lnames2 = load(file = importFile)
