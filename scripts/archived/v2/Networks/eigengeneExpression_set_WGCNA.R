#!/usr/bin/env Rscript

# script to identify set-specific and consensos modules
# usage: Rscript eigengeneExpression_set_WGCNA.R workingDir subsetTag minModuleSize
# usage ex: Rscript eigengeneExpression_set_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/ensembl/GCA_021134715.1/biostatistics/NetworkAnalysis/WGCN_tolerance_WGCNA OLYM 60
# usage ex: Rscript eigengeneExpression_set_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/ensembl/GCA_021134715.1/biostatistics/NetworkAnalysis/WGCN_tolerance_WGCNA tol 60
# usage ex: Rscript eigengeneExpression_set_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/ensembl/GCA_021134715.1/biostatistics/NetworkAnalysis/WGCN_tolerance_WGCNA nTol 60
# alternate usage ex: Rscript eigengeneExpression_set_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/ensembl/GCA_021134715.1/biostatistics/NetworkAnalysis/WGCN_tolerance_WGCNA Y05 60
# alternate usage ex: Rscript eigengeneExpression_set_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/ensembl/GCA_021134715.1/biostatistics/NetworkAnalysis/WGCN_tolerance_WGCNA Y023 60
# alternate usage ex: Rscript eigengeneExpression_set_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/ensembl/GCA_021134715.1/biostatistics/NetworkAnalysis/WGCN_tolerance_WGCNA E05 60
# alternate usage ex: Rscript eigengeneExpression_set_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/ensembl/GCA_021134715.1/biostatistics/NetworkAnalysis/WGCN_tolerance_WGCNA R2 60

#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)

#Set working directory
workingDir = args[1];
#workingDir="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/ensembl/GCA_021134715.1/biostatistics/NetworkAnalysis/WGCN_tolerance_WGCNA"
setwd(workingDir)

# load the WGCNA package
library(WGCNA)

# the following setting is important, do not omit.
options(stringsAsFactors = FALSE)

# retrieve subset tag
subsetTag <- args[2]
#subsetTag <- "Y05"

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = as.numeric(args[3])
#minModuleSize = 60

# load the subset data and rename them so that 
# they do not conflict with the consensus data
importFile <- paste(subsetTag, minModuleSize, sep="_")
importFile <- paste(importFile, "networkConstruction-stepByStep.RData", sep="-")
lnames = load(file = importFile)
lnames

# transpose each subset
datExpr0 = data = as.data.frame(t(MEs))
names(datExpr0) = rownames(MEs)
rownames(datExpr0) = names(MEs)

# add a column for the row names
datExpr0 <- cbind(gene = rownames(datExpr0), datExpr0)
rownames(datExpr0) <- NULL

# export the expression data as a csv file
exportFile <- paste(subsetTag, minModuleSize, sep="_")
exportFile <- paste(exportFile, "eigengeneExpression.csv", sep="_")
write.table(datExpr0, file=exportFile, sep=",", row.names=FALSE)
