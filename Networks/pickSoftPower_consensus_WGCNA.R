#!/usr/bin/env Rscript

# script to help pick a soft threshold values for a set of samples using WGNCA
# usage: Rscript pickSoftPower_consensus_WGCNA.R workingDir fileTag
# usage ex: Rscript pickSoftPower_consensus_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/ensembl/GCA_021134715.1/biostatistics/NetworkAnalysis/WGCN_genotypes_WGCNA genotype
# usage ex: Rscript pickSoftPower_consensus_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/ensembl/GCA_021134715.1/biostatistics/NetworkAnalysis/WGCN_tolerance_WGCNA tolerance

#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)

#Set working directory
workingDir = args[1];
#workingDir="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_genotype_WGCNA"
setwd(workingDir)

# retrieve file tage
fileTag <- args[2]

# load the WGCNA package
library(WGCNA)

# the following setting is important, do not omit.
options(stringsAsFactors = FALSE)

# Allow multi-threading within WGCNA. At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads()

# Load the data saved in the first part
importFile <- paste("Consensus-dataInput", fileTag, sep="-")
importFile <- paste(importFile, "RData", sep=".")
lnames = load(file = importFile)
#The variable lnames contains the names of loaded variables.
lnames
# Get the number of sets in the multiExpr structure.
nSets = checkSets(multiExpr)$nSets

# Choose a set of soft-thresholding powers
powers = c(seq(4,10,by=1), seq(12,36, by=2))
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets)
# Call the network topology analysis function for each set in turn
for (set in 1:nSets){
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, 
                                                     powerVector=powers,
                                                     verbose = 2)[[2]])
}
collectGarbage()

# Plot the results:
colors = c("black", "red", "blue", "green")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", 
             "Mean connectivity", 
             "Median connectivity",
             "Max connectivity")
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4)
for (set in 1:nSets){
  for (col in 1:length(plotCols)){
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE)
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE)
  }
}

# Plot the quantities in the chosen columns vs. the soft thresholding power
#sizeGrWindow(8, 6)
#par(mfcol = c(4,1))
#par(mar = c(4.2, 4.2 , 2.2, 0.5))
exportFile <- paste("ConsensusSoftPowers", fileTag, sep="_")
exportFile <- paste(exportFile, "png", sep=".")
png(file = exportFile, width = 12, height = 12, units="in", res=150)
par(mfrow=c(4,1))
par(mar = c(0, 4, 2, 0))
cex1 = 0.7
for (col in 1:length(plotCols)) for (set in 1:nSets){
  if (set==1){
    plot(powerTables[[set]]$data[,1], 
         -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",
         ylab=colNames[col],
         type="n", 
         ylim = ylim[, col],
         main = colNames[col])
    addGrid()
  }
  if (col==1){
    text(powerTables[[set]]$data[,1], 
         -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,
         cex=cex1,
         col=colors[set])
  } else
    text(powerTables[[set]]$data[,1], 
         powerTables[[set]]$data[,plotCols[col]],
         labels=powers,
         cex=cex1,
         col=colors[set])
  if (col==1){
    legend("bottomright", 
           legend = setLabels, 
           col = colors, 
           pch = 20)
  } else
    legend("topright", 
           legend = setLabels, 
           col = colors, 
           pch = 20)
}
dev.off()
