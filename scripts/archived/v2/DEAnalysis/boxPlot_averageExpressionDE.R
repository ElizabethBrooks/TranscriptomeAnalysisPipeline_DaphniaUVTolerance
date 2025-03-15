#!/usr/bin/env Rscript

#if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
#BiocManager::install('_______')

# load libraries
library(ggpubr)
#library(multcomp)
library(rcartocolor)

# Plotting Palettes
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
# https://github.com/Nowosad/rcartocolor
plotColors <- carto_pal(12, "Safe")
plotColorSubset <- c(plotColors[11], plotColors[6], plotColors[4], plotColors[5])

#Set working directory
workingDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Genotypes"
setwd(workingDir)

# import normalized gene count data for the Olympics
inFile <- "glmQLF_normalizedCounts_logTransformed.csv"
normList <- read.csv(file=inFile, row.names="gene")

#Subset the treatment and control samples
#SETVIS <- normList[,grepl("VIS",colnames(normList))]
#SETUV <- normList[,grepl("UV",colnames(normList))]
#Subset the tolerant and non-tolerant samples
SETY05 <- normList[,grepl("Y05",colnames(normList))]
SETY05_VIS <- SETY05[,grepl("VIS",colnames(SETY05))]
SETY05_UV <- SETY05[,grepl("UV",colnames(SETY05))]
SETY023 <- normList[,grepl("Y023",colnames(normList))]
SETY023_VIS <- SETY023[,grepl("VIS",colnames(SETY023))]
SETY023_UV <- SETY023[,grepl("UV",colnames(SETY023))]
SETE05 <- normList[,grepl("E05",colnames(normList))]
SETE05_VIS <- SETE05[,grepl("VIS",colnames(SETE05))]
SETE05_UV <- SETE05[,grepl("UV",colnames(SETE05))]
SETR2 <- normList[,grepl("R2",colnames(normList))]
SETR2_VIS <- SETR2[,grepl("VIS",colnames(SETR2))]
SETR2_UV <- SETR2[,grepl("UV",colnames(SETR2))]

#Calculate the means of each gene across samples
#meanVIS <- cbind(rowMeans(as.matrix(SETVIS)), rep("VIS", length(SETVIS)))
#meanUV <- cbind(rowMeans(as.matrix(SETUV)), rep("UV", length(SETUV)))
#Calculate the means of each gene across samples
meanY05_VIS <- data.frame(expression = rowMeans(as.matrix(SETY05_VIS)), 
                          genotype = rep("Y05", nrow(SETY05_VIS)), 
                          treatment = rep("VIS", nrow(SETY05_VIS)))
meanY05_UV <- data.frame(expression = rowMeans(as.matrix(SETY05_UV)), 
                         genotype = rep("Y05", nrow(SETY05_UV)), 
                         treatment = rep("UV", nrow(SETY05_UV)))
meanY023_VIS <- data.frame(expression = rowMeans(as.matrix(SETY023_VIS)), 
                           genotype = rep("Y023", nrow(SETY023_VIS)), 
                           treatment = rep("VIS", nrow(SETY023_VIS)))
meanY023_UV <- data.frame(expression = rowMeans(as.matrix(SETY023_UV)), 
                          genotype = rep("Y023", nrow(SETY023_UV)), 
                          treatment = rep("UV", nrow(SETY023_UV)))
meanE05_VIS <- data.frame(expression = rowMeans(as.matrix(SETE05_VIS)), 
                          genotype = rep("E05", nrow(SETE05_VIS)), 
                          treatment = rep("VIS", nrow(SETE05_VIS)))
meanE05_UV <- data.frame(expression = rowMeans(as.matrix(SETE05_UV)), 
                         genotype = rep("E05", nrow(SETE05_UV)), 
                         treatment = rep("UV", nrow(SETE05_UV)))
meanR2_VIS <- data.frame(expression = rowMeans(as.matrix(SETR2_VIS)), 
                         genotype = rep("R2", nrow(SETR2_VIS)), 
                         treatment = rep("VIS", nrow(SETR2_VIS)))
meanR2_UV <- data.frame(expression = rowMeans(as.matrix(SETR2_UV)), 
                        genotype = rep("R2", nrow(SETR2_UV)), 
                        treatment = rep("UV", nrow(SETR2_UV)))

# combine tables
meanList <- rbind(meanY05_VIS, meanY023_VIS, meanE05_VIS, meanR2_VIS, meanY05_UV, meanY023_UV, meanE05_UV, meanR2_UV)

# create a colored box plot
#exportFile <- "treatment_coloredBoxPlot.png"
#png(file=exportFile)
ggboxplot(data=meanList, x="treatment", y="expression", color="genotype",
          palette = plotColorSubset)
#dev.off()

