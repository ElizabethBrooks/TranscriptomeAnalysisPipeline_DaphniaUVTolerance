#!/usr/bin/env Rscript

# usage: Rscript ANOVA_OLYMGenotypes_aov.r workingDir countsFile factorGroupingFile
# usage Ex: Rscript ANOVA_OLYMGenotypes_aov.r /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/ensembl/GCA_021134715.1/biostatistics/NetworkAnalysis/WGCN_genotypes_WGCNA /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_OLYM_WGCNA/OLYM_60_eigengeneExpression_line.csv /Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVGenotypes/InputData/expDesign_OlympicsGenotypes.csv
# R script to perform statistical analysis of gene count tables using aov
# note: https://www.r-bloggers.com/2022/05/two-way-anova-example-in-r-quick-guide/
# https://conjugateprior.org/2013/01/formulae-in-r-anova/
# https://environmentalcomputing.net/statistics/linear-models/anova/anova-nested/

# turn off scientific notation
options(scipen=999)

# install packages
#install.packages("ggpubr")
#install.packages("multcomp")

# load libraries
library(ggpubr)
library(multcomp)
library(rcartocolor)

# Plotting Palettes
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
# https://github.com/Nowosad/rcartocolor
plotColors <- carto_pal(12, "Safe")
plotColorSubset <- c(plotColors[11], plotColors[6], plotColors[4], plotColors[5])

# retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)

# set working directory
workingDir = args[1];
setwd(workingDir)

# import expression data
inputTable <- read.csv(args[2], row.names="gene")

# transpose expression data
inputTable <- as.data.frame(t(inputTable))

# import grouping factor
targets <- read.csv(file=args[3], row.names="sample")

# setup data frame
expData <- merge(targets, inputTable, by = 'row.names') 

# convert treatment to a factor
expData$treatment <- factor(expData$treatment,
                             levels = c("UV", "VIS"),
                             labels = c("UV", "VIS"))

# convert tolerance to a factor
expData$tolerance <- factor(expData$tolerance,
                            levels = c("Tol", "NTol"),
                            labels = c("Tol", "NTol"))

# convert genotype to a factor
expData$genotype <- factor(expData$genotype,
                             levels = c("Y05", "Y023", "E05", "R2"),
                             labels = c("Y05", "Y023", "E05", "R2"))

# convert replicate to a factor
#expData$replicate <- factor(expData$replicate,
#                             levels = c("1", "2", "3"),
#                             labels = c("1", "2", "3"))


# view data
#print(expData)


# set row names
rownames(expData) <- expData$Row.names

# remove the Row.names column
#expData <- expData[,2:6]
expData <- expData[,-1]

# remove the replicate column
expData <- expData[,-4]

# check the structure of the expression data
#str(expData)

# make frequency tables
#table(expData$treatment, expData$genotype, expData$tolerance, expData$replicate)

# retrieve module column name
#modName <- colnames(expData)[5]
modName <- colnames(expData)[4]

# update module column name
#colnames(expData)[5] <- "expression"
colnames(expData)[4] <- "expression"


# view data
#print(expData)


# compute two-way anova
#affy.aov <- aov(expression ~ tolerance/genotype * treatment + Error(replicate), data = expData)
#affy.aov <- aov(expression ~ treatment * tolerance/genotype, data = expData)
#affy.aov <- aov(expression ~ tolerance/genotype * treatment, data = expData)
affy.aov <- aov(expression ~ treatment * (tolerance/genotype), data = expData)

# output summary statistics
sumAffy <- summary(affy.aov)

# write summary statistics to a file
exportFile <- paste(modName, "ANOVA_summary.csv", sep="_")
capture.output(sumAffy, file=exportFile)


# create a colored box plot
exportFile <- paste(modName, "coloredBoxPlot.png", sep="_")
png(file=exportFile)
ggboxplot(data=expData, x="treatment", y="expression", color="genotype",
          palette = plotColorSubset)
dev.off()

# box plot with two variable factors
#boxplot(expression ~ tolerance/genotype * treatment + Error(replicate), data=expData, frame=FALSE,
#boxplot(expression ~ treatment * tolerance/genotype, data=expData, frame=FALSE,
#boxplot(expression ~ tolerance/genotype * treatment, data=expData, frame=FALSE,
#exportFile <- paste(modName, "twoVariableFactorBoxPlot.png", sep="_")
#png(file=exportFile)
#boxplot(expression ~ treatment * (tolerance/genotype), data=expData, frame=FALSE,
#        col = plotColorSubset)
#dev.off()

# two way interaction plot
# two way interaction plot
exportFile <- paste(modName, "twoWayInteractionPlot.png", sep="_")
png(file=exportFile)
interaction.plot(x.factor = expData$treatment, trace.factor = expData$genotype,
                 response = expData$expression, fun = mean,
                 type = "b", legend = TRUE,
                 xlab = "Treatment", ylab="Expression",
                 pch=c(1,19), col = plotColorSubset)
dev.off()


## check the validity of ANOVA assumptions
# The data must be regularly distributed, and the variation between groups must be homogeneous

## examine the assumption of homogeneity of variance
# the residuals versus fits graphic are used to assess for variance homogeneity

# plot homogeneity of variances
exportFile <- paste(modName, "homogeneityPlot.png", sep="_")
png(file=exportFile)
plot(affy.aov, 1)
dev.off()

# examine the assumption of normality
# in a residuals’ normality plot the residuals quantiles are displayed against the normal distribution quantiles
# the residuals’ normal probability plot is used to confirm that the residuals are normally distributed
# the residuals’ normal probability plot should roughly follow a straight line

# normality plot
exportFile <- paste(modName, "normalityPlot.png", sep="_")
png(file=exportFile)
plot(affy.aov, 2)
dev.off()

# extract the residuals
aovRes <- residuals(object = affy.aov)

# run Shapiro-Wilk test
swTest <- shapiro.test(x = aovRes)

# write test statistics to a file
exportFile <- paste(modName, "shapiroTest_summary.csv", sep="_")
capture.output(swTest, file=exportFile)
