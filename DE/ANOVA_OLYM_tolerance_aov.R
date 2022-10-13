#!/usr/bin/env Rscript
# usage: Rscript ANOVA_OlympicsGenotypes_aov.r workingDir countsFile factorGroupingFile
# usage Ex: Rscript ANOVA_OlympicsGenotypes_aov.r /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCNA_DETolerance /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_OLYM_WGCNA/OLYM_60_eigengeneExpression_line.csv /Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_OlympicsTolerance.csv
# R script to perform statistical analysis of gene count tables using aov
# note: https://www.r-bloggers.com/2022/05/two-way-anova-example-in-r-quick-guide/

# turn off scientific notation
options(scipen=999)

# install packages
#install.packages("ggpubr")
#install.packages("multcomp")

# load libraries
library(ggpubr)
library(multcomp)

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
affyData <- merge(targets, inputTable, by = 'row.names') 

# convert treatment to a factor
affyData$treatment <- factor(affyData$treatment,
                             levels = c("VIS", "UV"),
                             labels = c("VIS", "UV"))
# convert tolerance to a factor
affyData$tolerance <- factor(affyData$tolerance,
                             levels = c("tol", "nTol"),
                             labels = c("tol", "nTol"))

# set row names
rownames(affyData) <- affyData$Row.names

# remove the Row.names column
affyData <- affyData[,2:4]

# check the structure of the expression data
#str(affyData)

# make frequency tables
#table(affyData$treatment, affyData$tolerance)

# retrieve module column name
modName <- colnames(affyData)[3]

# update module column name
colnames(affyData)[3] <- "expression"

# create a colored box plot
exportFile <- paste(modName, "coloredBoxPlot.png", sep="_")
png(file=exportFile)
ggboxplot(data=affyData, x="treatment", y="expression", color="tolerance",
          palette = c("#E69F00", "#009E73"))
dev.off()

# box plot with two variable factors
exportFile <- paste(modName, "twoVariableFactorBoxPlot.png", sep="_")
png(file=exportFile)
boxplot(expression ~ tolerance * treatment, data=affyData, frame=FALSE,
        col = c("#E69F00", "#009E73"))
dev.off()

# two way interaction plot
interaction.plot(x.factor = affyData$treatment, trace.factor = affyData$tolerance,
                 response = affyData$expression, fun = mean,
                 type = "b", legend = TRUE,
                 xlab = "Treatment", ylab="Expression",
                 pch=c(1,19), col = c("#E69F00", "#009E73"))

# compute two-way anova
affy.aov <- aov(expression ~ tolerance * treatment, data = affyData)

# output summary statistics
sumAffy <- summary(affy.aov)

# write summary statistics to a file
exportFile <- paste(modName, "ANOVA_summary.csv", sep="_")
capture.output(sumAffy, file=exportFile)

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
