#!/usr/bin/env Rscript
#Usage: Rscript glm_edgeR.r countsFile startColumn endColumn factorGroupingFile
#Usage Ex: Rscript glm_edgeR.r cleaned.csv 1 24 expDesign_Olympics_GRP1.csv
#R script to perform statistical analysis of gene count tables using edgeR two way ANOVA

countsTable <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/GeneCountsAnalyzed/genome_sortedName_samtoolsHisat2_run2_counted_htseq_run1/cleaned.csv", row.names="gene")[ ,1:24]
targets <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_Olympics_GRP1.csv", row.names="sample")

#Install edgeR, this should only need to be done once
#Since edgeR is already installed on the CRC this can be skipped if using the module
#bioLite("edgeR")
#install.packages("statmod")

#Load the edgeR library
library("edgeR")
library("statmod")

#Retrieve input file name of gene counts
#args = commandArgs(trailingOnly=TRUE)
#Test if there is one input argument
#if (length(args)!=3) {
#  stop("One file name and a range of columns must be supplied.n", call.=FALSE)
#}

#Import gene count data
#countsTable <- read.csv(file=args[1], row.names="gene")[ ,args[2]:args[3]]
#head(countsTable)
#Import grouping factor
#targets <- read.csv(file=args[4], row.names="sample")

#Setup a design matrix
group <- factor(paste(targets$treatment,targets$genotype,sep="."))
#cbind(targets,Group=group)
#Create DGE list object
list <- DGEList(counts=countsTable,group=group)
colnames(list) <- targets$sample

#Retain genes only if it is expressed at a minimum level
keep <- filterByExpr(list)
summary(keep)
list <- list[keep, , keep.lib.sizes=FALSE]

#Use TMM normalization to eliminate composition biases
# between libraries
list <- calcNormFactors(list)
#list$samples

#Verify TMM normalization using a MD plot
plotMD(cpm(list, log=TRUE), column=1)
abline(h=0, col="red", lty=2, lwd=2)

#Use a MDS plot to visualizes the differences
# between the expression profiles of different samples
points <- c(0,1,2,3,15,16,17,18)
colors <- rep(c("blue", "darkgreen", "red", "black"), 2)
plotMDS(list, col=colors[group], pch=points[group])
legend("topleft", legend=levels(group), pch=points, col=colors, ncol=2)
plotMDS(list, col=colors[group], pch=points[group])

#The experimental design is parametrized with a one-way layout, 
# where one coefficient is assigned to each group
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
#design

#Next, the NB dispersion is estimated
list <- estimateDisp(list, design, robust=TRUE)
list$common.dispersion
#Visualize the dispersion estimates with a BCV plot
plotBCV(list)

#Now, estimate and plot the QL dispersions
fit <- glmQLFit(list, design, robust=TRUE)
#head(fit$coefficients)
plotQLDisp(fit)

#Define a matrix of contrasts, where each column
# represents a contrast between two groups of interest
con.allPairs <- makeContrasts(
	#Pairwise
	E05.UVvsVIS = UV.E05 - VIS.E05,
	R2.UVvsVIS = UV.R2 - VIS.R2,
	Y023.UVvsVIS = UV.Y023 - VIS.Y023,
	Y05.UVvsVIS = UV.Y05 - VIS.Y05,
	levels=design)

#All pairs
test.allPairs <- glmLRT(fit, contrast=con.allPairs)
#topTags(test.allPairs)
summary(decideTests(test.allPairs))
plotMD(test.allPairs)
abline(h=c(-1, 1), col="blue")

#Pairwise E05.UVvsVIS
#qlf.E05.UVvsVIS <- glmQLFTest(fit, contrast=con[,"E05.UVvsVIS"])
#Test whether the differential expression is significant
treat.E05.UVvsVIS <- glmTreat(fit, contrast=con.allPairs[,"E05.UVvsVIS"], lfc=log2(1.2))
#topTags(treat.E05.UVvsVIS)
summary(decideTests(treat.E05.UVvsVIS))
plotMD(treat.E05.UVvsVIS)
abline(h=c(-1, 1), col="blue")

#Pairwise R2.UVvsVIS
#Test whether the differential expression is significant
treat.R2.UVvsVIS <- glmTreat(fit, contrast=con.allPairs[,"R2.UVvsVIS"], lfc=log2(1.2))
#topTags(treat.R2.UVvsVIS)
summary(decideTests(treat.R2.UVvsVIS))
plotMD(treat.R2.UVvsVIS)
abline(h=c(-1, 1), col="blue")

#Pairwise Y023.UVvsVIS
#Test whether the differential expression is significant
treat.Y023.UVvsVIS <- glmTreat(fit, contrast=con.allPairs[,"Y023.UVvsVIS"], lfc=log2(1.2))
#topTags(treat.Y023.UVvsVIS)
summary(decideTests(treat.Y023.UVvsVIS))
plotMD(treat.Y023.UVvsVIS)
abline(h=c(-1, 1), col="blue")

#Pairwise Y05.UVvsVIS
#Test whether the differential expression is significant
treat.Y05.UVvsVIS <- glmTreat(fit, contrast=con.allPairs[,"Y05.UVvsVIS"], lfc=log2(1.2))
#topTags(treat.Y05.UVvsVIS)
summary(decideTests(treat.Y05.UVvsVIS))
plotMD(treat.Y05.UVvsVIS)
abline(h=c(-1, 1), col="blue")

#ANOVA like comparisons of UV
anov.UV <- makeContrasts(UV.R2 - UV.E05,
  UV.Y023 - UV.E05,
  UV.Y05 - UV.E05,
  levels=design)
test.anov.UV <- glmLRT(fit, contrast=anov.UV)
#topTags(test.anov.UV)
summary(decideTests(test.anov.UV))
plotMD(test.anov.UV)
abline(h=c(-1, 1), col="blue")

#ANOVA like comparisons of VIS
anov.VIS <- makeContrasts(VIS.R2 - VIS.E05,
  VIS.Y023 - VIS.E05,
  VIS.Y05 - VIS.E05,
  levels=design)
test.anov.VIS <- glmLRT(fit, contrast=anov.VIS)
#topTags(test.anov.VIS)
summary(decideTests(test.anov.VIS))
plotMD(test.anov.VIS)
abline(h=c(-1, 1), col="blue")

#Test whether the average across all UV groups is equal to the average across
#all VIS groups, to examine the overall effect of treatment
con.UVvsVIS <- makeContrasts(UVvsVIS = (UV.E05 + UV.R2 + UV.Y023 + UV.Y05)/4
  - (VIS.E05 + VIS.R2 + VIS.Y023 + VIS.Y05)/4,
  levels=design)
test.anov.VIS <- glmLRT(fit, contrast=con.UVvsVIS)
#topTags(test.anov.VIS)
summary(decideTests(test.anov.VIS))
plotMD(test.anov.VIS)
abline(h=c(-1, 1), col="blue")

#Look at genes with significant expression across all UV groups
treat.anov.VIS <- glmTreat(fit, contrast=con.UVvsVIS, lfc=log2(1.2))
#topTags(treat.anov.VIS)
summary(decideTests(treat.anov.VIS))
plotMD(treat.anov.VIS)
abline(h=c(-1, 1), col="blue")
