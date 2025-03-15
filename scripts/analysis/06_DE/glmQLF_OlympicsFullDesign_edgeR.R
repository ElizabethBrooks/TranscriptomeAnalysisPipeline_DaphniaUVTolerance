#!/usr/bin/env Rscript
#Usage: Rscript glmQLF_OlympicsGenotypesTolerance_edgeR.r workingDir countsFile startColumn endColumn factorGroupingFile
# usage ex: Rscript glmQLF_OlympicsGenotypesTolerance_edgeR.r /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/ensembl/GCA_021134715.1/biostatistics/DEAnalysis/GenotypesTolerance /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/ensembl/GCA_021134715.1/biostatistics/geneCounts_merged_genome_counted_htseq_run1.txt 1 24 /Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_OlympicsGenotypesTolerance.csv
#Usage Ex: Rscript glmQLF_OlympicsGenotypesTolerance_edgeR.r /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/DEGenotypesTolerance /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/GeneCountsAnalyzed_KAP4/Formatted/cleaned.csv 1 24 /Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_OlympicsGenotypesTolerance.csv
#Usage Ex: Rscript glmQLF_OlympicsGenotypesTolerance_edgeR.r /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCNA_DEGenotypesTolerance /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_OLYM_WGCNA/OLYM_60_eigengeneExpression.csv 1 24 /Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_OlympicsGenotypesTolerance.csv
#R script to perform statistical analysis of gene count tables using edgeR GLM

# https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7873980/
# https://support.bioconductor.org/p/132926/
# https://support.bioconductor.org/p/106608/
# https://conjugateprior.org/2013/01/formulae-in-r-anova/

#Install edgeR and statmod, this should only need to be done once
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("edgeR")
#install.packages("statmod")

#Turn off scientific notation
options(scipen = 999)

#Load the edgeR library
library(edgeR)
library(statmod)
library(ggplot2)
library(ggrepel)
library(ggVennDiagram)
library(rcartocolor)

# Plotting Palettes
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
# https://github.com/Nowosad/rcartocolor
plotColors <- carto_pal(12, "Safe")
plotColorSubset <- c(plotColors[4], plotColors[5], plotColors[6])

#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)

#Set working directory
#workingDir = args[1]
workingDir="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/FullDesign"
setwd(workingDir)

#Import gene count data
#inputTable <- read.csv(file=args[2], row.names="gene")[ ,args[3]:args[4]]
inputTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/GeneCountsAnalyzed/Formatted/cleaned.csv", row.names="gene")[ ,1:24]

#Trim the data table
countsTable <- head(inputTable, - 5)

#Import grouping factor
#targets <- read.csv(file=args[5], row.names="sample")
targets <- read.csv(file="/Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_OlympicsFullDesign.csv", row.names="sample")

#Retrieve input FDR cutoff
#fdrCut=as.numeric(args[6])

#Setup a design matrix
group <- factor(paste(targets$treatment,targets$tolerance,targets$genotype,sep="."))
#cbind(targets,Group=group)
#Create DGE list object
list <- DGEList(counts=countsTable,group=group)
colnames(list) <- rownames(targets)

#Plot the library sizes before normalization
jpeg("glmQLF_plotBarsBefore.jpg")
barplot(list$samples$lib.size*1e-6, names=1:24, ylab="Library size (millions)")
dev.off()

#Retain genes only if it is expressed at a minimum level
keep <- filterByExpr(list)
summary(keep)
list <- list[keep, , keep.lib.sizes=FALSE]

#Use TMM normalization to eliminate composition biases
# between libraries
list <- calcNormFactors(list)
#list$samples
#Write normalized counts to file
normList <- cpm(list, normalized.lib.sizes=TRUE)
write.table(normList, file="glmQLF_normalizedCounts.csv", sep=",", row.names=TRUE, quote=FALSE)

#Write log transformed normalized counts to file
normListLog <- cpm(list, normalized.lib.sizes=TRUE, log=TRUE)
write.table(normListLog, file="glmQLF_normalizedCounts_logTransformed.csv", sep=",", row.names=TRUE, quote=FALSE)

#Verify TMM normalization using a MD plot
#Write plot to file
jpeg("glmQLF_plotMDBefore.jpg")
plotMD(cpm(list, log=TRUE), column=1)
abline(h=0, col=plotColorSubset[3], lty=2, lwd=2)
dev.off()

#Use a MDS plot to visualizes the differences
# between the expression profiles of different samples
points <- c(0,1,2,3,15,16,17,18)
colors <- rep(c(plotColors[4], plotColors[5], plotColors[6], plotColors[11]), 2)
#Write plot with legend to file
jpeg("glmQLF_plotMDS.jpg")
plotMDS(list, col=colors[group], pch=points[group])
legend("topleft", legend=levels(group), pch=points, col=colors, ncol=2)
dev.off()
#Write plot without legend to file
jpeg("glmQLF_plotMDS_noLegend.jpg")
plotMDS(list, col=colors[group], pch=points[group])
dev.off()

# Create a PCA plot with a legend
jpeg("glmQLF_plotPCA.jpg")
plotMDS(list, col=colors[group], pch=points[group], gene.selection="common")
legend("topleft", legend=levels(group), pch=points, col=colors, ncol=2)
dev.off()

# Create a PCA plot without a legend
jpeg("glmQLF_plotPCA_noLegend.jpg")
plotMDS(list, col=colors[group], pch=points[group], gene.selection="common")
dev.off()

##
#The experimental design is parametrized with a one-way layout, 
# where one coefficient is assigned to each group
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
#design

#Next, the NB dispersion is estimated
list <- estimateDisp(list, design, robust=TRUE)
#list$common.dispersion
#Visualize the dispersion estimates with a BCV plot
#Write plot to file
jpeg("glmQLF_plotBCV.jpg")
plotBCV(list)
dev.off()

#Now, estimate and plot the QL dispersions
fit <- glmQLFit(list, design, robust=TRUE)
#head(fit$coefficients)
#Write plot to file
jpeg("glmQLF_plotQLDisp.jpg")
plotQLDisp(fit)
dev.off()

# view column order
colnames(fit)


# testing constrasts
# treatment
con.treatment.test <- c(c(1,1,1,1,-1,-1,-1,-1))
treat.anov.test <- glmTreat(fit, contrast=con.treatment.test, lfc=log2(1.2))
summary(decideTests(treat.anov.test))
# tolerance
con.tolerance.test <- c(c(-1,-1,1,1,-1,-1,1,1))
tol.anov.test <- glmTreat(fit, contrast=con.tolerance.test, lfc=log2(1.2))
summary(decideTests(tol.anov.test))
# interaction
con.interaction.test <- c(c(-1,-1,1,1),-c(-1,-1,1,1))
int.anov.test <- glmTreat(fit, contrast=con.interaction.test, lfc=log2(1.2))
summary(decideTests(int.anov.test))

# testing nested (averaged) constrasts
# treatment
con.treatment.nest <- c(c(0.5,0.5,0.5,0.5,-0.5,-0.5,-0.5,-0.5))
treat.anov.nest <- glmTreat(fit, contrast=con.treatment.nest, lfc=log2(1.2))
summary(decideTests(treat.anov.nest))
# tolerance
con.tolerance.nest <- c(c(-0.5,-0.5,0.5,0.5,-0.5,-0.5,0.5,0.5))
tol.anov.nest <- glmTreat(fit, contrast=con.tolerance.nest, lfc=log2(1.2))
summary(decideTests(tol.anov.nest))
# interaction
con.interaction.nest <- c(c(-0.5,-0.5,0.5,0.5),-c(-0.5,-0.5,0.5,0.5))
int.anov.nest <- glmTreat(fit, contrast=con.interaction.nest, lfc=log2(1.2))
summary(decideTests(int.anov.nest))

# testing explicit contrasts
con.all.test <- makeContrasts(treatment = (UV.NTol.E05 - UV.NTol.R2 - UV.Tol.Y023 - UV.Tol.Y05)
                               - (VIS.NTol.E05 - VIS.NTol.R2 - VIS.Tol.Y023 - VIS.Tol.Y05),
                              tolerance = (UV.NTol.E05 - UV.NTol.R2 - VIS.NTol.E05 - VIS.NTol.R2)
                               - (UV.Tol.Y023 - UV.Tol.Y05 - VIS.Tol.Y023 - VIS.Tol.Y05),
                              interaction = ((UV.NTol.E05 - UV.NTol.R2 - UV.Tol.Y023 - UV.Tol.Y05)
                               - (VIS.NTol.E05 - VIS.NTol.R2 - VIS.Tol.Y023 - VIS.Tol.Y05))
                               - ((UV.NTol.E05 - UV.NTol.R2 - VIS.NTol.E05 - VIS.NTol.R2)
                               - (UV.Tol.Y023 - UV.Tol.Y05 - VIS.Tol.Y023 - VIS.Tol.Y05)),
                               levels=design)
# treatment
treat.anov <- glmTreat(fit, contrast=con.all.test[,"treatment"], lfc=log2(1.2))
summary(decideTests(treat.anov))
# tolerance
tol.anov <- glmTreat(fit, contrast=con.all.test[,"tolerance"], lfc=log2(1.2))
summary(decideTests(tol.anov))
# interaction
int.anov <- glmTreat(fit, contrast=con.all.test[,"interaction"], lfc=log2(1.2))
summary(decideTests(int.anov))

# testing explicit nested contrasts
con.all.nest <- makeContrasts(treatment = (UV.NTol.E05 + UV.NTol.R2 + UV.Tol.Y023 + UV.Tol.Y05)/4
                               - (VIS.NTol.E05 + VIS.NTol.R2 + VIS.Tol.Y023 + VIS.Tol.Y05)/4,
                              tolerance = (UV.NTol.E05 + UV.NTol.R2 + VIS.NTol.E05 + VIS.NTol.R2)/4
                               - (UV.Tol.Y023 + UV.Tol.Y05 + VIS.Tol.Y023 + VIS.Tol.Y05)/4,
                              interaction = ((UV.NTol.E05 + UV.NTol.R2 + UV.Tol.Y023 + UV.Tol.Y05)/4
                               - (VIS.NTol.E05 + VIS.NTol.R2 + VIS.Tol.Y023 + VIS.Tol.Y05)/4)
                               - ((UV.NTol.E05 + UV.NTol.R2 + VIS.NTol.E05 + VIS.NTol.R2)/4
                               - (UV.Tol.Y023 + UV.Tol.Y05 + VIS.Tol.Y023 + VIS.Tol.Y05)/4),
                              levels=design)
# treatment
treat.anov <- glmTreat(fit, contrast=con.all.nest[,"treatment"], lfc=log2(1.2))
summary(decideTests(treat.anov))
# tolerance
tol.anov <- glmTreat(fit, contrast=con.all.nest[,"tolerance"], lfc=log2(1.2))
summary(decideTests(tol.anov))
# interaction
int.anov <- glmTreat(fit, contrast=con.all.nest[,"interaction"], lfc=log2(1.2))
summary(decideTests(int.anov))
# test interaction
test.int.anov <- glmTreat(fit, contrast=c(con.all.nest[,"treatment"]-con.all.nest[,"tolerance"]), lfc=log2(1.2))
summary(decideTests(test.int.anov))


# testing nested interaction formulas
# requires specifying a baseline genotype
# first, make sure that the placebo is the reference level
# targets$treatment <- relevel(factor(targets$treatment), ref="VIS")
# then form the design matrix
# http://www.sthda.com/english/wiki/two-way-anova-test-in-r
# https://conjugateprior.org/2013/01/formulae-in-r-anova/
# https://support.bioconductor.org/p/68092/
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7873980/
# setup the design matrix
# invalid designs
#design.test <- model.matrix(~ (treatment * tolerance/genotype) + Error(replicate), data=targets)
#design.test <- model.matrix(~ treatment * tolerance/genotype, data=targets)
#design.test <- model.matrix(~ treatment * (tolerance/genotype), data=targets)
#design.test <- model.matrix(~ treatment * (tolerance + tolerance:genotype), data=targets)
#design.test <- model.matrix(~ (treatment * tolerance/genotype) + (tolerance + tolerance:genotype), data=targets)
#design.test <- model.matrix(~ (treatment * genotype) + (tolerance + tolerance:genotype), data=targets)
#design.test <- model.matrix(~ tolerance/(treatment * genotype), data=targets)
#design.test <- model.matrix(~ tolerance + tolerance:(treatment * genotype), data=targets)
#design.test <- model.matrix(~ tolerance + tolerance:(treatment + genotype + treatment:genotype), data=targets)
#design.test <- model.matrix(~ treatment * tolerance/genotype, data=targets)
#design.test <- model.matrix(~ (treatment + tolerance + treatment:tolerance) + (tolerance + tolerance:genotype), data=targets)
#design.test <- model.matrix(~ (treatment * tolerance) + (tolerance + tolerance:genotype), data=targets)
#design.test <- model.matrix(~ 0 + treatment + tolerance + genotype + replicate, data=targets)
#design.test <- model.matrix(~ 0 + treatment + (tolerance/genotype), data=targets)
#design.test <- model.matrix(~ 0 + treatment + tolerance + genotype, data=targets)
# valid designs
#design.test <- model.matrix(~ treatment * tolerance, data=targets)
#design.test <- model.matrix(~ treatment + tolerance + treatment:tolerance, data=targets)
#design.test <- model.matrix(~ treatment * genotype, data=targets)
#design.test <- model.matrix(~ 0 + treatment + tolerance + replicate, data=targets)
#design.test <- model.matrix(~ 0 + genotype + treatment + genotype:treatment, data=targets)
#design.test <- model.matrix(~ 0 + genotype + treatment, data=targets)

#design.test <- model.matrix(~ 0 + tolerance * treatment, data=targets)

# fit the design
#fit.test <- glmQLFit(list, design.test, robust=TRUE)
# view the coefficient names 
#colnames(fit.test)
# treatment (UV vs VIS)
#glm.treat.test <- glmTreat(fit.test, coef=4, lfc=log2(1.2))
#summary(decideTests(glm.treat.test))


#Test whether the average across all treatment groups is equal to the average across
#all VIS groups, to examine the overall effect of treatment
con.treatment <- makeContrasts(treatment = (UV.E05 + UV.R2 + UV.Y023 + UV.Y05)/4
  - (VIS.E05 + VIS.R2 + VIS.Y023 + VIS.Y05)/4,
  levels=design)

#Look at genes expressed across all treatment groups using QL F-test
#test.anov.treatment <- glmQLFTest(fit, contrast=con.treatment)
#summary(decideTests(test.anov.treatment))
#Write plot to file
#jpeg("glmQLF_2WayANOVA_treatment_plotMD.jpg")
#plotMD(test.anov.treatment)
#abline(h=c(-1, 1), col="blue")
#dev.off()
#Write tags table of DE genes to file
#tagsTblANOVA <- topTags(test.anov.treatment, n=nrow(test.anov.treatment$table), adjust.method="fdr")$table
#tagsTblANOVA.keep <- tagsTblANOVA$FDR <= fdrCut
#tagsTblANOVA.out <- tagsTblANOVA[tagsTblANOVA.keep,]
#write.table(tagsTblANOVA, file="glmQLF_2WayANOVA_treatment_topTags.csv", sep=",", row.names=TRUE, quote=FALSE)

#Look at genes with significant expression across all treatment groups
treat.anov.treatment <- glmTreat(fit, contrast=con.treatment, lfc=log2(1.2))
summary(decideTests(treat.anov.treatment))
#Write plot to file
jpeg("glmQLF_2WayANOVA_treatment_plotMD_LFC1.2.jpg")
plotMD(treat.anov.treatment)
abline(h=c(-1, 1), col="blue")
dev.off()
#Write tags table of DE genes to file
tagsTblANOVA.filtered <- topTags(treat.anov.treatment, n=nrow(treat.anov.treatment$table), adjust.method="fdr")$table
#tagsTblANOVA.filtered.keep <- tagsTblANOVA.filtered$FDR <= fdrCut
#tagsTblANOVA.filtered.out <- tagsTblANOVA.filtered[tagsTblANOVA.filtered.keep,]
write.table(tagsTblANOVA.filtered, file="glmQLF_2WayANOVA_treatment_topTags_LFC1.2.csv", sep=",", row.names=TRUE, quote=FALSE)


#Test whether the average across all tolerant groups is equal to the average across
#all not tolerant groups, to examine the overall effect of tolerance
con.tolerance <- makeContrasts(tolerance = (UV.Y05 + VIS.Y05 + UV.Y023 + VIS.Y023)/4
  - (UV.E05 + VIS.E05 + UV.R2 + VIS.R2)/4,
  levels=design)

#Look at genes expressed across all tolerance groups using QL F-test
#test.anov.tolerance <- glmQLFTest(fit, contrast=con.tolerance)
#summary(decideTests(test.anov.tolerance))
#Write plot to file
#jpeg("glmQLF_2WayANOVA_tolerance_plotMD.jpg")
#plotMD(test.anov.tolerance)
#abline(h=c(-1, 1), col="blue")
#dev.off()
#Write tags table of DE genes to file
#tagsTblANOVAtolerance <- topTags(test.anov.tolerance, n=nrow(test.anov.tolerance$table), adjust.method="fdr")$table
#tagsTblANOVAtolerance.keep <- tagsTblANOVAtolerance$FDR <= fdrCut
#tagsTblANOVAtolerance.out <- tagsTblANOVAtolerance[tagsTblANOVAtolerance.keep,]
#write.table(tagsTblANOVAtolerance, file="glmQLF_2WayANOVA_tolerance_topTags.csv", sep=",", row.names=TRUE, quote=FALSE)

#Look at genes with significant expression across all tolerance groups
treat.anov.tolerance <- glmTreat(fit, contrast=con.tolerance, lfc=log2(1.2))
summary(decideTests(treat.anov.tolerance))
#Write plot to file
jpeg("glmQLF_2WayANOVA_tolerance_plotMD_LFC1.2.jpg")
plotMD(treat.anov.tolerance)
abline(h=c(-1, 1), col="blue")
dev.off()
#Write tags table of DE genes to file
tagsTblANOVAtolerance.filtered <- topTags(treat.anov.tolerance, n=nrow(treat.anov.tolerance$table), adjust.method="fdr")$table
#tagsTblANOVAtolerance.filtered.keep <- tagsTblANOVAtolerance.filtered$FDR <= fdrCut
#tagsTblANOVAtolerance.filtered.out <- tagsTblANOVAtolerance.filtered[tagsTblANOVAtolerance.filtered.keep,]
write.table(tagsTblANOVAtolerance.filtered, file="glmQLF_2WayANOVA_tolerance_topTags_LFC1.2.csv", sep=",", row.names=TRUE, quote=FALSE)


#Test whether there is an interaction effect
con.Inter <- makeContrasts(Inter = ((UV.E05 + UV.R2 + UV.Y023 + UV.Y05)/4
  - (VIS.E05 + VIS.R2 + VIS.Y023 + VIS.Y05)/4)
  - ((UV.Y05 + VIS.Y05 + UV.Y023 + VIS.Y023)/4
  - (UV.E05 + VIS.E05 + UV.R2 + VIS.R2)/4),
  levels=design)

#Look at genes expressed across all interaction groups using QL F-test
#test.anov.Inter <- glmQLFTest(fit, contrast=con.Inter)
#summary(decideTests(test.anov.Inter))
#Write plot to file
#jpeg("glmQLF_2WayANOVA_interaction_plotMD.jpg")
#plotMD(test.anov.Inter)
#abline(h=c(-1, 1), col="blue")
#dev.off()
#Write tags table of DE genes to file
#tagsTblANOVAInter <- topTags(test.anov.Inter, n=nrow(test.anov.Inter$table), adjust.method="fdr")$table
#tagsTblANOVAInter.keep <- tagsTblANOVAInter$FDR <= fdrCut
#tagsTblANOVAInter.out <- tagsTblANOVAInter[tagsTblANOVAInter.keep,]
#write.table(tagsTblANOVAInter, file="glmQLF_2WayANOVA_interaction_topTags.csv", sep=",", row.names=TRUE, quote=FALSE)

#Look at genes with significant expression
treat.anov.Inter <- glmTreat(fit, contrast=con.Inter, lfc=log2(1.2))
summary(decideTests(treat.anov.Inter))
#Write plot to file
jpeg("glmQLF_2WayANOVA_interaction_plotMD_LFC1.2.jpg")
plotMD(treat.anov.Inter)
abline(h=c(-1, 1), col="blue")
dev.off()
#Generate table of DE genes
tagsTblANOVAInter.filtered <- topTags(treat.anov.Inter, n=nrow(treat.anov.Inter$table), adjust.method="fdr")$table
#tagsTblANOVAInter.filtered.keep <- tagsTblANOVAInter.filtered$FDR <= fdrCut
#tagsTblANOVAInter.filtered.out <- tagsTblANOVAInter.filtered[tagsTblANOVAInter.filtered.keep,]
write.table(tagsTblANOVAInter.filtered, file="glmQLF_2WayANOVA_interaction_topTags_LFC1.2.csv", sep=",", row.names=TRUE, quote=FALSE)
