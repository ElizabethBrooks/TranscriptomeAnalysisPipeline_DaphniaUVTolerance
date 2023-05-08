#!/usr/bin/env Rscript
#Usage: Rscript glmQLF_edgeR.r workingDir countsFile startColumn endColumn factorGroupingFile
#Usage Ex: Rscript glmQLF_edgeR.r /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/DEGenotypes cleaned.csv 1 24 expDesign_OlympicsTolerance.csv
#R script to perform statistical analysis of gene count tables using edgeR GLM

#Install edgeR and statmod, this should only need to be done once
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("edgeR")
#install.packages("statmod")
#install.packages("ggrepel")

#Turn off scientific notation
options(scipen = 999)

#Load the edgeR library
library(edgeR)
library(statmod)
library(ghibli)
library(ggplot2)
library(ggrepel)
library(ggVennDiagram)

#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)

#Set working directory
#workingDir = args[1];
workingDir="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Tolerance"
#workingDir="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCNA_DEGenotypes"
setwd(workingDir)

#Import gene count data
#inputTable <- read.csv(file=args[2], row.names="gene")[ ,args[3]:args[4]]
inputTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/GeneCountsAnalyzed/Formatted/cleaned.csv", row.names="gene")[ ,1:24]
#inputTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_OLYM_WGCNA/OLYM_60_eigengeneExpression.csv", row.names="gene")[ ,1:24]

#Trim the data table
countsTable <- head(inputTable, - 5)

#Import grouping factor
#targets <- read.csv(file=args[5], row.names="sample")
targets <- read.csv(file="/Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_OlympicsTolerance.csv", row.names="sample")
#targets <- read.csv(file="/Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_OlympicsTolerance.csv", row.names="sample")

#Retrieve input FDR cutoff
#fdrCut=as.numeric(args[5])

# Plotting Palettes
# retrieve the vector of colors associated with PonyoMedium
ghibli_colors <- ghibli_palette("PonyoMedium", type = "discrete")
# vector with a subset of colors associated with PonyoMedium
ghibli_subset <- c(ghibli_colors[3], ghibli_colors[6], ghibli_colors[4])

#Setup a design matrix
group <- factor(paste(targets$treatment,targets$tolerance,sep="."))
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
abline(h=0, col=ghibli_colors[4], lty=2, lwd=2)
dev.off()

#Use a MDS plot to visualizes the differences
# between the expression profiles of different samples
points <- c(0,1,15,16)
colors <- rep(c(ghibli_colors[3], ghibli_colors[4]), 2)
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


#Test whether the average across all UV groups is equal to the average across
#all VIS groups, to examine the overall effect of treatment
con.treatment <- makeContrasts(treatment = (UV.NTol + UV.Tol)/4
  - (VIS.NTol + VIS.Tol)/4,
  levels=design)

#Look at genes with significant expression across all UV groups
treat.anov.treatment <- glmTreat(fit, contrast=con.treatment, lfc=log2(1.2))
summary(decideTests(treat.anov.treatment))
#Write plot to file
jpeg("glmQLF_2WayANOVA_treatment_plotMD_LFC1.2.jpg")
plotMD(treat.anov.treatment)
abline(h=c(-1, 1), col=ghibli_colors[3])
dev.off()
#Write tags table of DE genes to file
tagsTblANOVATreatment <- topTags(treat.anov.treatment, n=nrow(treat.anov.treatment$table), adjust.method="fdr")$table
write.table(tagsTblANOVATreatment, file="glmQLF_2WayANOVA_treatment_topTags_LFC1.2.csv", sep=",", row.names=TRUE, quote=FALSE)

# add column for identifying direction of DE gene expression
tagsTblANOVATreatment$topDE <- "NA"
# identify significantly up DE genes
tagsTblANOVATreatment$topDE[tagsTblANOVATreatment$logFC > 1 & tagsTblANOVATreatment$FDR < 0.05] <- "UP"
# identify significantly down DE genes
tagsTblANOVATreatment$topDE[tagsTblANOVATreatment$logFC < -1 & tagsTblANOVATreatment$FDR < 0.05] <- "DOWN"
# create volcano plot
jpeg("glmQLF_2WayANOVA_treatment_volcano_LFC1.2.jpg")
ggplot(data=tagsTblANOVATreatment, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset, breaks = c("Up", "Down"))
dev.off()
# create volcano plot with labels
labelSetTreatment <- tagsTblANOVATreatment[tagsTblANOVATreatment$topDE == "UP" | tagsTblANOVATreatment$topDE == "DOWN",]
jpeg("glmQLF_2WayANOVA_treatment_volcanoLabeled_LFC1.2.jpg")
ggplot(data=tagsTblANOVATreatment, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  ggrepel::geom_text_repel(data = labelSetTreatment, aes(label = row.names(labelSetTreatment))) +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset, breaks = c("Up", "Down"))
dev.off()
# identify significantly DE genes by FDR
tagsTblANOVATreatment.glm_keep <- tagsTblANOVATreatment$FDR < 0.05
# create filtered results table of DE genes
tagsTblANOVATreatment.filtered <- tagsTblANOVATreatment[tagsTblANOVATreatment.glm_keep,]


#Test whether the average across all tolerant groups is equal to the average across
#all not tolerant groups, to examine the overall effect of tolerance
con.tolerance <- makeContrasts(tolerance = (UV.Tol + VIS.Tol)/4
  - (UV.NTol + VIS.NTol)/4,
  levels=design)

#Look at genes with significant expression across all UV groups
treat.anov.tolerance <- glmTreat(fit, contrast=con.tolerance, lfc=log2(1.2))
summary(decideTests(treat.anov.tolerance))
#Write plot to file
jpeg("glmQLF_2WayANOVA_tolerance_plotMD_LFC1.2.jpg")
plotMD(treat.anov.tolerance)
abline(h=c(-1, 1), col=ghibli_colors[3])
dev.off()
#Write tags table of DE genes to file
tagsTblANOVATolerance <- topTags(treat.anov.tolerance, n=nrow(treat.anov.tolerance$table), adjust.method="fdr")$table
write.table(tagsTblANOVATolerance, file="glmQLF_2WayANOVA_tolerance_topTags_LFC1.2.csv", sep=",", row.names=TRUE, quote=FALSE)

# add column for identifying direction of DE gene expression
tagsTblANOVATolerance$topDE <- "NA"
# identify significantly up DE genes
tagsTblANOVATolerance$topDE[tagsTblANOVATolerance$logFC > 1 & tagsTblANOVATolerance$FDR < 0.05] <- "UP"
# identify significantly down DE genes
tagsTblANOVATolerance$topDE[tagsTblANOVATolerance$logFC < -1 & tagsTblANOVATolerance$FDR < 0.05] <- "DOWN"
# create volcano plot
jpeg("glmQLF_2WayANOVA_tolerance_volcano_LFC1.2.jpg")
ggplot(data=tagsTblANOVATolerance, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset, breaks = c("Up", "Down"))
dev.off()
# create volcano plot with labels
labelSetTolerance <- tagsTblANOVATolerance[tagsTblANOVATolerance$topDE == "UP" | tagsTblANOVATolerance$topDE == "DOWN",]
jpeg("glmQLF_2WayANOVA_tolerance_volcanoLabeled_LFC1.2.jpg")
ggplot(data=tagsTblANOVATolerance, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  ggrepel::geom_text_repel(data = labelSetTolerance, aes(label = row.names(labelSetTolerance)), max.overlaps=20) +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset, breaks = c("Up", "Down"))
dev.off()
# identify significantly DE genes by FDR
tagsTblANOVATolerance.glm_keep <- tagsTblANOVATolerance$FDR < 0.05
# create filtered results table of DE genes
tagsTblANOVATolerance.filtered <- tagsTblANOVATolerance[tagsTblANOVATolerance.glm_keep,]


#Test whether there is an interaction effect
con.Inter <- makeContrasts(Inter = ((UV.NTol + UV.Tol)/2
  - (VIS.NTol + VIS.Tol)/2)
  - ((UV.Tol + VIS.Tol)/2
  - (UV.NTol + VIS.NTol)/2),
  levels=design)

#Look at genes with significant expression
treat.anov.Inter <- glmTreat(fit, contrast=con.Inter, lfc=log2(1.2))
summary(decideTests(treat.anov.Inter))
#Write plot to file
jpeg("glmQLF_2WayANOVA_interaction_plotMD_LFC1.2.jpg")
plotMD(treat.anov.Inter)
abline(h=c(-1, 1), col=ghibli_colors[3])
dev.off()
#Generate table of DE genes
tagsTblANOVAInter <- topTags(treat.anov.Inter, n=nrow(treat.anov.Inter$table), adjust.method="fdr")$table
write.table(tagsTblANOVAInter, file="glmQLF_2WayANOVA_interaction_topTags_LFC1.2.csv", sep=",", row.names=TRUE, quote=FALSE)

# add column for identifying direction of DE gene expression
tagsTblANOVAInter$topDE <- "NA"
# identify significantly up DE genes
tagsTblANOVAInter$topDE[tagsTblANOVAInter$logFC > 1 & tagsTblANOVAInter$FDR < 0.05] <- "UP"
# identify significantly down DE genes
tagsTblANOVAInter$topDE[tagsTblANOVAInter$logFC < -1 & tagsTblANOVAInter$FDR < 0.05] <- "DOWN"
# create volcano plot
jpeg("glmQLF_2WayANOVA_interaction_volcano_LFC1.2.jpg")
ggplot(data=tagsTblANOVAInter, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset, breaks = c("Up", "Down"))
dev.off()
# create volcano plot with labels
labelSetInteraction <- tagsTblANOVAInter[tagsTblANOVAInter$topDE == "UP" | tagsTblANOVAInter$topDE == "DOWN",]
jpeg("glmQLF_2WayANOVA_interaction_volcanoLabeled_LFC1.2.jpg")
ggplot(data=tagsTblANOVAInter, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  ggrepel::geom_text_repel(data = labelSetInteraction, aes(label = row.names(labelSetInteraction)), max.overlaps=100) +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset, breaks = c("Up", "Down"))
dev.off()
# identify significantly DE genes by FDR
tagsTblANOVAInter.glm_keep <- tagsTblANOVAInter$FDR < 0.05
# create filtered results table of DE genes
tagsTblANOVAInter.filtered <- tagsTblANOVAInter[tagsTblANOVAInter.glm_keep,]


# GLM Results Exploration
# retrieve set of DE gene names for hours contrast
geneSet_treatment <- rownames(tagsTblANOVATreatment.filtered)
# retrieve set of DE gene names for hours contrast
geneSet_tolerance <- rownames(tagsTblANOVATolerance.filtered)
# retrieve set of DE gene names for interaction contrast
geneSet_interaction <- rownames(tagsTblANOVAInter.filtered)
# create combined glm_list of DE gene names
glm_list_venn <- list(treatment = geneSet_treatment, 
                      tolerance = geneSet_tolerance,
                      interaction = geneSet_interaction)
# create venn diagram
jpeg("glmQLF_2WayANOVA_venn_LFC1.2.jpg")
ggVennDiagram(glm_list_venn, label_alpha=0.25, category.names = c("treatment","tolerance","interaction")) +
  scale_colour_discrete(type = ghibli_subset)
dev.off()

