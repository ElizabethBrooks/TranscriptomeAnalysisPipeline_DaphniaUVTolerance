#!/usr/bin/env Rscript
#Usage: Rscript glmQLF_edgeR_pvalues.r countsFile startColumn endColumn factorGroupingFile FDR
#Usage Ex: Rscript glmQLF_edgeR_pvalues.r cleaned.csv 1 24 expDesign_Olympics_GRP1.csv 0.10
#R script to perform statistical analysis of gene count tables using edgeR GLM

#Install edgeR and statmod, this should only need to be done once
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("edgeR")
#install.packages("statmod")

#Load the edgeR library
library("edgeR")
library("statmod")

#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)

#Import gene count data
countsTable <- read.csv(file=args[1], row.names="gene")[ ,args[2]:args[3]]
#head(countsTable)
#Import grouping factor
targets <- read.csv(file=args[4], row.names="sample")
#Retrieve input FDR cutoff
fdrCut=as.numeric(args[5])

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

#Use a MDS plot to visualizes the differences
# between the expression profiles of different samples
points <- c(0,1,2,3,15,16,17,18)
colors <- rep(c("blue", "darkgreen", "red", "black"), 2)

#The experimental design is parametrized with a one-way layout, 
# where one coefficient is assigned to each group
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

#Next, the NB dispersion is estimated
list <- estimateDisp(list, design, robust=TRUE)

#Now, estimate and plot the QL dispersions
fit <- glmQLFit(list, design, robust=TRUE)

#Test whether the average across all UV groups is equal to the average across
#all VIS groups, to examine the overall effect of treatment
con.UVvsVIS <- makeContrasts(UVvsVIS = (UV.E05 + UV.R2 + UV.Y023 + UV.Y05)/4
  - (VIS.E05 + VIS.R2 + VIS.Y023 + VIS.Y05)/4,
  levels=design)

#Look at genes expressed across all UV groups using QL F-test
test.anov.UVVIS <- glmQLFTest(fit, contrast=con.UVvsVIS)
#Write tags table of DE genes to file
tagsTblANOVA <- topTags(test.anov.UVVIS, n=nrow(test.anov.UVVIS$table))$table
tagsTblANOVA.keep <- tagsTblANOVA$FDR <= fdrCut
tagsTblANOVA.out <- tagsTblANOVA[tagsTblANOVA.keep,]
write.table(tagsTblANOVA.out, file="glmQLF_2WayANOVA_UVvsVIS_topTags.csv", sep=",", row.names=TRUE)

#Test whether the average across all tolerant groups is equal to the average across
#all not tolerant groups, to examine the overall effect of tolerance
con.TvsN <- makeContrasts(TvsN = (UV.Y05 + VIS.Y05 + UV.E05 + VIS.E05)/4
  - (UV.Y023 + VIS.Y023 + UV.R2 + VIS.R2)/4,
  levels=design)

#Look at genes expressed across all UV groups using QL F-test
test.anov.TN <- glmQLFTest(fit, contrast=con.TvsN)
#Write tags table of DE genes to file
tagsTblANOVATN <- topTags(test.anov.TN, n=nrow(test.anov.TN$table))$table
tagsTblANOVATN.keep <- tagsTblANOVATN$FDR <= fdrCut
tagsTblANOVATN.out <- tagsTblANOVATN[tagsTblANOVATN.keep,]
write.table(tagsTblANOVATN.out, file="glmQLF_2WayANOVA_TvsN_topTags.csv", sep=",", row.names=TRUE)

#Test whether there is an interaction effect
con.Inter <- makeContrasts(Inter = ((UV.E05 + UV.R2 + UV.Y023 + UV.Y05)/4
  - (VIS.E05 + VIS.R2 + VIS.Y023 + VIS.Y05)/4)
  - ((UV.Y05 + VIS.Y05 + UV.E05 + VIS.E05)/4
  - (UV.Y023 + VIS.Y023 + UV.R2 + VIS.R2)/4),
  levels=design)

#Look at genes expressed across all UV groups using QL F-test
test.anov.Inter <- glmQLFTest(fit, contrast=con.Inter)
#Write tags table of DE genes to file
tagsTblANOVAInter <- topTags(test.anov.Inter, n=nrow(test.anov.Inter$table))$table
tagsTblANOVAInter.keep <- tagsTblANOVAInter$FDR <= fdrCut
tagsTblANOVAInter.out <- tagsTblANOVAInter[tagsTblANOVAInter.keep,]
write.table(tagsTblANOVAInter.out, file="glmQLF_2WayANOVA_interaction_topTags.csv", sep=",", row.names=TRUE)