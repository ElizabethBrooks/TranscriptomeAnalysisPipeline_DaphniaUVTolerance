#!/usr/bin/env Rscript
#Usage: Rscript glmLRT_edgeR.r countsFile startColumn endColumn factorGroupingFile
#Usage Ex: Rscript glmLRT_edgeR.r cleaned.csv 1 24 expDesign_Olympics_GRP1.csv
#R script to perform statistical analysis of gene count tables using edgeR two way ANOVA

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
#Test if there is one input argument
if (length(args)!=4) {
  stop("Two file names and a range of columns must be supplied.n", call.=FALSE)
}

#Import gene count data
countsTable <- read.csv(file=args[1], row.names="gene")[ ,args[2]:args[3]]
#head(countsTable)
#Import grouping factor
targets <- read.csv(file=args[4], row.names="sample")

#Setup a design matrix
group <- factor(paste(targets$treatment,targets$genotype,sep="."))
#cbind(targets,Group=group)
#Create DGE list object
list <- DGEList(counts=countsTable,group=group)
colnames(list) <- targets$sample

#Plot the library sizes before normalization
jpeg("glmLRT_plotBarsBefore.jpg")
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
write.table(normList, file="glmLRT_normalizedCounts.csv", sep=",", row.names=TRUE)

#Verify TMM normalization using a MD plot
#Write plot to file
jpeg("glmLRT_plotMDBefore.jpg")
plotMD(cpm(list, log=TRUE), column=1)
abline(h=0, col="red", lty=2, lwd=2)
dev.off()

#Use a MDS plot to visualizes the differences
# between the expression profiles of different samples
points <- c(0,1,2,3,15,16,17,18)
colors <- rep(c("blue", "darkgreen", "red", "black"), 2)
#Write plot with legend to file
jpeg("glmLRT_plotMDS.jpg")
plotMDS(list, col=colors[group], pch=points[group])
legend("topleft", legend=levels(group), pch=points, col=colors, ncol=2)
dev.off()
#Write plot without legend to file
jpeg("glmLRT_plotMDS_noLegend.jpg")
plotMDS(list, col=colors[group], pch=points[group])
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
jpeg("glmLRT_plotBCV.jpg")
plotBCV(list)
dev.off()

#Now, estimate the dispersions
fit <- glmFit(list, design, robust=TRUE)
#head(fit$coefficients)

#Define a matrix of contrasts, where each column
# represents a contrast between two groups of interest
con.allPairs <- makeContrasts(
  #Pairwise
  E05.UVvsVIS = UV.E05 - VIS.E05,
  R2.UVvsVIS = UV.R2 - VIS.R2,
  Y023.UVvsVIS = UV.Y023 - VIS.Y023,
  Y05.UVvsVIS = UV.Y05 - VIS.Y05,
  levels=design)

#All pairs using likelihood ratio test
test.allPairs <- glmLRT(fit, contrast=con.allPairs)
summary(decideTests(test.allPairs))
#Write plot to file
jpeg("glmLRT_allPairwise_plotMD.jpg")
plotMD(test.allPairs)
abline(h=c(-1, 1), col="blue")
dev.off()
#Write tags table of DE genes to file
tagsTblAllPairwise <- topTags(test.allPairs, n=nrow(test.allPairs$table))$table
tagsTblAllPairwise.keep <- tagsTblAllPairwise$FDR <= 0.05
tagsTblAllPairwise.out <- tagsTblAllPairwise[tagsTblAllPairwise.keep,]
write.table(tagsTblAllPairwise.out, file="glmLRT_allPairwise_topTags.csv", sep=",", row.names=TRUE)

#Pairwise E05.UVvsVIS
#Test whether the differential expression is significant
treat.E05.UVvsVIS <- glmTreat(fit, contrast=con.allPairs[,"E05.UVvsVIS"], lfc=log2(1.2))
summary(decideTests(treat.E05.UVvsVIS))
#Write plot to file
jpeg("glmLRT_E05Pairwise_plotMD_filtered.jpg")
plotMD(treat.E05.UVvsVIS)
abline(h=c(-1, 1), col="blue")
dev.off()
#Write tags table of DE genes to file
tagsTblE05Pairwise <- topTags(treat.E05.UVvsVIS, n=nrow(treat.E05.UVvsVIS$table))$table
tagsTblE05Pairwise.keep <- tagsTblE05Pairwise$FDR <= 0.05
tagsTblE05Pairwise.out <- tagsTblE05Pairwise[tagsTblE05Pairwise.keep,]
write.table(tagsTblE05Pairwise.out, file="glmLRT_E05Pairwise_topTags_filtered.csv", sep=",", row.names=TRUE)

#Pairwise R2.UVvsVIS
#Test whether the differential expression is significant
treat.R2.UVvsVIS <- glmTreat(fit, contrast=con.allPairs[,"R2.UVvsVIS"], lfc=log2(1.2))
summary(decideTests(treat.R2.UVvsVIS))
#Write plot to file
jpeg("glmLRT_R2Pairwise_plotMD_filtered.jpg")
plotMD(treat.R2.UVvsVIS)
abline(h=c(-1, 1), col="blue")
dev.off()
#Write tags table of DE genes to file
tagsTblR2Pairwise <- topTags(treat.R2.UVvsVIS, n=nrow(treat.R2.UVvsVIS$table))$table
tagsTblR2Pairwise.keep <- tagsTblR2Pairwise$FDR <= 0.05
tagsTblR2Pairwise.out <- tagsTblR2Pairwise[tagsTblR2Pairwise.keep,]
write.table(tagsTblR2Pairwise.out, file="glmLRT_R2Pairwise_topTags_filtered.csv", sep=",", row.names=TRUE)

#Pairwise Y023.UVvsVIS
#Test whether the differential expression is significant
treat.Y023.UVvsVIS <- glmTreat(fit, contrast=con.allPairs[,"Y023.UVvsVIS"], lfc=log2(1.2))
summary(decideTests(treat.Y023.UVvsVIS))
#Write plot to file
jpeg("glmLRT_Y023Pairwise_plotMD_filtered.jpg")
plotMD(treat.Y023.UVvsVIS)
abline(h=c(-1, 1), col="blue")
dev.off()
#Write tags table of DE genes to file
tagsTblY023Pairwise <- topTags(treat.Y023.UVvsVIS, n=nrow(treat.Y023.UVvsVIS$table))$table
tagsTblY023Pairwise.keep <- tagsTblY023Pairwise$FDR <= 0.05
tagsTblY023Pairwise.out <- tagsTblY023Pairwise[tagsTblY023Pairwise.keep,]
write.table(tagsTblY023Pairwise.out, file="glmLRT_Y023Pairwise_topTags_filtered.csv", sep=",", row.names=TRUE)

#Pairwise Y05.UVvsVIS
#Test whether the differential expression is significant
treat.Y05.UVvsVIS <- glmTreat(fit, contrast=con.allPairs[,"Y05.UVvsVIS"], lfc=log2(1.2))
summary(decideTests(treat.Y05.UVvsVIS))
#Write plot to file
jpeg("glmLRT_Y05Pairwise_plotMD_filtered.jpg")
plotMD(treat.Y05.UVvsVIS)
abline(h=c(-1, 1), col="blue")
dev.off()
#Write tags table of DE genes to file
tagsTblY05Pairwise <- topTags(treat.Y05.UVvsVIS, n=nrow(treat.Y05.UVvsVIS$table))$table
tagsTblY05Pairwise.keep <- tagsTblY05Pairwise$FDR <= 0.05
tagsTblY05Pairwise.out <- tagsTblY05Pairwise[tagsTblY05Pairwise.keep,]
write.table(tagsTblY05Pairwise.out, file="glmLRT_Y05Pairwise_topTags_filtered.csv", sep=",", row.names=TRUE)

#ANOVA like comparisons of UV using likelihood ratio test
anov.UV <- makeContrasts(UV.R2 - UV.E05,
  UV.Y023 - UV.E05,
  UV.Y05 - UV.E05,
  levels=design)
#Look at genes using likelihood ratio test
test.anov.UV <- glmLRT(fit, contrast=anov.UV)
summary(decideTests(test.anov.UV))
#Write plot to file
jpeg("glmLRT_UV1WayANOVA_plotMD.jpg")
plotMD(test.anov.UV)
abline(h=c(-1, 1), col="blue")
dev.off()
#Write tags table of DE genes to file
tagsTblUVANOVA <- topTags(test.anov.UV, n=nrow(test.anov.UV$table))$table
tagsTblUVANOVA.keep <- tagsTblUVANOVA$FDR <= 0.05
tagsTblUVANOVA.out <- tagsTblUVANOVA[tagsTblUVANOVA.keep,]
write.table(tagsTblUVANOVA.out, file="glmLRT_UV1WayANOVA_topTags.csv", sep=",", row.names=TRUE)

#ANOVA like comparisons of VIS
anov.VIS <- makeContrasts(VIS.R2 - VIS.E05,
  VIS.Y023 - VIS.E05,
  VIS.Y05 - VIS.E05,
  levels=design)
#Look at genes using likelihood ratio test
test.anov.VIS <- glmLRT(fit, contrast=anov.VIS)
summary(decideTests(test.anov.VIS))
#Write plot to file
jpeg("glmLRT_VIS1WayANOVA_plotMD.jpg")
plotMD(test.anov.VIS)
abline(h=c(-1, 1), col="blue")
dev.off()
#Write tags table of DE genes to file
tagsTblVISANOVA <- topTags(test.anov.VIS, n=nrow(test.anov.VIS$table))$table
tagsTblVISANOVA.keep <- tagsTblVISANOVA$FDR <= 0.05
tagsTblVISANOVA.out <- tagsTblVISANOVA[tagsTblVISANOVA.keep,]
write.table(tagsTblVISANOVA.out, file="glmLRT_VIS1WayANOVA_topTags.csv", sep=",", row.names=TRUE)

#Test whether the average across all UV groups is equal to the average across
#all VIS groups, to examine the overall effect of treatment
con.UVvsVIS <- makeContrasts(UVvsVIS = (UV.E05 + UV.R2 + UV.Y023 + UV.Y05)/4
  - (VIS.E05 + VIS.R2 + VIS.Y023 + VIS.Y05)/4,
  levels=design)

#Look at genes expressed across all UV groups using likelihood ratio test
test.anov.UVVIS <- glmLRT(fit, contrast=con.UVvsVIS)
summary(decideTests(test.anov.UVVIS))
#Write plot to file
jpeg("glmLRT_2WayANOVA_plotMD.jpg")
plotMD(test.anov.UVVIS)
abline(h=c(-1, 1), col="blue")
dev.off()
#Write tags table of DE genes to file
tagsTblANOVA <- topTags(test.anov.UVVIS, n=nrow(test.anov.UVVIS$table))$table
tagsTblANOVA.keep <- tagsTblANOVA$FDR <= 0.05
tagsTblANOVA.out <- tagsTblANOVA[tagsTblANOVA.keep,]
write.table(tagsTblANOVA.out, file="glmLRT_2WayANOVA_topTags.csv", sep=",", row.names=TRUE)

#Look at genes with significant expression across all UV groups
treat.anov.UVVIS <- glmTreat(fit, contrast=con.UVvsVIS, lfc=log2(1.2))
summary(decideTests(treat.anov.UVVIS))
#Write plot to file
jpeg("glmLRT_2WayANOVA_plotMD_filtered.jpg")
plotMD(treat.anov.UVVIS)
abline(h=c(-1, 1), col="blue")
dev.off()
#Write tags table of DE genes to file
tagsTblANOVA.filtered <- topTags(treat.anov.UVVIS, n=nrow(treat.anov.UVVIS$table))$table
tagsTblANOVA.filtered.keep <- tagsTblANOVA.filtered$FDR <= 0.05
tagsTblANOVA.filtered.out <- tagsTblANOVA.filtered[tagsTblANOVA.filtered.keep,]
write.table(tagsTblANOVA.filtered.out, file="glmLRT_2WayANOVA_topTags_filtered.csv", sep=",", row.names=TRUE)