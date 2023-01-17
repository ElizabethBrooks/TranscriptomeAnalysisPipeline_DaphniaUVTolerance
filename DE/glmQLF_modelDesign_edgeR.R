#!/usr/bin/env Rscript
#Usage: Rscript glmLRT_edgeR.r countsFile startColumn endColumn factorGroupingFile FDR
#Usage Ex: Rscript glmLRT_edgeR.r cleaned.csv 1 24 expDesign_Olympics_GRP1.csv 0.10
#R script to perform statistical analysis of gene count tables using edgeR two way ANOVA

#Set working directory
#workingDir = args[1];
workingDir="/Users/bamflappy/PfrenderLab/DEA_PA42_v4.1/glmQLFAnalysis_modelDesign_FDR0.10"
setwd(workingDir); 

#Install edgeR and statmod, this should only need to be done once
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("edgeR")
#install.packages("statmod")

#Load the edgeR library
library("edgeR")
library("statmod")

#Retrieve input file name of gene counts
#args = commandArgs(trailingOnly=TRUE)
#Test if there is one input argument
#if (length(args)!=5) {
#  stop("Two file names and a range of columns must be supplied.n", call.=FALSE)
#}

#Import gene count data
#countsTable <- read.csv(file=args[1], row.names="gene")[ ,args[2]:args[3]]
inputTable <- read.csv(file="/Users/bamflappy/PfrenderLab/PA42_v4.1/geneCounts_cleaned_PA42_v4.1.csv", row.names="gene")[ ,1:24]

#Trim the data table
countsTable <- head(inputTable, - 5)

#Import grouping factor
#targets <- read.csv(file=args[4], row.names="sample")
#targets <- read.csv(file="/Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_binned_Olympics.csv", row.names="sample")
targets <- read.csv(file="/Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_Olympics.csv", row.names="sample")
#Retrieve input FDR cutoff
#fdrCut=as.numeric(args[5])
#fdrCut=0.10

#Setup a design matrix
group <- factor(paste(targets$treatment,targets$genotype,sep="."))
#cbind(targets,Group=group)
#Create DGE list object
list <- DGEList(counts=countsTable,group=group)
colnames(list) <- rownames(targets)

#Plot the library sizes before normalization
#jpeg("glmQLF_plotBarsBefore.jpg")
barplot(list$samples$lib.size*1e-6, names=1:24, ylab="Library size (millions)")
#dev.off()

#Retain genes only if it is expressed at a minimum level
keep <- filterByExpr(list)
summary(keep)
list <- list[keep, , keep.lib.sizes=FALSE]

#Use TMM normalization to eliminate composition biases
# between libraries
list <- calcNormFactors(list)
#list$samples
#Write normalized counts to file
#normList <- cpm(list, normalized.lib.sizes=TRUE)
#write.table(normList, file="glmQLF_modelDesign_normalizedCounts.csv", sep=",", row.names=TRUE)

#Verify TMM normalization using a MD plot
#Write plot to file
#jpeg("glmQLF_plotMDBefore.jpg")
plotMD(cpm(list, log=TRUE), column=1)
abline(h=0, col="red", lty=2, lwd=2)
#dev.off()

#Use a MDS plot to visualizes the differences
# between the expression profiles of different samples
points <- c(0,1,2,3,15,16,17,18)
colors <- rep(c("blue", "darkgreen", "red", "black"), 2)
#Write plot with legend to file
#jpeg("glmLRT_plotMDS.jpg")
plotMDS(list, col=colors[group], pch=points[group])
legend("topleft", legend=levels(group), pch=points, col=colors, ncol=2)
#dev.off()
#Write plot without legend to file
#jpeg("glmQLF_plotMDS_noLegend.jpg")
plotMDS(list, col=colors[group], pch=points[group])
#dev.off()


#Define each treatment combination as a group
#The experimental design is parametrized with a one-way layout, 
# where one coefficient is assigned to each group
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
#design


#Nested interaction formula
#Make sure that VIS is the reference level
#targets$treatment <- relevel(factor(targets$treatment), ref="VIS")
#Form the design matrix to consider all the levels of genotype for each treatment
#design <- model.matrix(~treatment + treatment:genotype, data=targets)


#Factorial model
#Interaction at any time
##design <- model.matrix(~treatment * genotype, data=targets)
#design <- model.matrix(~treatment + genotype + treatment:genotype, data=targets)


#Next, the NB dispersion is estimated
list <- estimateDisp(list, design, robust=TRUE)
#list$common.dispersion
#Visualize the dispersion estimates with a BCV plot
#Write plot to file
#jpeg("glmQLF_plotBCV.jpg")
plotBCV(list)
#dev.off()

#Fit the design model
fit <- glmQLFit(list, design)
colnames(fit)


#Make the desired contrasts
my.contrasts <- makeContrasts(
  treatment.E05 = UV.E05-VIS.E05,
  treatment.R2 = UV.R2-VIS.R2,
  treatment.Y05 = UV.Y05-VIS.Y05,
  treatment.Y023 = UV.Y023-VIS.Y023,
  treatment.Y023 = UV.Y023-VIS.Y023,
  levels=design)

#Test E05 contrast baseline differences between the UV and the VIS
qlf.E05 <- glmQLFTest(fit, contrast=my.contrasts[,"treatment.E05"])
summary(decideTests(qlf.E05))
plotMD(qlf.E05)
#Test E05 contrast baseline differences between the UV and the VIS
qlf.R2 <- glmQLFTest(fit, contrast=my.contrasts[,"treatment.R2"])
summary(decideTests(qlf.R2))
plotMD(qlf.R2)
#Test E05 contrast baseline differences between the UV and the VIS
qlf.Y05 <- glmQLFTest(fit, contrast=my.contrasts[,"treatment.Y05"])
summary(decideTests(qlf.Y05))
plotMD(qlf.Y05)
#Test E05 contrast baseline differences between the UV and the VIS
qlf.Y023 <- glmQLFTest(fit, contrast=my.contrasts[,"treatment.Y023"])
summary(decideTests(qlf.Y023))
plotMD(qlf.Y023)

