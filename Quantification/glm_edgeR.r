#!/usr/bin/env Rscript
#Usage: Rscript glm_edgeR.r countsFile factorGroupingFile
#Usage Ex: Rscript glm_edgeR.r geneCounts_merged_counted_htseqTophat2_run1_fullset_cleaned.csv expDesign_Olympics_GRP1.csv
#R script to perform statistical analysis of gene count tables using edgeR two way ANOVA
#Install edgeR, this should only need to be done once
#Since edgeR is already installed on the CRC this can be skipped if using the module
#bioLite("edgeR")

#Load the edgeR library
library("edgeR")
#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)
#Test if there is one input argument
if (length(args)!=3) {
  stop("One file name and a range of columns must be supplied.n", call.=FALSE)
}

#Import gene count data
countsTable <- read.csv(file=args[1], row.names="gene")
head(countsTable)
#Import grouping factor
targets <- read.csv(file=args[2], row.names="sample")
#Setup a design matrix
group <- factor(paste(targets$Treatment,targets$Tolerance,sep="."))
#cbind(targets,Group=group)
#Create DGE list object
list <- DGEList(counts=countsTable,group=group)

#There is no purpose in analysing genes that are not expressed in either 
# experimental condition, so genes are first filtered on expression levels
keep <- filterByExpr(list)
table(keep)
list <- list[keep, , keep.lib.sizes=FALSE]
#Calculate normalized factors
list <- calcNormFactors(list)
#Write normalized counts to file
normList <- cpm(list, normalized.lib.sizes=TRUE)
write.table(normList, file="stats_normalizedCounts.csv", sep=",", row.names=TRUE)
#View normalization factors
list$samples
dim(list)

#Produce a matrix of pseudo-counts
#Estimate common dispersion and tagwise dispersions
list <- estimateDisp(list)
list$common.dispersion
#View dispersion estimates and biological coefficient of variation
jpeg("plotBCV.jpg")
plotBCV(list)
dev.off()

#Groups
#UV.TolerantY05 - C.G1
#UV.TolerantE05 - C.G2
#UV.NotTolerantR2 - C.G3
#UV.NotTolerantY023 - C.G4
#VIS.TolerantY05 - T.G1
#VIS.TolerantE05 - T.G2
#VIS.NotTolerantR2 - T.G3
#VIS.NotTolerantY023 - T.G4

#Create the design matrix
#Each tolerance for each treatment is a group
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
tested <- glmQLFit(list, design)

#View top DE genes
topTags(tested)
#Create results table of DE genes
resultsTbl <- topTags(tested, n=nrow(tested$table))$table
#Output resulting table
write.table(resultsTbl, file="stats_glm.csv", sep=",", row.names=TRUE)