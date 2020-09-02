#!/usr/bin/env Rscript
#Usage: Rscript twoWayANOVA_OlympicsGrp1_edgeR.r countsFile factorGroupingFile
#Usage Ex: Rscript twoWayANOVA_OlympicsGrp1_edgeR.r geneCounts_merged_counted_htseqTophat2_run1_fullset_cleaned.csv expDesign_Olympics_GRP1.csv
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
Group <- factor(paste(targets$Treatment,targets$Tolerance,sep="."))
cbind(targets,Group=Group)

#Create the design matrix
#Each tolerance for each treatment is a group
design <- model.matrix(~0+Group)
colnames(design) <- levels(Group)
fit <- glmQLFit(y, design)

#Contrasts for comparisons
my.contrasts <- makeContrasts(
	UV.NTvsT = UV.NotTolerant-UV.Tolerant,
	VIS.NTvsT = VIS.NotTolerant-VIS.Tolerant,
	UVvsVIS.NT = UV.NotTolerant-VIS.NotTolerant,
	UVvsVIS.T = (UV.Tolerant-UV.NotTolerant)-(VIS.Tolerant-VIS.NotTolerant),
	levels = design)

#To find genes responding to UV at low tolerance
tested <- glmQLFTest(fit, contrast=my.contrasts[,"UV.NTvsT"])
#To find genes responding to VIS at low tolerance
tested <- glmQLFTest(fit, contrast=my.contrasts[,"VIS.NTvsT"])
#To find genes with baseline differences between UV and VIS at low tolerance
tested <- glmQLFTest(fit, contrast=my.contrasts[,"UVvsVIS.NT"])
#To find genes with baseline differences between UV and VIS at high tolerance
tested <- glmQLFTest(fit, contrast=my.contrasts[,"UVvsVIS.T"])

#View the top DE genes
topTags(tested)
#Create results table of DE genes
resultsTbl <- topTags(tested, n=nrow(tested$table))$table
#Output resulting table
write.table(resultsTbl, file="stats_exactTest.csv", sep=",", row.names=TRUE)