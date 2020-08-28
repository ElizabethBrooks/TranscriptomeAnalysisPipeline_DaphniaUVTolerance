#!/usr/bin/env Rscript
#Usage: Rscript twoWayANOVA_OlympicsGrp1_edgeR.r countsFile factorGroupingFile
#Usage Ex: Rscript twoWayANOVA_OlympicsGrp1_edgeR.r geneCounts_merged_counted_htseqTophat2_run1_fullset_cleaned.csv DMelUV_ExpDesign_Olympics_GRP1.csv
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
	UV.TsNT = UV.NotTolerant-UV.Tolerant,
	VIS.TvsNT = VIS.NotTolerant-VIS.Tolerant,
	UVvsVIS.NT = UV.NotTolerant-VIS.NotTolerant,
	UVvsVIS.T = (UV.Tolerant-UV.NotTolerant)-(VIS.Tolerant-VIS.NotTolerant),
	levels=design)