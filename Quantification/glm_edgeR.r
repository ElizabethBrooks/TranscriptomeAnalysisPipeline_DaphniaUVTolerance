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
group <- factor(paste(targets$Treatment,targets$Genotype,sep="."))
#cbind(targets,Group=group)
#Create DGE list object
list <- DGEList(counts=countsTable,group=group)
colnames(list) <- targets$Sample

#Retain genes only if it is expressed at a minimum level
keep <- filterByExpr(list)
summary(keep)
list <- list[keep, , keep.lib.sizes=FALSE]

#Use TMM normalization to eliminate composition biases
# between libraries
list <- calcNormFactors(list)
list$samples

#Verify TMM normalization using a MD plot
plotMD(cpm(list, log=TRUE), column=1)
abline(h=0, col="red", lty=2, lwd=2)

#Use a MDS plot to visualizes the differences
# between the expression profiles of different samples
points <- c(0,1,2,3,15,16,17,18)
colors <- rep(c("blue", "darkgreen", "red", "black"), 2)
plotMDS(y, col=colors[group], pch=points[group])
legend("topleft", legend=levels(group), pch=points, col=colors, ncol=2)

#The experimental design is parametrized with a one-way layout, 
# where one coefficient is assigned to each group
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
design
attr(,"assign")
attr(,"contrasts")
attr(,"contrasts")$group

#Next, the NB dispersion is estimated
list <- estimateDisp(y, design, robust=TRUE)
list$common.dispersion
#Visualize the dispersion estimates with a BCV plot
plotBCV(list)

#Now, estimate and plot the QL dispersions
fit <- glmQLFit(list, design, robust=TRUE)
head(fit$coefficients)
plotQLDisp(fit)

#Grouping factors
#UV.E05
#UV.R2
#UV.Y023
#UV.Y05
#VIS.E05
#VIS.R2
#VIS.Y023
#VIS.Y05

#Define a matrix of contrasts, where each column
# represents a contrast between two groups of interest
con <- makeContrasts(
	#Control
	UV.E05vsR2 = UV.E05 - UV.R2, #A-B
	UV.E05vsY05 = UV.E05 - UV.Y05, #A-D
	UV.R2vsY023 = UV.R2 - UV.Y023, #B-C
	UV.Y023vsY05 = UV.Y023 - UV.Y05, #C-D
	#Treatment
	VIS.E05vsR2 = VIS.E05 - VIS.R2, #A-B
	VIS.E05vsY05 = VIS.E05 - VIS.Y05, #A-D
	VIS.R2vsY023 = VIS.R2 - VIS.Y023, #B-C
	VIS.Y023vsY05 = VIS.Y023 - VIS.Y05, #C-D
	#Pairwise
	E05.UVvsVIS = UV.E05 - VIS.E05,
	R2.UVvsVIS = UV.R2 - VIS.R2,
	Y023.UVvsVIS = UV.Y023 - VIS.Y023,
	Y05.UVvsVIS = UV.Y05 - VIS.Y05,
	
	levels=design)