#!/usr/bin/env Rscript
#Usage: Rscript statistics_edgeR.r countsFile.csv startColPos endColPos
#Usage Ex: Rscript statistics_edgeR.r ../GeneCounts_Merged/merged_counts_legacy_cleaned.csv 1 6
#R script to perform statistical analysis of gene count tables using edgeR
#Install edgeR, this should only need to be done once
#Since edgeR is already installed on the CRC this can be skipped if using the module
#bioLite("edgeR")
#Load the edgeR library

#GLM Method
#Import gene count data
countsTable <- read.csv(file=args[1], row.names="gene")[ ,args[2]:args[3]]
head(countsTable)
#Add grouping factor
group <- factor(c(rep("ctrl",3),rep("treat",3)))
#Create DGE list object
list <- DGEList(counts=countsTable,group=group)

#There is no purpose in analysing genes that are not expressed in either 
# experimental condition, so genes are first filtered on expression levels
keep <- filterByExpr(list)
table(keep)
list <- list[keep, , keep.lib.sizes=FALSE]
#Calculate normalized factors
list <- calcNormFactors(list)
#View normalization factors
list$samples

#An MDS plot shows the relative similarities of the six samples
#Distances on an MDS plot of a DGEList object correspond to leading log-fold-change
# between each pair of samples, and can be used to view batch and treatment effects
plotMDS(list, col=rep(1:2, each=3))
#Draw a heatmap of individual RNA-seq samples using moderated log-counts-per-million
logcpm <- cpm(list, log=TRUE)

#Identify GO terms and KEGG pathways that are over-represented in
# group 2 (treatment) compared to group 1 (control)
qlf <- glmQLFTest(fit, coef="treat")
go <- goana(qlf, species="Dm")
topGO(go, sort="up")
keg <- kegga(qlf, species="Dm")
topKEGG(keg, sort="up")

#Design the model
design <- model.matrix(~group)
#Estimate common dispersion, trended dispersions, and tagwise dispersions
list <- estimateDisp(list,design)
#Perform quasi-likelihood F-tests
#The fit has three parameters. The first is the baseline level of group 1,
# and the second and third are the 2 vs 1 and 3 vs 1 differences
fit <- glmQLFit(list,design)
qlf <- glmQLFTest(fit,coef=2)
topTags(qlf)