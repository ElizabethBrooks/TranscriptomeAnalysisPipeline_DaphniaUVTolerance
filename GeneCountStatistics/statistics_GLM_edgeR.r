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

#Plot the library sizes before normalization
barplot(list$samples$lib.size*1e-6, names=1:36, ylab="Library size (millions)")
#There is no purpose in analysing genes that are not expressed in either 
# experimental condition, so genes are first filtered on expression levels
keep <- filterByExpr(list)
table(keep)
list <- list[keep, , keep.lib.sizes=FALSE]
#Calculate normalized factors
list <- calcNormFactors(list)
#View normalization factors
list$samples
dim(list)
#Plot the library sizes after normalization
barplot(list$samples$lib.size*1e-6, names=1:36, ylab="Library size (millions)")

#Draw a MDS plot to show the relative similarities of the samples
# and to view batch and treatment effects
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
#Compute predictive log2-foldchanges (logFC) 
logFC <- predFC(list,design,prior.count=1,dispersion=0.05)
cor(logFC[,4:6])
#Estimate common dispersion, trended dispersions, and tagwise dispersions
list <- estimateDisp(list,design)
#Estimate the genewise dispersion estimates over all genes, allowing for a
# possible abundance trend
plotBCV(list)
#The fit has three parameters. The first is the baseline level of group 1,
# and the second and third are the 2 vs 1 and 3 vs 1 differences
fit <- glmQLFit(list,design)
#Estimate Quasi-Likelihood dispersions, then visualize
plotQLDisp(fit)

#Now conduct QL F-tests for the treatment effect and show the top genes
#By default, the test is for the last coefficient in the design matrix
qlf <- glmQLFTest(fit)
topTags(qlf)
#Look at the individual counts-per-million for the top genes.
top <- rownames(topTags(qlf))
cpm(list)[top,]
#Summarize the total number of genes significantly up-regulated or
# down-regulated at 5% FDR
summary(decideTests(qlf))
#Plot all the logFCs against average count size, highlighting the DE genes
#The blue lines indicate 2-fold up or down
plotMD(qlf)
abline(h=c(-1,1), col="blue")

#Perform quasi-likelihood F-tests to test for significant differential
# expression in each gene, and show the top genes
qlf <- glmQLFTest(fit,coef=2)
topTags(qlf)
#Calculte the false discovery rate
FDR <- p.adjust(qlf$table$PValue, method="BH")
sum(FDR < 0.05)