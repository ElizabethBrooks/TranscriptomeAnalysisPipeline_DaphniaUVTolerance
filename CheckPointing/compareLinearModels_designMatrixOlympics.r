#Install edgeR and statmod, this should only need to be done once
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("edgeR")
#install.packages("statmod")

#Turn off scientific notation
#options(scipen = 999)

#Load the edgeR library
library("edgeR")
library("statmod")
library("limma")

#Import gene count data
countsTable <- read.csv(file=args[1], row.names="gene")[ ,args[2]:args[3]]
#head(countsTable)
#Import grouping factor
targets <- read.csv(file=args[4], row.names="sample")
#Retrieve input FDR cutoff
#fdrCut=as.numeric(args[5])

#Setup a design matrix
group <- factor(paste(targets$treatment,targets$genotype,sep="."))
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
write.table(normList, file="glmQLF_normalizedCounts.csv", sep=",", row.names=TRUE)

#Verify TMM normalization using a MD plot
#Write plot to file
jpeg("glmQLF_plotMDBefore.jpg")
plotMD(cpm(list, log=TRUE), column=1)
abline(h=0, col="red", lty=2, lwd=2)
dev.off()

#Use a MDS plot to visualizes the differences
# between the expression profiles of different samples
points <- c(0,1,2,3,15,16,17,18)
colors <- rep(c("blue", "darkgreen", "red", "black"), 2)
#Write plot with legend to file
jpeg("glmQLF_plotMDS.jpg")
plotMDS(list, col=colors[group], pch=points[group])
legend("topleft", legend=levels(group), pch=points, col=colors, ncol=2)
dev.off()
#Write plot without legend to file
jpeg("glmQLF_plotMDS_noLegend.jpg")
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
con.UVvsVIS <- makeContrasts(UVvsVIS = (UV.E05 + UV.R2 + UV.Y023 + UV.Y05)/4
                             - (VIS.E05 + VIS.R2 + VIS.Y023 + VIS.Y05)/4,
                             levels=design)

#Look at genes with significant expression across all UV groups
treat.anov.UVVIS <- glmTreat(fit, contrast=con.UVvsVIS, lfc=log2(1.2))
summary(decideTests(treat.anov.UVVIS))
#Write plot to file
jpeg("glmQLF_2WayANOVA_UVvsVIS_plotMD_LFC1.2.jpg")
plotMD(treat.anov.UVVIS)
abline(h=c(-1, 1), col="blue")
dev.off()
#Write tags table of DE genes to file
tagsTblANOVA.filtered <- topTags(treat.anov.UVVIS, n=nrow(treat.anov.UVVIS$table), adjust.method="fdr")$table
#tagsTblANOVA.filtered.keep <- tagsTblANOVA.filtered$FDR <= fdrCut
#tagsTblANOVA.filtered.out <- tagsTblANOVA.filtered[tagsTblANOVA.filtered.keep,]
#write.table(tagsTblANOVA.filtered, file="glmQLF_2WayANOVA_UVvsVIS_topTags_LFC1.2.csv", sep=",", row.names=TRUE)


#Test whether the average across all tolerant groups is equal to the average across
#all not tolerant groups, to examine the overall effect of tolerance
con.TvsN <- makeContrasts(TvsN = (UV.Y05 + VIS.Y05 + UV.Y023 + VIS.Y023)/4
                          - (UV.E05 + VIS.E05 + UV.R2 + VIS.R2)/4,
                          levels=design)

#Look at genes with significant expression across all UV groups
treat.anov.TN <- glmTreat(fit, contrast=con.TvsN, lfc=log2(1.2))
summary(decideTests(treat.anov.TN))
#Write plot to file
jpeg("glmQLF_2WayANOVA_TvsN_plotMD_LFC1.2.jpg")
plotMD(treat.anov.TN)
abline(h=c(-1, 1), col="blue")
dev.off()
#Write tags table of DE genes to file
tagsTblANOVATN.filtered <- topTags(treat.anov.TN, n=nrow(treat.anov.TN$table), adjust.method="fdr")$table
#tagsTblANOVATN.filtered.keep <- tagsTblANOVATN.filtered$FDR <= fdrCut
#tagsTblANOVATN.filtered.out <- tagsTblANOVATN.filtered[tagsTblANOVATN.filtered.keep,]
#write.table(tagsTblANOVATN.filtered, file="glmQLF_2WayANOVA_TvsN_topTags_LFC1.2.csv", sep=",", row.names=TRUE)


#Test whether there is an interaction effect
con.Inter <- makeContrasts(Inter = ((UV.E05 + UV.R2 + UV.Y023 + UV.Y05)/4
                                    - (VIS.E05 + VIS.R2 + VIS.Y023 + VIS.Y05)/4)
                           - ((UV.Y05 + VIS.Y05 + UV.Y023 + VIS.Y023)/4
                              - (UV.E05 + VIS.E05 + UV.R2 + VIS.R2)/4),
                           levels=design)

#Look at genes with significant expression
treat.anov.Inter <- glmTreat(fit, contrast=con.Inter, lfc=log2(1.2))
summary(decideTests(treat.anov.Inter))
#Write plot to file
jpeg("glmQLF_2WayANOVA_interaction_plotMD_LFC1.2.jpg")
plotMD(treat.anov.Inter)
abline(h=c(-1, 1), col="blue")
dev.off()
#Generate table of DE genes
tagsTblANOVAInter.filtered <- topTags(treat.anov.Inter, n=nrow(treat.anov.Inter$table), adjust.method="fdr")$table
#tagsTblANOVAInter.filtered.keep <- tagsTblANOVAInter.filtered$FDR <= fdrCut
#tagsTblANOVAInter.filtered.out <- tagsTblANOVAInter.filtered[tagsTblANOVAInter.filtered.keep,]
#write.table(tagsTblANOVAInter.filtered, file="glmQLF_2WayANOVA_interaction_topTags_LFC1.2.csv", sep=",", row.names=TRUE)

#Import gene count data
#Row names are set to the uniprot IDs
# to match the gene IDs of the MSigDB KEGG DNA repair gene sets
#countsTable <- read.csv(file=args[1])
countsTable <- read.csv(file="/Users/bamflappy/PfrenderLab/PA42_v4.1/PA42_v4.1_normalizedLogCountsOlympics_uniprot.csv")

#Create a subset of the input counts table containing only the gene counts
counts <- countsTable[3:26]

#Import grouping factor
#targets <- read.csv(file=args[2], row.names="sample")
targets <- read.csv(file="/Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_camera_Olympics.csv", row.names="sample")

#Setup a design matrix
tolerance <- targets$treatment
treatment <- targets$tolerance

#The experimental design is parametrized with a one-way layout, 
# where one coefficient is assigned to each group
design <- model.matrix(~ 0 + tolerance + treatment + tolerance:treatment)
colnames(design) <- c("(Intercept)","NTol.Tol","UV.VIS","Interaction")

#Fit the linear model using the design matrix
fit <- eBayes(lmFit(counts, design))
#Output DEA results
#Tolerance contrast results
outFit2 <- topTable(fit, coef=2, adjust.method="BH", sort.by="P")
#write.table(outFit2, file="eBayes_tolerance_topTable.csv", sep=",", row.names=TRUE)
#Treatment contrast results
outFit3 <- topTable(fit, coef=3, adjust.method="BH", sort.by="P")
#write.table(outFit3, file="eBayes_treatment_topTable.csv", sep=",", row.names=TRUE)
#Interaction contrast results
outFit4 <- topTable(fit, coef=4, adjust.method="BH", sort.by="P")
#write.table(outFit4, file="eBayes_interaction_topTable.csv", sep=",", row.names=TRUE)
