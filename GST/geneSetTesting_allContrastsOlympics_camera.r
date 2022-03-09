#!/usr/bin/env Rscript
#Usage: Rscript geneSetTest_camera.r countsFile factorGroupingFile
#Usage Ex: Rscript geneSetTest_camera.r PA42_v4.1_normalizedCountsOlympics_uniprot.csv expDesign_camera_Olympics.csv
#R script to perform gene set enrichment testing using camera

#Set working directory
setwd("/Users/bamflappy/PfrenderLab/dMelUV/GSTA_PA42_v4.1")

# Load libraries
library("limma")

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
tolerance <- targets$tolerance
treatment <- targets$treatment

#The experimental design is parametrized with a one-way layout, 
# where one coefficient is assigned to each group
design <- model.matrix(~ 0 + tolerance + treatment + tolerance:treatment)
colnames(design) <- c("(Intercept)","NTol.Tol","UV.VIS","Interaction")

#Import gene sets
set1 <- read.table(file="/Users/bamflappy/PfrenderLab/PA42_v4.1/geneSet1_KEGG_BASE_EXCISION_REPAIR.txt")
set2 <- read.table(file="/Users/bamflappy/PfrenderLab/PA42_v4.1/geneSet2_KEGG_HOMOLOGOUS_RECOMBINATION.txt")
set3 <- read.table(file="/Users/bamflappy/PfrenderLab/PA42_v4.1/geneSet3_KEGG_MISMATCH_REPAIR.txt")
set4 <- read.table(file="/Users/bamflappy/PfrenderLab/PA42_v4.1/geneSet4_KEGG_NON_HOMOLOGOUS_END_JOINING.txt")
set5 <- read.table(file="/Users/bamflappy/PfrenderLab/PA42_v4.1/geneSet5_KEGG_NUCLEOTIDE_EXCISION_REPAIR.txt")
set6 <- read.table(file="/Users/bamflappy/PfrenderLab/PA42_v4.1/geneSet6_KEGG_P53_SIGNALING_PATHWAY.txt")

#Vectors containing indices of genes corresponding to gene sets
index1 <- which(countsTable$sprot %in% set1$V1)
index2 <- which(countsTable$sprot %in% set2$V1)
index3 <- which(countsTable$sprot %in% set3$V1)
index4 <- which(countsTable$sprot %in% set4$V1)
index5 <- which(countsTable$sprot %in% set5$V1)
index6 <- which(countsTable$sprot %in% set6$V1)

#Individual gene set tests to check against batch and PR results
#What is this one doing?
camera(counts, index1, design)
camera(counts, index2, design)
camera(counts, index3, design)
camera(counts, index4, design)
camera(counts, index5, design)
camera(counts, index6, design)
#Test with the tolerance contrast
camera(counts, index1, design, contrast=2)
camera(counts, index2, design, contrast=2)
camera(counts, index3, design, contrast=2)
camera(counts, index4, design, contrast=2)
camera(counts, index5, design, contrast=2)
camera(counts, index6, design, contrast=2)
#Test with the treatment contrast
camera(counts, index1, design, contrast=3)
camera(counts, index2, design, contrast=3)
camera(counts, index3, design, contrast=3)
camera(counts, index4, design, contrast=3)
camera(counts, index5, design, contrast=3)
camera(counts, index6, design, contrast=3)
#Test with the interaction contrast
camera(counts, index1, design, contrast=4)
camera(counts, index2, design, contrast=4)
camera(counts, index3, design, contrast=4)
camera(counts, index4, design, contrast=4)
camera(counts, index5, design, contrast=4)
camera(counts, index6, design, contrast=4)

#Batch gene set tests
#Inter-gene correlation will be estimated for each tested set
camera(counts, list(set1=index1,set2=index2,set3=index3,set4=index4,set5=index5,set6=index6), design, inter.gene.cor=NA)
#Test with the tolerance contrast
camera(counts, list(set1=index1,set2=index2,set3=index3,set4=index4,set5=index5,set6=index6), design, contrast=2, inter.gene.cor=NA)
#Test with the treatment contrast
camera(counts, list(set1=index1,set2=index2,set3=index3,set4=index4,set5=index5,set6=index6), design, contrast=3, inter.gene.cor=NA)
#Test with the interaction contrast
camera(counts, list(set1=index1,set2=index2,set3=index3,set4=index4,set5=index5,set6=index6), design, contrast=4, inter.gene.cor=NA)

#Specifying an inter-gene correlation of 0.01
camera(counts, list(set1=index1,set2=index2,set3=index3,set4=index4,set5=index5,set6=index6), design, inter.gene.cor=0.01)
#Test with the tolerance contrast
camera(counts, list(set1=index1,set2=index2,set3=index3,set4=index4,set5=index5,set6=index6), design, contrast=2, inter.gene.cor=0.01)
#Test with the treatment contrast
camera(counts, list(set1=index1,set2=index2,set3=index3,set4=index4,set5=index5,set6=index6), design, contrast=3, inter.gene.cor=0.01)
#Test with the interaction contrast
camera(counts, list(set1=index1,set2=index2,set3=index3,set4=index4,set5=index5,set6=index6), design, contrast=4, inter.gene.cor=0.01)

#Specifying an inter-gene correlation of 0.05
camera(counts, list(set1=index1,set2=index2,set3=index3,set4=index4,set5=index5,set6=index6), design, inter.gene.cor=0.05)
#Test with the tolerance contrast
camera(counts, list(set1=index1,set2=index2,set3=index3,set4=index4,set5=index5,set6=index6), design, contrast=2, inter.gene.cor=0.05)
#Test with the treatment contrast
camera(counts, list(set1=index1,set2=index2,set3=index3,set4=index4,set5=index5,set6=index6), design, contrast=3, inter.gene.cor=0.05)
#Test with the interaction contrast
camera(counts, list(set1=index1,set2=index2,set3=index3,set4=index4,set5=index5,set6=index6), design, contrast=4, inter.gene.cor=0.05)

#Fit the linear model using the design matrix
fit <- eBayes(lmFit(counts, design))
#Output DEA results
#Tolerance contrast results
outFit2 <- topTable(fit, coef=2, adjust.method="BH", sort.by="P")
write.table(outFit2, file="eBayes_tolerance_topTable.csv", sep=",", row.names=TRUE)
#Treatment contrast results
outFit3 <- topTable(fit, coef=3, adjust.method="BH", sort.by="P")
write.table(outFit3, file="eBayes_treatment_topTable.csv", sep=",", row.names=TRUE)
#Interaction contrast results
outFit4 <- topTable(fit, coef=4, adjust.method="BH", sort.by="P")
write.table(outFit4, file="eBayes_interaction_topTable.csv", sep=",", row.names=TRUE)

#Pre-ranked (PR) gene set test version
#Using moderated t-statistic for NTol.Tol group
# corresponding to the tolerance effect
cameraPR(fit$t[,2], list(set1=index1,set2=index2,set3=index3,set4=index4,set5=index5,set6=index6))
#Using moderated t-statistic for UV.VIS group
# corresponding to the treatment effect
cameraPR(fit$t[,3], list(set1=index1,set2=index2,set3=index3,set4=index4,set5=index5,set6=index6))
#Using moderated t-statistic for Interaction group
# corresponding to the interaction effect
cameraPR(fit$t[,4], list(set1=index1,set2=index2,set3=index3,set4=index4,set5=index5,set6=index6))
#Using moderated F-statistic
cameraPR(fit$F, list(set1=index1,set2=index2,set3=index3,set4=index4,set5=index5,set6=index6))
