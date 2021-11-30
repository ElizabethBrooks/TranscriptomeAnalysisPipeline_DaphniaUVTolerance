#!/usr/bin/env Rscript
#Usage: Rscript geneSetTest_camera.r countsFile factorGroupingFile
#Usage Ex: Rscript geneSetTest_camera.r PA42_v4.1_normalizedCountsOlympics_uniprot.csv expDesign_binned_Olympics.csv
#R script to perform gene set enrichment testing using camera

# Load libraries
library("limma")

#Import gene count data
#Row names are set to the uniprot IDs
# to match the gene IDs of the MSigDB KEGG DNA repair gene sets
#countsTable <- read.csv(file=args[1], row.names="sprot")
countsTableIn <- read.csv(file="/Users/bamflappy/PfrenderLab/PA42_v4.1/PA42_v4.1_normalizedCountsOlympics_uniprot.csv", row.names="sprot")

#Remove the column with Daphnia gene IDs
countsTable <- countsTableIn[2:25]

#Import grouping factor
#targets <- read.csv(file=args[2], row.names="sample")
targets <- read.csv(file="/Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_binned_Olympics.csv", row.names="sample")

#Setup a design matrix
group <- factor(paste(targets$treatment,targets$genotype,sep="."))

#The experimental design is parametrized with a one-way layout, 
# where one coefficient is assigned to each group
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

#Import gene sets
set1 <- read.csv(file="/Users/bamflappy/PfrenderLab/PA42_v4.1/geneSet1_KEGG_BASE_EXCISION_REPAIR.txt")
set2 <- read.csv(file="/Users/bamflappy/PfrenderLab/PA42_v4.1/geneSet2_KEGG_HOMOLOGOUS_RECOMBINATION.txt")
set3 <- read.csv(file="/Users/bamflappy/PfrenderLab/PA42_v4.1/geneSet3_KEGG_MISMATCH_REPAIR.txt")
set4 <- read.csv(file="/Users/bamflappy/PfrenderLab/PA42_v4.1/geneSet4_KEGG_NON_HOMOLOGOUS_END_JOINING.txt")
set5 <- read.csv(file="/Users/bamflappy/PfrenderLab/PA42_v4.1/geneSet5_KEGG_NUCLEOTIDE_EXCISION_REPAIR.txt")
set6 <- read.csv(file="/Users/bamflappy/PfrenderLab/PA42_v4.1/geneSet6_KEGG_P53_SIGNALING_PATHWAY.txt")

#Vectors containing indices of genes corresponding to gene sets
index1 <- countsTable[rownames(countsTable) %in% set1,]
index2 <- countsTable[rownames(countsTable) %in% set2,]
#etc...

#Individual gene set tests to check against batch and PR results
camera(countsTable, index1, design)
camera(countsTable, index2, design)
#etc...

#Batch gene set tests
#Inter-gene correlation will be estimated for each tested set
camera(countsTable, list(set1=index1,set2=index2), design, inter.gene.cor=NA)
#Specifying an inter-gene correlation of 0.01
camera(countsTable, list(set1=index1,set2=index2), design, inter.gene.cor=0.01)
#Specifying an inter-gene correlation of 0.05
camera(countsTable, list(set1=index1,set2=index2), design, inter.gene.cor=0.05)

#Pre-ranked (PR) gene set test version
#Fit the linear model using the design matrix
fit <- eBayes(lmFit(countsTable, design))
#Using moderated F-statistic
cameraPR(fit$F, list(set1=index1,set2=index2))
