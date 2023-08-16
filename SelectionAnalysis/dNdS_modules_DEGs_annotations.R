#!/usr/bin/env Rscript

# load libraries
library(tidyr)

# turn off scientific notation
#options(scipen = 999)

# turn on scientific notation
#options(scipen = 0)

# set the working directory
workingDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/SelectionTests/Genotypes"
setwd(workingDir)

# retrieve subset tag
set <- "OLYM"

# set the minimum module size
minModSize <- "30"

# retrieve WGCNA directory
inDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Genotypes"

# set the full subset tag name
tag <- paste(set, minModSize, sep="_")

# Load the expression and trait data saved in the first part
importFile <- paste(set, "dataInput.RData", sep="-")
importFile <- paste(inDir, importFile, sep="/")
lnames1 = load(file = importFile)

# Load network data saved in the second part
importFile <- paste(tag, "networkConstruction-stepByStep.RData", sep="-")
importFile <- paste(inDir, importFile, sep="/")
lnames2 = load(file = importFile)

# create list of module colors mapped to numbers
numMods <- length(unique(moduleColors))
colorTable <- data.frame(
  color = unique(moduleColors),
  number = seq(from = 1, to = numMods, by = 1)
)

# initialize module data frame
resultsTable <- data.frame(
  gene = character(),
  color = character(),
  number = numeric()
)

# match gene IDs with module colors
for(i in 1:numMods){
  gene <- names(datExpr)[moduleColors==colorTable[i,1]]
  color <- rep(colorTable[i,1], length(gene))
  number <- rep(colorTable[i,2], length(gene))
  moduleData <- cbind(gene, color, number)
  resultsTable <- rbind(resultsTable, moduleData)
}


# retrieve dN dS values
dNdSTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/SelectionTests/Pulex_Olympics_kaksResults.csv", row.names="geneID")

# remove NAs
dNdSSubset <- na.omit(dNdSTable)

# subset dN dS values to remove outliers
#dNdSSubset <- dNdSSubset[dNdSSubset$dNdS < 99,]


# merging data frames
# https://sparkbyexamples.com/r-programming/r-join-data-frames-with-examples/#full-outer-join

# add geneID column
geneID <- row.names(dNdSSubset)
dNdSSubset <- cbind(geneID,dNdSSubset)
colnames(resultsTable)[1] ="geneID"

# full outer join data frames
resultsTable <- merge(x = resultsTable, y = dNdSSubset, 
                      by = "geneID", all=TRUE)

# set color and numner tags for NAs
resultsTable$color <- resultsTable$color %>% replace_na('None')
resultsTable$number <- resultsTable$number %>% replace_na('None')

# remove NAs
resultsTable <- na.omit(resultsTable)

# add column for identifying mode of selection
resultsTable$Selection[resultsTable$dNdS == 99] <- "Error"
resultsTable$Selection[resultsTable$dNdS > 1 & resultsTable$dNdS < 99] <- "Positive"
resultsTable$Selection[resultsTable$dNdS < 1] <- "Negative"


# retrieve DEGs
interactionTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Genotypes/glmQLF_2WayANOVA_interaction_topTags_LFC1.2.csv", row.names="gene")
treatmentTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Genotypes/glmQLF_2WayANOVA_treatment_topTags_LFC1.2.csv", row.names="gene")
genotypesTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Genotypes/glmQLF_2WayANOVA_tolerance_topTags_LFC1.2.csv", row.names="gene")

# keep only sig
interactionSig <- interactionTable[interactionTable$FDR < 0.05,]
treatmentSig <- treatmentTable[treatmentTable$FDR < 0.05,]
genotypesSig <- genotypesTable[genotypesTable$FDR < 0.05,]

# add effect tags
interactionSig$Effect <- "Interaction"
treatmentSig$Effect <- "Treatment"
genotypesSig$Effect <- "Tolerance"

# add geneID column
geneID <- row.names(dNdSSubset)
dNdSSubset <- cbind(geneID,dNdSSubset)
geneID <- row.names(interactionSig)
interactionSig <- cbind(geneID,interactionSig)
geneID <- row.names(treatmentSig)
treatmentSig <- cbind(geneID,treatmentSig)
geneID <- row.names(genotypesSig)
genotypesSig <- cbind(geneID,genotypesSig)

# keep necessary columns
dNdSSubset <- dNdSSubset[,c("geneID","dN","dS","dNdS")]
treatmentSubset <- treatmentSig[,c("geneID","logFC","FDR","Effect")]
genotypesSubset <- genotypesSig[,c("geneID","logFC","FDR","Effect")]
interactionSubset <- interactionSig[,c("geneID","logFC","FDR","Effect")]


# import annotation data
annotationTable <- read.delim(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/email/KAP4_NCBI_functional_annotation.txt", row.names="Gene_locus")

# add geneID column
geneID <- row.names(annotationTable)
annotationTable <- cbind(geneID,annotationTable)


# merging data frames
# https://sparkbyexamples.com/r-programming/r-join-data-frames-with-examples/#full-outer-join

# combine all tables
effectTable <- rbind(interactionSubset, treatmentSubset, genotypesSubset)

# full outer join
fullTable <- merge(x = resultsTable, y = effectTable, 
                   by = "geneID", all=TRUE)


# force scientific notation
formatC(x, format = "e", digits = 2) 

# convert columns to character type
fullTable$logFC <- as.character(fullTable$logFC)
fullTable$FDR <- as.character(fullTable$FDR)

# set tags for NAs
fullTable$Effect <- fullTable$Effect %>% replace_na('None')
fullTable$logFC <- fullTable$logFC %>% replace_na('None')
fullTable$FDR <- fullTable$FDR %>% replace_na('None')

# remove NAs
#plotTable <- na.omit(plotTable)



