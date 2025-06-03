
# install pacakges, if necessary
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("biomaRt")

# load libraries
library("biomaRt")
#listMarts()

# specify the hsapiens mart
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

# check out filters
filters = listFilters(ensembl)
View(filters)

# check out attributes
attributes = listAttributes(ensembl)
View(attributes)

# import normalized gene count data
#inputTable <- read.csv(file=args[1], row.names="gene")[ ,args[2]:args[3]]
#inputTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/DEGenotypes/glmQLF_normalizedCounts.csv", row.names="gene", header=TRUE)[ ,1:24]
#View(inputTable)

# retrieve the Ensembl IDs
ensIDs=c("ENST00000369985")

# map the Ensembl to the EntrezGene IDs 
getBM(attributes=c("ensembl_transcript_id", "entrezgene_id"), 
      filters = "ensembl_transcript_id", 
      values = ensIDs, 
      mart = ensembl)
