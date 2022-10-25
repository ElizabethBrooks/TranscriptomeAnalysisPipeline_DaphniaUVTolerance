
# load libraries
library("biomaRt")

# specify the hsapiens mart
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

# check out filters
filters = listFilters(ensembl)
View(filters)

# check out attributes
attributes = listAttributes(ensembl)
View(attributes)

# retrieve the Ensembl IDs
ensIDs=c("XM_046580782.1")

# map the Ensembl to the EntrezGene IDs 
getBM(attributes=c("ensembl_transcript_id", "entrezgene_id"), 
      filters = "ensembl_transcript_id", 
      values = ensIDs, 
      mart = ensembl)
