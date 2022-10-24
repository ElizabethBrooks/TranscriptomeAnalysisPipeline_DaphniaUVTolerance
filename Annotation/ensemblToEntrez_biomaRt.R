
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
ensIDs=c("XM_046601202.1", "XR_006879527.1", "XR_006878626.1", "XR_006878982.1", "XM_046601405.1", "XM_046601407.1", "XM_046601408.1", "XM_046601117.1", "XR_006878744.1", "XR_006878741.1")

# map the Ensembl to the EntrezGene IDs 
getBM(attributes=c("ensembl_transcript_id", "entrezgene_id"), 
      filters = "ensembl_transcript_id", 
      values = ensIDs, 
      mart = ensembl)
