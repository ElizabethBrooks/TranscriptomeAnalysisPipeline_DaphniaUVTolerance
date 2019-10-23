#R script to compare two gene count matrices
#Retrieve gene count tables
gCount1 = read.csv("../../GeneCounts_Merged/merged_counts_legacy_tagged.csv", sep=",", row.names=1)
gCount2 = read.csv("../../GeneCounts_Merged/merged_counts_subset_tagged.csv", sep=",", row.names=1)
#Convert gene count tables to matrices
gCount1M <- as.matrix(gCount1)
gCount2M <- as.matrix(gCount2)
#Determine matrix equality
eq <- gCount1M==gCount2M
#Get the percentage of non equal values
round(sum(!eq)/length(eq)*100, 2)