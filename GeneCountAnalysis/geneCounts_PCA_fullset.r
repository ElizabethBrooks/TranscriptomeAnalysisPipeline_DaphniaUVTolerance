#R script to compare two gene count matrices
#{ggfortify} lets {ggplot2} know how to interpret PCA objects
library(ggfortify)
#Retrieve gene count tables with gene IDs
gCount0 = read.csv("../../GeneCounts_Merged/merged_counts_fullset_transposed.csv", sep=",", row.names=1)
#Plot principal componants of PCA performed with prcomp
autoplot(prcomp(gCount0))