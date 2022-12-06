#if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
#BiocManager::install('______')


#Load the libraries 
library(filesstrings)
#library(topGO)
library(edgeR)
#library(GO.db)
#library(reshape2)
#library(ggplot2)
#library(Rgraphviz)

#Import gene count data & specify genotype
countsTable <- read.csv(file= "geneCounts_cleaned_PA42_v4.1.csv", 
                        row.names="gene")[ ,1:6]
genotype <- unlist(strsplit(colnames(countsTable)[1], '_'))[1] #this is extremely messy but it should work for genotype of any length

#Add grouping factor
group <- factor(c(rep("ctrl",3),rep("treat",3)))
#Create DGE list object
list <- DGEList(counts=countsTable,group=group)

#There is no purpose in analyzing genes that are not expressed in either 
# experimental condition, so genes are first filtered on expression levels
keep <- filterByExpr(list)
list <- list[keep, , keep.lib.sizes=FALSE]
#Calculate normalized factors
list <- calcNormFactors(list)

#Produce a matrix of pseudo-counts
#Estimate common dispersion and tagwise dispersions
list <- estimateDisp(list)

#Perform an exact test for treat vs ctrl
tested <- exactTest(list, pair=c("ctrl", "treat"))

#Create results table of DE genes
resultsTbl <- topTags(tested, n=nrow(tested$table))$table

#write csv file of DGE results table
write.csv(tested$table, file = paste(genotype, '_DGE_results.csv', sep = ''))
file.move(paste(getwd(), '/', genotype, '_DGE_results.csv', sep = ''), "/Users/bryanmichalek/Documents/Notre Dame/Spring 2021/Pfrender/DGE_results", overwrite = TRUE)
