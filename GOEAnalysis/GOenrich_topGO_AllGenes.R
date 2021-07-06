#if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
#BiocManager::install('_______')


#Load the libraries
library(filesstrings)
library(topGO)
library(edgeR)
library(GO.db)
library(reshape2)
library(ggplot2)
library(Rgraphviz)

#Import gene count data & specify genotype
countsTable <- read.csv(file= "/home/mae/Documents/RNASeq_Workshop_ND/geneCounts_cleaned_PA42_v4.1.csv", 
                        row.names="gene")[ ,31:36]
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

#GO enrichment
#Read in custom GO annotations
GOmaps <- readMappings(file="/home/mae/Documents/RNASeq_Workshop_ND/gene2GO_PA42_v4.1_transcripts.map",  sep='\t',  IDsep=',')

#Create named list of all genes (gene universe) and p-values. The gene universe is set to be
#the list of all genes contained in the gene2GO list of annotated genes.
DGE_results_table <- tested$table
list_genes <- as.numeric(DGE_results_table$PValue)
list_genes <- setNames(list_genes, rownames(DGE_results_table))
list_genes_filtered <- list_genes[names(list_genes) %in% names(GOmaps)]

#Create function to return list of interesting DE genes (0 == not significant, 1 == significant)
get_interesting_DE_genes <- function(geneUniverse){
  interesting_DE_genes <- rep(0, length(geneUniverse))
  for(i in 1:length(geneUniverse)){
    if(geneUniverse[i] < 0.05){
      interesting_DE_genes[i] = 1
    }
  }
  interesting_DE_genes <- setNames(interesting_DE_genes, names(geneUniverse))
  return(interesting_DE_genes)
}

#Create topGOdata objects for enrichment analysis (1 for each ontology)
BP_GO_data <- new('topGOdata', ontology = 'BP', allGenes = list_genes_filtered, 
                  geneSel = get_interesting_DE_genes, nodeSize = 10, annot = annFUN.gene2GO, 
                  gene2GO = GOmaps)
MF_GO_data <- new('topGOdata', ontology = 'MF', allGenes = list_genes_filtered, 
                  geneSel = get_interesting_DE_genes, nodeSize = 10, annot = annFUN.gene2GO, 
                  gene2GO = GOmaps)
CC_GO_data <- new('topGOdata', ontology = 'CC', allGenes = list_genes_filtered, 
                  geneSel = get_interesting_DE_genes, nodeSize = 10, annot = annFUN.gene2GO, 
                  gene2GO = GOmaps)

#Summary functions
#numGenes(BP_GO_data)
#length(sigGenes(BP_GO_data))
#numGenes(MF_GO_data)
#length(sigGenes(MF_GO_data))
#numGenes(CC_GO_data)
#length(sigGenes(CC_GO_data))

#PerformGO enrichment on topGOdata object
BP_GO_results <- runTest(BP_GO_data, statistic = 'Fisher')
MF_GO_results <- runTest(MF_GO_data, statistic = 'Fisher')
CC_GO_results <- runTest(CC_GO_data, statistic = 'Fisher')

#If you want to see names of GO terms (can filter for only significant ones if you want...etc.)
#head(names(BP_GO_results@score))
#geneData(BP_GO_results)

#Visualization/plot stuff
    #Store p-values as named list... ('score(x)' or 'x@score' returns named list of p-val's 
    #where names are the GO terms)
pval_BP_GO <- score(BP_GO_results)
pval_MF_GO <- score(MF_GO_results)
pval_CC_GO <- score(CC_GO_results)

    #plot histogram to see range of p-values
par(mfrow=c(3, 1))
hist(pval_BP_GO, 35, xlab = 'p-values', main = 'Range of BP GO term p-values')
hist(pval_MF_GO, 35, xlab = 'p-values', main = 'Range of BP GO term p-values')
hist(pval_CC_GO, 35, xlab = 'p-values', main = 'Range of BP GO term p-values')
dev.off()

    #GenTable to get statistics on GO terms
list_BP_GO_terms <- usedGO(BP_GO_data)
list_MF_GO_terms <- usedGO(MF_GO_data)
list_CC_GO_terms <- usedGO(CC_GO_data)

BP_GO_results_table <- GenTable(BP_GO_data, weightFisher = BP_GO_results, orderBy = 'weightFisher', 
                                topNodes = length(list_BP_GO_terms))
MF_GO_results_table <- GenTable(MF_GO_data, weightFisher = MF_GO_results, orderBy = 'weightFisher', 
                                topNodes = length(list_MF_GO_terms))
CC_GO_results_table <- GenTable(CC_GO_data, weightFisher = CC_GO_results, orderBy = 'weightFisher', 
                                topNodes = length(list_CC_GO_terms))

    #If you just want the significant GO terms
BP_sigGO_results_table <- BP_GO_results_table[BP_GO_results_table$weightFisher <= 0.05, ]
MF_sigGO_results_table <- MF_GO_results_table[MF_GO_results_table$weightFisher <= 0.05, ]
CC_sigGO_results_table <- CC_GO_results_table[CC_GO_results_table$weightFisher <= 0.05, ]

    #showGroupDensity example for most significant GO term
BP_topSigGO_ID <- BP_GO_results_table[1, 'GO.ID']
MF_topSigGO_ID <- MF_GO_results_table[1, 'GO.ID']
CC_topSigGO_ID <- CC_GO_results_table[1, 'GO.ID']

pdf(file = paste(genotype, 'All_Genes','TopSigGO_Density.pdf', sep = '_'))
showGroupDensity(BP_GO_data, whichGO = BP_topSigGO_ID, ranks = TRUE)
showGroupDensity(MF_GO_data, whichGO = MF_topSigGO_ID, ranks = TRUE)
showGroupDensity(CC_GO_data, whichGO = CC_topSigGO_ID, ranks = TRUE)
dev.off()

    #printGraph to plot subgraphs induced by the most significant GO terms...saves to a file
printGraph(BP_GO_data, BP_GO_results, firstSigNodes = 5, 
           fn.prefix = paste(genotype, 'all_genes','BP_GO', sep = '_'), useInfo = 'all', pdfSW = TRUE)
printGraph(MF_GO_data, MF_GO_results, firstSigNodes = 5, 
           fn.prefix = paste(genotype, 'all_genes', 'MF_GO', sep = '_'), useInfo = 'all', pdfSW = TRUE)
printGraph(CC_GO_data, CC_GO_results, firstSigNodes = 5, 
           fn.prefix = paste(genotype, 'all_genes', 'CC_GO', sep = '_'), useInfo = 'all', pdfSW = TRUE)


#Move some files around just to help myself stay organized
#file.move(paste(getwd(), '/', genotype, "_All_Genes_TopSigGO_Density.pdf", sep = ''), paste("/Users/bryanmichalek/Documents/Notre Dame/Spring 2021/Pfrender/GO_Custom_Annotations_Results/", genotype, "_all_genes", sep = ''))
#file.move(paste(getwd(), '/', genotype, "_all_genes_BP_GO_weight01_5_all.pdf", sep = ''), paste("/Users/bryanmichalek/Documents/Notre Dame/Spring 2021/Pfrender/GO_Custom_Annotations_Results/", genotype, "_all_genes", sep = ''))
#file.move(paste(getwd(), '/', genotype, "_all_genes_MF_GO_weight01_5_all.pdf", sep = ''), paste("/Users/bryanmichalek/Documents/Notre Dame/Spring 2021/Pfrender/GO_Custom_Annotations_Results/", genotype, "_all_genes", sep = ''))
#file.move(paste(getwd(), '/', genotype, "_all_genes_CC_GO_weight01_5_all.pdf", sep = ''), paste("/Users/bryanmichalek/Documents/Notre Dame/Spring 2021/Pfrender/GO_Custom_Annotations_Results/", genotype, "_all_genes", sep = ''))

#showGroupDensity for UV tolerance assocaited GO terms
CC_DNA_repair_complex <- "GO:1990391"
BP_DNA_integrity_checkpoint <- "GO:0031570"
BP_Response_UV <- "GO:0009411"
BP_Mitotic_cell_cycle_checkpoint <- "GO:0007093"
BP_Cellular_response_DNA_damage_stimulus <- "GO:0006974"
MF_Single_stranded_DNA_binding <- "GO:0003697"
MF_Damaged_DNA_binding <- "GO:0003684"

pdf(file="/home/mae/Documents/RNASeq_Workshop_ND/GSEA_UVT_OlympicsTolerance.pdf")
showGroupDensity(CC_GO_data, whichGO = CC_DNA_repair_complex, ranks = TRUE)
showGroupDensity(BP_GO_data, whichGO = BP_DNA_integrity_checkpoint, ranks = TRUE)
showGroupDensity(BP_GO_data, whichGO = BP_Response_UV, ranks = TRUE)
showGroupDensity(BP_GO_data, whichGO = BP_Mitotic_cell_cycle_checkpoint, ranks = TRUE)
showGroupDensity(BP_GO_data, whichGO = BP_Cellular_response_DNA_damage_stimulus, ranks = TRUE)
showGroupDensity(MF_GO_data, whichGO = MF_Single_stranded_DNA_binding, ranks = TRUE)
showGroupDensity(MF_GO_data, whichGO = MF_Damaged_DNA_binding, ranks = TRUE)
dev.off()