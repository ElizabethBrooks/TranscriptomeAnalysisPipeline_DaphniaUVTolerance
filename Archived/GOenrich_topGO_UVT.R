#if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
#BiocManager::install('_______')

#Load the libraries
library(filesstrings)
library(topGO)
library(GO.db)
library(reshape2)
library(ggplot2)
library(Rgraphviz)

#Import gene count data for the Olympics
countsTable <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/DEGs_UVBvsCntrl_Tcast_pvalues.csv")
#countsTable <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/Dmel_Svetec_2016_uvResponsiveGenes_pvalues.csv")


#GO enrichment
#Read in custom GO annotations
GOmaps <- readMappings(file="/home/mae/Documents/RNASeq_Workshop_ND/gene2GO_PA42_v4.1_transcripts.map",  sep="\t",  IDsep=",")

#Create named list of all genes (gene universe) and p-values. The gene universe is set to be
#the list of all genes contained in the gene2GO list of annotated genes.
list_genes <- as.numeric(countsTable$pval)
list_genes <- setNames(list_genes, rownames(countsTable))
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
