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
library(statmod)

#Import gene count data for the Olympics
countsTable <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/geneCounts_cleaned_PA42_v4.1.csv", row.names="gene")[ ,1:24]
#Import grouping factor
targets <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/expDesign_Olympics.csv", row.names="sample")
#Set FDR cutoff
fdrCut=0.10

#Setup a design matrix
group <- factor(paste(targets$treatment,targets$genotype,sep="."))
#Create DGE list object
list <- DGEList(counts=countsTable,group=group)
colnames(list) <- targets$sample

#Retain genes only if it is expressed at a minimum level
keep <- filterByExpr(list)
list <- list[keep, , keep.lib.sizes=FALSE]

#Use TMM normalization to eliminate composition biases between libraries
list <- calcNormFactors(list)
#Write normalized counts to file
normList <- cpm(list, normalized.lib.sizes=TRUE)

#The experimental design is specified with a one-way layout, 
# where one coefficient is assigned to each group
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

#Estimate the NB dispersion
list <- estimateDisp(list, design, robust=TRUE)
#Estimate and plot the QL dispersions
fit <- glmQLFit(list, design, robust=TRUE)

#Test whether the average across all UV groups is equal to the average across
#all VIS groups, to examine the overall effect of treatment
con.UVvsVIS <- makeContrasts(UVvsVIS = (UV.E05 + UV.R2 + UV.Y023 + UV.Y05)/4
                             - (VIS.E05 + VIS.R2 + VIS.Y023 + VIS.Y05)/4,
                             levels=design)
#Look at genes expressed across all UV groups using QL F-test
test.anov.UVVIS <- glmQLFTest(fit, contrast=con.UVvsVIS)

#Test whether the average across all tolerant groups is equal to the average across
#all not tolerant groups, to examine the overall effect of tolerance
con.TvsN <- makeContrasts(TvsN = (UV.Y05 + VIS.Y05 + UV.E05 + VIS.E05)/4
                          - (UV.Y023 + VIS.Y023 + UV.R2 + VIS.R2)/4,
                          levels=design)
#Look at genes expressed across all UV groups using QL F-test
test.anov.TN <- glmQLFTest(fit, contrast=con.TvsN)

#Test whether there is an interaction effect
con.Inter <- makeContrasts(Inter = ((UV.E05 + UV.R2 + UV.Y023 + UV.Y05)/4
                                    - (VIS.E05 + VIS.R2 + VIS.Y023 + VIS.Y05)/4)
                           - ((UV.Y05 + VIS.Y05 + UV.E05 + VIS.E05)/4
                              - (UV.Y023 + VIS.Y023 + UV.R2 + VIS.R2)/4),
                           levels=design)
#Look at genes expressed across all UV groups using QL F-test
test.anov.Inter <- glmQLFTest(fit, contrast=con.Inter)


#GO enrichment
#Read in custom GO annotations
GOmaps <- readMappings(file="/home/mae/Documents/RNASeq_Workshop_ND/gene2GO_PA42_v4.1_transcripts.map",  sep="\t",  IDsep=",")

#Create named list of all genes (gene universe) and p-values. The gene universe is set to be
#the list of all genes contained in the gene2GO list of annotated genes.
DGE_results_table <- test.anov.UVVIS$table
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

#showGroupDensity for UV tolerance assocaited GO terms
CC_DNA_repair_complex <- "GO:1990391"
BP_DNA_integrity_checkpoint <- "GO:0031570"
BP_Response_UV <- "GO:0009411"
BP_Mitotic_cell_cycle_checkpoint <- "GO:0007093"
BP_Cellular_response_DNA_damage_stimulus <- "GO:0006974"
MF_Single_stranded_DNA_binding <- "GO:0003697"
MF_Damaged_DNA_binding <- "GO:0003684"

#pdf(file="/home/mae/Documents/RNASeq_Workshop_ND/All_Genes_TopSigGO_Density_OlympicsInteraction.pdf")
showGroupDensity(CC_GO_data, whichGO = CC_DNA_repair_complex, ranks = TRUE)
showGroupDensity(BP_GO_data, whichGO = BP_DNA_integrity_checkpoint, ranks = TRUE)
showGroupDensity(BP_GO_data, whichGO = BP_Response_UV, ranks = TRUE)
showGroupDensity(BP_GO_data, whichGO = BP_Mitotic_cell_cycle_checkpoint, ranks = TRUE)
showGroupDensity(BP_GO_data, whichGO = BP_Cellular_response_DNA_damage_stimulus, ranks = TRUE)
showGroupDensity(MF_GO_data, whichGO = MF_Single_stranded_DNA_binding, ranks = TRUE)
showGroupDensity(MF_GO_data, whichGO = MF_Damaged_DNA_binding, ranks = TRUE)
#dev.off()
