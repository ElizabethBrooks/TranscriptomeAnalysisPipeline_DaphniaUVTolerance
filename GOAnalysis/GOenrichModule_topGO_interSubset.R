#if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
#BiocManager::install("")

#Load the libraries
#library(filesstrings)
library(topGO)
library(edgeR)
#library(GO.db)
#library(reshape2)
library(ggplot2)
library(Rgraphviz)
#library(statmod)

#Set working directory
workingDir="/Users/bamflappy/PfrenderLab/WGCNA_PA42_v4.1"
setwd(workingDir); 

# Load the expression and trait data saved in the first part
#lnames1 = load(file = "PA42_v4.1_dataInputInter.RData");
# Load network data saved in the second part.
lnames2 = load(file = "PA42_v4.1_networkConstructionInter_auto_threshold8_signedNowick.RData");


#Import gene count data for the Olympics
countsTable <- read.csv(file="/Users/bamflappy/PfrenderLab/PA42_v4.1/geneCounts_cleaned_PA42_v4.1.csv", row.names="gene")[ ,1:24]
#Import grouping factor
targets <- read.csv(file="/Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_Olympics.csv", row.names="sample")


#Setup a design matrix
group <- factor(paste(targets$treatment,targets$genotype,sep="."))
#Create DGE list object
list <- DGEList(counts=countsTable,group=group)
colnames(list) <- rownames(targets)

#Retain genes only if it is expressed at a minimum level
keep <- filterByExpr(list)
list <- list[keep, , keep.lib.sizes=FALSE]

#Use TMM normalization to eliminate composition biases between libraries
list <- calcNormFactors(list)
normList <- cpm(list, normalized.lib.sizes=TRUE)

#The experimental design is specified with a one-way layout, 
# where one coefficient is assigned to each group
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

#Estimate the NB dispersion
list <- estimateDisp(list, design, robust=TRUE)
#Estimate and plot the QL dispersions
fit <- glmQLFit(list, design, robust=TRUE)

#Test whether there is an interaction effect
con.Inter <- makeContrasts(Inter = ((UV.E05 + UV.R2 + UV.Y023 + UV.Y05)/4
                                    - (VIS.E05 + VIS.R2 + VIS.Y023 + VIS.Y05)/4)
                           - ((UV.Y05 + VIS.Y05 + UV.E05 + VIS.E05)/4
                              - (UV.Y023 + VIS.Y023 + UV.R2 + VIS.R2)/4),
                           levels=design)
#Look at genes expressed across all treatment groups using QL F-test
test.anov.Inter <- glmQLFTest(fit, contrast=con.Inter)
summary(decideTests(test.anov.Inter))


#GO enrichment
#Read in custom GO annotations
GOmaps <- readMappings(file="/Users/bamflappy/PfrenderLab/PA42_v4.1/gene2GO_PA42_v4.1_transcripts.map",  sep="\t",  IDsep=",")

#Match gene to its geneID, module color, and module number
moduleList <- data.frame(geneID=names(moduleLabels), color=moduleColors, number=unname(moduleLabels), stringsAsFactors=FALSE)
DGE_results_table <- test.anov.Inter$table
DGE_results_table$geneID <- rownames(DGE_results_table)
DGE_results_table <- merge(DGE_results_table, moduleList, by = c("geneID"), all = TRUE)

#Replace NAs with an unused module number
highest <- max(unique(unname(moduleLabels)))+1
DGE_results_table[is.na(DGE_results_table)] <- highest

#Create named list of all genes (gene universe) and p-values. The gene universe is set to be
#the list of all genes contained in the gene2GO list of annotated genes.
list_genes <- as.numeric(DGE_results_table$number)
list_genes <- setNames(list_genes, DGE_results_table$geneID)
list_genes_filtered <- list_genes[names(list_genes) %in% names(GOmaps)]

#Create data frame to hold the module number, color, total genes, significant genes, 
# and top significant BP, MF, and CC GO terms
numMods <- length(unique(unname(moduleLabels)))+1
moduleBPResults = data.frame(matrix(ncol = 10, nrow = numMods))
colnames(moduleBPResults) <- c("number","color","total","sig","topID","topTerm","annotated","topSig","topExpected","topWeight")
moduleMFResults = data.frame(matrix(ncol = 10, nrow = numMods))
colnames(moduleMFResults) <- c("number","color","total","sig","topID","topTerm","annotated","topSig","topExpected","topWeight")
moduleCCResults = data.frame(matrix(ncol = 10, nrow = numMods))
colnames(moduleCCResults) <- c("number","color","total","sig","topID","topTerm","annotated","topSig","topExpected","topWeight")

#Loop through each module
lowest <- min(unique(unname(moduleLabels)))
var <- 0
for(j in lowest:highest){
  #Create function to return list of interesting DE genes (0 == not significant, 1 == significant)
  get_interesting_DE_genes <- function(geneUniverse){
    interesting_DE_genes <- rep(0, length(geneUniverse))
    for(i in 1:length(geneUniverse)){
      if(geneUniverse[i] == j){
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
  
  #PerformGO enrichment on topGOdata object
  BP_GO_results <- runTest(BP_GO_data, statistic = 'Fisher')
  MF_GO_results <- runTest(MF_GO_data, statistic = 'Fisher')
  CC_GO_results <- runTest(CC_GO_data, statistic = 'Fisher')
  
  #GenTable to get statistics on GO terms
  list_BP_GO_terms <- usedGO(BP_GO_data)
  list_MF_GO_terms <- usedGO(MF_GO_data)
  list_CC_GO_terms <- usedGO(CC_GO_data)
  
  #Create summary data table
  BP_GO_results_table <- GenTable(BP_GO_data, weightFisher = BP_GO_results, orderBy = 'weightFisher', 
                                  topNodes = length(list_BP_GO_terms))
  MF_GO_results_table <- GenTable(MF_GO_data, weightFisher = MF_GO_results, orderBy = 'weightFisher', 
                                  topNodes = length(list_MF_GO_terms))
  CC_GO_results_table <- GenTable(CC_GO_data, weightFisher = CC_GO_results, orderBy = 'weightFisher', 
                                  topNodes = length(list_CC_GO_terms))
  
  #Set row number
  var <- var+1
  
  #Summary BP functions
  moduleBPResults$number[var] <- j
  if(j < highest){
    moduleBPResults$color[var] <- head(moduleList[moduleList$number %in% j,2],1)
  }
  moduleBPResults$total[var] <- numGenes(BP_GO_data)
  moduleBPResults$sig[var] <- length(sigGenes(BP_GO_data))
  moduleBPResults$topID[var] <- BP_GO_results_table[1,1]
  moduleBPResults$topTerm[var] <- BP_GO_results_table[1,2]
  moduleBPResults$annotated[var] <- BP_GO_results_table[1,3]
  moduleBPResults$topSig[var] <- BP_GO_results_table[1,4]
  moduleBPResults$topExpected[var] <- BP_GO_results_table[1,5]
  moduleBPResults$topWeight[var] <- BP_GO_results_table[1,6]
  
  #Summary MF functions
  moduleMFResults$number[var] <- j
  if(j < highest){
    moduleMFResults$color[var] <- head(moduleList[moduleList$number %in% j,2],1)
  }
  moduleMFResults$total[var] <- numGenes(MF_GO_data)
  moduleMFResults$sig[var] <- length(sigGenes(MF_GO_data))
  moduleMFResults$topID[var] <- MF_GO_results_table[1,1]
  moduleMFResults$topTerm[var] <- MF_GO_results_table[1,2]
  moduleMFResults$annotated[var] <- MF_GO_results_table[1,3]
  moduleMFResults$topSig[var] <- MF_GO_results_table[1,4]
  moduleMFResults$topExpected[var] <- MF_GO_results_table[1,5]
  moduleMFResults$topWeight[var] <- MF_GO_results_table[1,6]
  
  #Summary CC functions
  moduleCCResults$number[var] <- j
  if(j < highest){
    moduleCCResults$color[var] <- head(moduleList[moduleList$number %in% j,2],1)
  }
  moduleCCResults$total[var] <- numGenes(CC_GO_data)
  moduleCCResults$sig[var] <- length(sigGenes(CC_GO_data))
  moduleCCResults$topID[var] <- CC_GO_results_table[1,1]
  moduleCCResults$topTerm[var] <- CC_GO_results_table[1,2]
  moduleCCResults$annotated[var] <- CC_GO_results_table[1,3]
  moduleCCResults$topSig[var] <- CC_GO_results_table[1,4]
  moduleCCResults$topExpected[var] <- CC_GO_results_table[1,5]
  moduleCCResults$topWeight[var] <- CC_GO_results_table[1,6]
}


#Write the resulting tables to files
write.table(moduleBPResults, file="GOAnalysis/moduleTopGO_BPResults_interSubset_signedNowick.csv", sep=",", row.names=FALSE)
write.table(moduleMFResults, file="GOAnalysis/moduleTopGO_MFResults_interSubset_signedNowick.csv", sep=",", row.names=FALSE)
write.table(moduleCCResults, file="GOAnalysis/moduleTopGO_CCResults_interSubset_signedNowick.csv", sep=",", row.names=FALSE)
