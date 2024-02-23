#!/usr/bin/env Rscript

# R script to perform statistical & annotation analysis of gene count tables using edgeR GLM

# https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7873980/
# https://support.bioconductor.org/p/132926/
# https://support.bioconductor.org/p/106608/

# install edgeR and statmod, this should only need to be done once
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("edgeR")
#install.packages("statmod")

# turn off scientific notation
options(scipen = 999)

# load the edgeR library
library(edgeR)
library(statmod)
library(ggplot2)
library(ggrepel)
library(ggVennDiagram)
library(rcartocolor)
library(topGO)
library(eulerr)
library(stringr)
library(tidyr)
library(gplots)

# plotting palettes
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
# https://github.com/Nowosad/rcartocolor
plotColors <- carto_pal(12, "Safe")
plotColorSubset <- c(plotColors[4], plotColors[5], plotColors[6])

# retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)

# set working directory
workingDir="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Genotypes"

# import gene count data
inputTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/GeneCountsAnalyzed/Formatted/cleaned.csv", row.names="gene")[ ,1:24]

# import annotations data
annotTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/email/geneToAnnotation_map.csv", row.names="gene")

# trim the data tables of quantification stats lines
countsTable <- head(inputTable, - 5)
annotations <- head(annotTable, - 5)

# add gene ID to empty cells
for(i in 1:length(annotations)){
  annotations[annotations$annotation == "",] <- rownames(annotations)[i]
  annotations[annotations$annotation == " ",] <- rownames(annotations)[i]
}

# import grouping factor
targets <- read.csv(file="/Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_OlympicsGenotypes.csv", row.names="sample")

# setup a design matrix
group <- factor(paste(targets$treatment,targets$genotype,sep="."))
#cbind(targets,Group=group)
#Create DGE list object
list <- DGEList(counts=countsTable,group=group)
colnames(list) <- rownames(targets)

# retain genes only if it is expressed at a minimum level
keep <- filterByExpr(list)
summary(keep)
list <- list[keep, , keep.lib.sizes=FALSE]

# use TMM normalization to eliminate composition biases
# between libraries
list <- calcNormFactors(list)
#list$samples
#Write normalized counts to file
normList <- cpm(list, normalized.lib.sizes=TRUE)

# write log transformed normalized counts to file
normListLog <- cpm(list, normalized.lib.sizes=TRUE, log=TRUE)


##
#The experimental design is parametrized with a one-way layout, 
# where one coefficient is assigned to each group
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
#design

#Next, the NB dispersion is estimated
list <- estimateDisp(list, design, robust=TRUE)
#list$common.dispersion

#Now, estimate and plot the QL dispersions
fit <- glmQLFit(list, design, robust=TRUE)
#head(fit$coefficients)

# view column order
colnames(fit)


# testing explicit nested contrasts
con.all.nest <- makeContrasts(treatment = (UV.E05 + UV.R2 + UV.Y023 + UV.Y05)/4
                              - (VIS.E05 + VIS.R2 + VIS.Y023 + VIS.Y05)/4,
                              tolerance = (UV.Y023 + UV.Y05 + VIS.Y023 + VIS.Y05)/4
                              - (UV.E05 + UV.R2 + VIS.E05 + VIS.R2)/4,
                              levels=design)
# treatment
treat.anov.treatment <- glmTreat(fit, contrast=con.all.nest[,"treatment"], lfc=log2(1.2))
summary(decideTests(treat.anov.treatment))
# tolerance
treat.anov.tolerance <- glmTreat(fit, contrast=con.all.nest[,"tolerance"], lfc=log2(1.2))
summary(decideTests(treat.anov.tolerance))
# interaction
treat.anov.Inter <- glmTreat(fit, contrast=c(con.all.nest[,"treatment"]-con.all.nest[,"tolerance"]), lfc=log2(1.2))
summary(decideTests(treat.anov.Inter))


# export tables of DE genes
#Write tags table of DE genes to file
tagsTblANOVATreatment <- topTags(treat.anov.treatment, n=nrow(treat.anov.treatment$table), adjust.method="fdr")$table

#Write tags table of DE genes to file
tagsTblANOVATolerance <- topTags(treat.anov.tolerance, n=nrow(treat.anov.tolerance$table), adjust.method="fdr")$table

#Generate table of DE genes
tagsTblANOVAInter <- topTags(treat.anov.Inter, n=nrow(treat.anov.Inter$table), adjust.method="fdr")$table


# add column for identifying direction of DE gene expression
tagsTblANOVATreatment$topDE <- "NA"
# identify significantly up DE genes
tagsTblANOVATreatment$topDE[tagsTblANOVATreatment$logFC > 1 & tagsTblANOVATreatment$FDR < 0.05] <- "UP"
# identify significantly down DE genes
tagsTblANOVATreatment$topDE[tagsTblANOVATreatment$logFC < -1 & tagsTblANOVATreatment$FDR < 0.05] <- "DOWN"
# create volcano plot with labels
labelSetTreatment <- tagsTblANOVATreatment[tagsTblANOVATreatment$topDE == "UP" | tagsTblANOVATreatment$topDE == "DOWN",]
# identify significantly DE genes by FDR
tagsTblANOVATreatment.glm_keep <- tagsTblANOVATreatment$FDR < 0.05
# create filtered results table of DE genes
tagsTblANOVATreatment.filtered <- tagsTblANOVATreatment[tagsTblANOVATreatment.glm_keep,]

# add column for identifying direction of DE gene expression
tagsTblANOVATolerance$topDE <- "NA"
# identify significantly up DE genes
tagsTblANOVATolerance$topDE[tagsTblANOVATolerance$logFC > 1 & tagsTblANOVATolerance$FDR < 0.05] <- "UP"
# identify significantly down DE genes
tagsTblANOVATolerance$topDE[tagsTblANOVATolerance$logFC < -1 & tagsTblANOVATolerance$FDR < 0.05] <- "DOWN"
# create volcano plot with labels
labelSetTolerance <- tagsTblANOVATolerance[tagsTblANOVATolerance$topDE == "UP" | tagsTblANOVATolerance$topDE == "DOWN",]
# identify significantly DE genes by FDR
tagsTblANOVATolerance.glm_keep <- tagsTblANOVATolerance$FDR < 0.05
# create filtered results table of DE genes
tagsTblANOVATolerance.filtered <- tagsTblANOVATolerance[tagsTblANOVATolerance.glm_keep,]

# add column for identifying direction of DE gene expression
tagsTblANOVAInter$topDE <- "NA"
# identify significantly up DE genes
tagsTblANOVAInter$topDE[tagsTblANOVAInter$logFC > 1 & tagsTblANOVAInter$FDR < 0.05] <- "UP"
# identify significantly down DE genes
tagsTblANOVAInter$topDE[tagsTblANOVAInter$logFC < -1 & tagsTblANOVAInter$FDR < 0.05] <- "DOWN"
# create volcano plot with labels
labelSetInteraction <- tagsTblANOVAInter[tagsTblANOVAInter$topDE == "UP" | tagsTblANOVAInter$topDE == "DOWN",]
# identify significantly DE genes by FDR
tagsTblANOVAInter.glm_keep <- tagsTblANOVAInter$FDR < 0.05
# create filtered results table of DE genes
tagsTblANOVAInter.filtered <- tagsTblANOVAInter[tagsTblANOVAInter.glm_keep,]

# venn diagram
# retrieve set of DE gene names for hours contrast
geneSet_treatment <- rownames(tagsTblANOVATreatment.filtered)
# retrieve set of DE gene names for hours contrast
geneSet_tolerance <- rownames(tagsTblANOVATolerance.filtered)
# retrieve set of DE gene names for interaction contrast
geneSet_interaction <- rownames(tagsTblANOVAInter.filtered)
# create combined glm_list of DE gene names
glm_list_venn <- list(Treatment = c(geneSet_treatment), 
                      Tolerance = c(geneSet_tolerance),
                      Interaction = c(geneSet_interaction))
# create venn diagram
##jpeg("glmQLF_2WayANOVA_venn_LFC1.2.jpg")
ggVennDiagram(glm_list_venn, label_alpha=0.25, category.names = c("Treatment","Tolerance","Interaction")) +
  scale_colour_discrete(type = plotColorSubset)
##dev.off()
# create venn lists
vennList <- gplots::venn(glm_list_venn, show.plot = FALSE)
# retrieve intersections
listAtt <- attributes(vennList)$intersections
listAtt

# heatmaps
# heatmap data
logcounts = cpm(list, log=TRUE)
# view DGE genes
# subset counts and annotations tables by DE gene set
DGESubset_treatment <- tagsTblANOVATreatment.filtered[!grepl("NA", tagsTblANOVATreatment.filtered$topDE),]
annotationsList_treatment <- subset(annotations,
                          grepl(
                            paste0(rownames(DGESubset_treatment), collapse = "|"),
                            rownames(annotations),
                            ignore.case = TRUE
                          )
)
logcountsSubset_treatment <- subset(logcounts,
                          grepl(
                            paste0(rownames(DGESubset_treatment), collapse = "|"),
                            rownames(logcounts),
                            ignore.case = TRUE
                          )
)
##jpeg("glmQLF_treatment_heatmap.jpg")
heatmap(logcountsSubset_treatment, main= "Heatmap of Treatment Effect DGE", margins = c(8, 8), labRow = annotationsList_treatment$annotation)
##dev.off()

# view tolerance DGE genes
# subset counts table by DE gene set
DGESubset_tolerance <- tagsTblANOVATolerance.filtered[!grepl("NA", tagsTblANOVATolerance.filtered$topDE),]
annotationsList_tolerance <- subset(annotations,
                          grepl(
                            paste0(rownames(DGESubset_tolerance), collapse = "|"),
                            rownames(annotations),
                            ignore.case = TRUE
                          )
)
logcountsSubset_tolerance <- subset(logcounts,
                          grepl(
                            paste0(rownames(DGESubset_tolerance), collapse = "|"),
                            rownames(logcounts),
                            ignore.case = TRUE
                          )
)
#jpeg("glmQLF_tolerance_heatmap.jpg")
heatmap(logcountsSubset_tolerance, main= "Heatmap of Tolerance Effect DGE", margins = c(8, 8), labRow = annotationsList_tolerance$annotation)
#dev.off()

# view interaction DGE genes
# subset counts table by DE gene set
DGESubset_interaction <- tagsTblANOVAInter.filtered[!grepl("NA", tagsTblANOVAInter.filtered$topDE),]
annotationsList_interaction <- subset(annotations,
                          grepl(
                            paste0(rownames(DGESubset_interaction), collapse = "|"),
                            rownames(annotations),
                            ignore.case = TRUE
                          )
)
logcountsSubset_interaction <- subset(logcounts,
                          grepl(
                            paste0(rownames(DGESubset_interaction), collapse = "|"),
                            rownames(logcounts),
                            ignore.case = TRUE
                          )
)
#jpeg("glmQLF_interaction_heatmap.jpg")
heatmap(logcountsSubset_interaction, main= "Heatmap of Interaction Effect DGE", margins = c(8, 8), labRow = annotationsList_interaction$annotation)
#dev.off()

# view intersection DGE genes
# $`Treatment:Tolerance:Interaction`
# [1] "gene-LOC124205968" "gene-LOC124208098" "gene-LOC124193441" "gene-LOC124204390" "gene-LOC124210001"
# [6] "gene-LOC124198796" "gene-LOC124204460" "gene-LOC124208917" "gene-LOC124196372" "gene-LOC124203684"
# [11] "gene-LOC124200027" "gene-LOC124209663" "gene-LOC124190261" "gene-LOC124204974"
# subset counts table by DE gene set
DGESubset_intersection <- c("gene-LOC124205968", "gene-LOC124208098", "gene-LOC124193441", "gene-LOC124204390", "gene-LOC124210001",
                           "gene-LOC124198796", "gene-LOC124204460", "gene-LOC124208917", "gene-LOC124196372", "gene-LOC124203684",
                           "gene-LOC124200027", "gene-LOC124209663", "gene-LOC124190261", "gene-LOC124204974")
logcountsSubset_intersection <- subset(logcounts,
                                      grepl(
                                        paste0(DGESubset_intersection, collapse = "|"),
                                        rownames(logcounts),
                                        ignore.case = TRUE
                                      )
)
##jpeg("glmQLF_intersection_heatmap.jpg")
heatmap(logcountsSubset_intersection, main= "Heatmap of Intersection of Effects DGE", margins = c(8, 1))
##dev.off()


# interesting GO
# set the working directory
workingDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Genotypes/GOAnalysis_OLYM_30"
setwd(workingDir)

# import modules GO data
modulesNames <- list.files(path="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Genotypes/GOAnalysis_OLYM_30/", pattern = "_BP_sigGO_terms.csv.flt$")
modulesNames <- str_remove(modulesNames, "_BP_sigGO_terms.csv.flt")
modulesFiles <- list.files(path="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Genotypes/GOAnalysis_OLYM_30/", pattern = "_BP_sigGO_terms.csv.flt$", full.names = TRUE)
modulesGO <- lapply(modulesFiles, read.csv)
names(modulesGO) <- modulesNames

# retrieve subset tag
set <- "OLYM"

# set the minimum module size
minModSize <- "30"

# retrieve WGCNA directory
inDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Genotypes"

# set the full subset tag name
tag <- paste(set, minModSize, sep="_")

# Load the expression and trait data saved in the first part
importFile <- paste(set, "dataInput.RData", sep="-")
importFile <- paste(inDir, importFile, sep="/")
lnames1 = load(file = importFile)

# Load network data saved in the second part
importFile <- paste(tag, "networkConstruction-stepByStep.RData", sep="-")
importFile <- paste(inDir, importFile, sep="/")
lnames2 = load(file = importFile)

# create list of module colors mapped to numbers
numMods <- length(unique(moduleColors))
colorTable <- data.frame(
  color = unique(moduleColors),
  number = seq(from = 1, to = numMods, by = 1)
)

# initialize module data frame
resultsTable <- data.frame(
  gene = character(),
  color = character(),
  number = numeric()
)

# match gene IDs with module colors
for(i in 1:numMods){
  gene <- names(datExpr)[moduleColors==colorTable[i,1]]
  color <- rep(colorTable[i,1], length(gene))
  number <- rep(colorTable[i,2], length(gene))
  moduleData <- cbind(gene, color, number)
  resultsTable <- rbind(resultsTable, moduleData)
}

# setup data frame of module sizes
modSizes <- data.frame(
  number = numeric(),
  color = character(),
  size = numeric()
)

# retrieve module sizes
for(i in 1:numMods){
  number <- head(resultsTable[resultsTable$number == i,3],1)
  color <- head(resultsTable[resultsTable$number == i,2],1)
  size <- nrow(resultsTable[resultsTable$number == i,])
  moduleData <- cbind(number, color, size)
  modSizes <-rbind(modSizes, moduleData)
}


# set color and numner tags for NAs
resultsTable$color <- resultsTable$color %>% replace_na('None')
resultsTable$number <- resultsTable$number %>% replace_na('0')

# remove NAs
resultsTable <- na.omit(resultsTable)

# update the number of modules
numMods <- numMods + 1

# update the table of module colors
colorTable[nrow(colorTable) + 1,] <- c("None","0")

# read in GO mappings
GOmaps <- readMappings(file = "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/email/geneToGO_tagged_map.txt")

# create named list of all genes (gene universe) and p-values. The gene universe is set to be
# the list of all genes contained in the gene2GO list of annotated genes.
list_genes <- as.numeric(resultsTable$number)
list_genes <- setNames(list_genes, resultsTable$gene)
list_genes_filtered <- list_genes[names(list_genes) %in% names(GOmaps)]

## GO enrichment
#Create function to return list of interesting DE genes (0 == not significant, 1 == significant)
get_interesting_DE_genes <- function(geneUniverse){
  interesting_DE_genes <- rep(0, length(geneUniverse))
  for(i in 1:length(geneUniverse)){
    if(geneUniverse[i] == "salmon4"){
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

# retrieve geneIDs associated with all GO terms
# https://support.bioconductor.org/p/29775/
allGO_salmon4 = genesInTerm(BP_GO_data)

# import DEGs
treatmentTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Genotypes/glmQLF_2WayANOVA_treatment_topTags_LFC1.2.csv", row.names="gene")

# create function to return list of interesting DE genes (0 == not significant, 1 == significant)
get_interesting_DE_genes_treatment <- function(geneUniverse){
  interesting_DE_genes <- rep(0, length(geneUniverse))
  for(i in 1:length(geneUniverse)){
    if(geneUniverse[i] < 0.05){
      interesting_DE_genes[i] = 1
    }
  }
  interesting_DE_genes <- setNames(interesting_DE_genes, names(geneUniverse))
  return(interesting_DE_genes)
}

# create named list of all genes (gene universe) and p-values
# the gene universe is set to be the list of all genes contained in the gene2GO list of annotated genes
list_genes_treatment <- as.numeric(treatmentTable$FDR)
list_genes_treatment <- setNames(list_genes_treatment, rownames(treatmentTable))
list_genes_treatment_filtered <- list_genes_treatment[names(list_genes_treatment) %in% names(GOmaps)]

# create topGOdata objects for enrichment analysis
BP_GO_data_treatment <- new('topGOdata', ontology = 'BP', allGenes = list_genes_treatment_filtered, 
                            geneSel = get_interesting_DE_genes_treatment, nodeSize = 10, annot = annFUN.gene2GO, 
                            gene2GO = GOmaps)

# retrieve geneIDs associated with all GO terms
# https://support.bioconductor.org/p/29775/
allGO_treatment = genesInTerm(BP_GO_data_treatment)


# setup plotting data
# We also searched for photoreactive repair (GO:0000719) and deoxyribodipyrimidine photo-lyase activity (GO:0003904) GO terms.
MMR <- unlist(data.frame(unlist(allGO_treatment["GO:0006298"])), use.names = FALSE)
cellOxidative <- unlist(data.frame(unlist(allGO_treatment["GO:0034599"])), use.names = FALSE)
#PR <- unlist(data.frame(unlist(allGO_treatment["GO:0000719"])), use.names = FALSE)
#PDR <- unlist(data.frame(unlist(allGO_treatment["GO:0006290"])), use.names = FALSE)
PRR <- unlist(data.frame(unlist(allGO_treatment["GO:0006301"])), use.names = FALSE)
salmon4ModulesGO <- c(resultsTable[resultsTable$color == "salmon4",]$gene)
names(salmon4ModulesGO) <- rep("salmon4", length(salmon4ModulesGO))
lightyellowModulesGO <- c(resultsTable[resultsTable$color == "lightyellow",]$gene)
names(lightyellowModulesGO) <- rep("lightyellow", length(lightyellowModulesGO))

# Heatmap of GO:0006298 DGE
# subset counts table by DE gene set
DGESubset_MMR <- tagsTblANOVATreatment[rownames(tagsTblANOVATreatment) %in% MMR,]
annotationsList_MMR <- subset(annotations,
                                    grepl(
                                      paste0(rownames(DGESubset_MMR), collapse = "|"),
                                      rownames(annotations),
                                      ignore.case = TRUE
                                    )
)
logcountsSubset_MMR <- subset(logcounts,
                                                  grepl(
                                                    paste0(rownames(DGESubset_MMR), collapse = "|"),
                                                    rownames(logcounts),
                                                    ignore.case = TRUE
                                                  )
)
#jpeg("glmQLF_MMR_heatmap.jpg")
heatmap(logcountsSubset_MMR, main= "Heatmap of GO:0006298 DGE", margins = c(8, 8), labRow = annotationsList_MMR$annotation)
#dev.off()

# Heatmap of GO:0034599 DGE
# subset counts table by DE gene set
DGESubset_cellOxidative <- tagsTblANOVATreatment[rownames(tagsTblANOVATreatment) %in% cellOxidative,]
annotationsList_cellOxidative <- subset(annotations,
                              grepl(
                                paste0(rownames(DGESubset_cellOxidative), collapse = "|"),
                                rownames(annotations),
                                ignore.case = TRUE
                              )
)
logcountsSubset_cellOxidative <- subset(logcounts,
                                        grepl(
                                          paste0(rownames(DGESubset_cellOxidative), collapse = "|"),
                                          rownames(logcounts),
                                          ignore.case = TRUE
                                        )
)
#jpeg("glmQLF_cellOxidative_heatmap.jpg")
heatmap(logcountsSubset_cellOxidative, main= "Heatmap of GO:0034599 DGE", margins = c(8, 8), labRow = annotationsList_cellOxidative$annotation)
#dev.off()

# Heatmap of GO:0000719 DGE
# subset counts table by DE gene set
#DGESubset_PR <- tagsTblANOVATreatment[rownames(tagsTblANOVATreatment) %in% PR,]
#logcountsSubset_PR <- subset(logcounts,
#                                        grepl(
#                                          paste0(rownames(DGESubset_PR), collapse = "|"),
#                                          rownames(logcounts),
#                                          ignore.case = TRUE
#                                        )
#)
##jpeg("glmQLF_PR_heatmap.jpg")
#heatmap(logcountsSubset_PR, main= "Heatmap of GO:0000719 DGE", margins = c(8, 1))
##dev.off()

# Heatmap of GO:0006290 DGE
# subset counts table by DE gene set
#DGESubset_PDR <- tagsTblANOVATreatment[rownames(tagsTblANOVATreatment) %in% PDR,]
#logcountsSubset_PDR <- subset(logcounts,
#                                        grepl(
#                                          paste0(rownames(DGESubset_PDR), collapse = "|"),
#                                          rownames(logcounts),
#                                          ignore.case = TRUE
#                                        )
#)
##jpeg("glmQLF_PDR_heatmap.jpg")
#heatmap(logcountsSubset_PDR, main= "Heatmap of GO:0006290 DGE", margins = c(8, 1))
##dev.off()

# Heatmap of GO:0006301 DGE
# subset counts table by DE gene set
DGESubset_PRR <- tagsTblANOVATreatment[rownames(tagsTblANOVATreatment) %in% PRR,]
annotationsList_PRR <- subset(annotations,
                                        grepl(
                                          paste0(rownames(DGESubset_PRR), collapse = "|"),
                                          rownames(annotations),
                                          ignore.case = TRUE
                                        )
)
logcountsSubset_PRR <- subset(logcounts,
                                        grepl(
                                          paste0(rownames(DGESubset_PRR), collapse = "|"),
                                          rownames(logcounts),
                                          ignore.case = TRUE
                                        )
)
#jpeg("glmQLF_PRR_heatmap.jpg")
heatmap(logcountsSubset_PRR, main= "Heatmap of GO:0006301 DGE", margins = c(8, 8), labRow = annotationsList_PRR$annotation)
#dev.off()

# Heatmap of Salmon4 Module Genes DGE
#resultsTable_salmon4 <- resultsTable[grepl("salmon4", resultsTable$color),]
# subset counts table by DE gene set
DGESubset_salmon4 <- tagsTblANOVATreatment[rownames(tagsTblANOVATreatment) %in% resultsTable_salmon4$gene,]
annotationsList_salmon4 <- subset(annotations,
                                  grepl(
                                    paste0(rownames(DGESubset_salmon4), collapse = "|"),
                                    rownames(annotations),
                                    ignore.case = TRUE
                                  )
)
logcountsSubset_salmon4 <- subset(logcounts,
                                  grepl(
                                    paste0(rownames(DGESubset_salmon4), collapse = "|"),
                                    rownames(logcounts),
                                    ignore.case = TRUE
                                  )
)
#jpeg("glmQLF_salmon4_heatmap.jpg")
heatmap(logcountsSubset_salmon4, main= "Heatmap of Salmon4 Module Genes DGE", margins = c(8, 8), labRow = annotationsList_salmon4$annotation)
#dev.off()

# Heatmap of Lightyellow Module Genes DGE
resultsTable_lightyellow <- resultsTable[grepl("lightyellow", resultsTable$color),]
# subset counts table by DE gene set
DGESubset_lightyellow <- tagsTblANOVATreatment[rownames(tagsTblANOVATreatment) %in% resultsTable_lightyellow$gene,]
annotationsList_lightyellow <- subset(annotations,
                                            grepl(
                                              paste0(rownames(DGESubset_lightyellow), collapse = "|"),
                                              rownames(annotations),
                                              ignore.case = TRUE
                                            )
)
logcountsSubset_lightyellow <- subset(logcounts,
                                      grepl(
                                        paste0(rownames(DGESubset_lightyellow), collapse = "|"),
                                        rownames(logcounts),
                                        ignore.case = TRUE
                                      )
)
#jpeg("glmQLF_lightyellow_heatmap.jpg")
heatmap(logcountsSubset_lightyellow, main= "Heatmap of Lightyellow Module Genes DGE", margins = c(8,8), labRow = annotationsList_lightyellow$annotation)
#dev.off()

# Heatmap of salmon4 & Treatment Effect DGE
resultsTable_salmon4 <- resultsTable[grepl("salmon4", resultsTable$color),]
# subset counts table by DE gene set
DGESubset_treatment_salmon4 <- tagsTblANOVATreatment.filtered[rownames(tagsTblANOVATreatment.filtered) %in% resultsTable_salmon4$gene,]
annotationsList_treatment_salmon4 <- subset(annotations,
                              grepl(
                                paste0(rownames(DGESubset_treatment_salmon4), collapse = "|"),
                                rownames(annotations),
                                ignore.case = TRUE
                              )
)
logcountsSubset_treatment_salmon4 <- subset(logcounts,
                                    grepl(
                                      paste0(rownames(DGESubset_treatment_salmon4), collapse = "|"),
                                      rownames(logcounts),
                                      ignore.case = TRUE
                                    )
)
#jpeg("glmQLF_treatment_salmon4_heatmap.jpg")
heatmap(logcountsSubset_treatment_salmon4, main= "Heatmap of salmon4 & Treatment Effect DGE", margins = c(8, 8), labRow = annotationsList_treatment_salmon4$annotation)
#dev.off()

# Heatmap of GO:0034599 DGE in salmon4 & Treatment Effect
# subset counts table by DE gene set
DGESubset_cellOxidative_treatment_salmon4 <- DGESubset_treatment_salmon4[rownames(DGESubset_treatment_salmon4) %in% cellOxidative,]
logcountsSubset_cellOxidative_treatment_salmon4 <- subset(logcounts,
                                            grepl(
                                              paste0(rownames(DGESubset_cellOxidative_treatment_salmon4), collapse = "|"),
                                              rownames(logcounts),
                                              ignore.case = TRUE
                                            )
)
##jpeg("glmQLF_cellOxidative_treatment_salmon4_heatmap.jpg")
heatmap(logcountsSubset_cellOxidative_treatment_salmon4, main= "Heatmap of GO:0034599 DGE in the salmon4 & Treatment Effect", margins = c(8, 1))
##dev.off()

# Heatmap of GO:0034599 DGE in the Treatment Effect
# subset counts table by DE gene set
DGESubset_cellOxidative_treatment <- tagsTblANOVATreatment.filtered[rownames(tagsTblANOVATreatment.filtered) %in% cellOxidative,]
logcountsSubset_cellOxidative_treatment <- subset(logcounts,
                                                          grepl(
                                                            paste0(rownames(DGESubset_cellOxidative_treatment), collapse = "|"),
                                                            rownames(logcounts),
                                                            ignore.case = TRUE
                                                          )
)
##jpeg("glmQLF_cellOxidative_treatment_heatmap.jpg")
heatmap(logcountsSubset_cellOxidative_treatment, main= "Heatmap of GO:0034599 DGE in the Treatment Effect", margins = c(8, 1))
##dev.off()

# Heatmap of GO:0034599 DGE in salmon4
# subset counts table by DE gene set
DGESubset_cellOxidative_salmon4 <- resultsTable_salmon4[resultsTable_salmon4$gene %in% cellOxidative,]
logcountsSubset_cellOxidative_salmon4 <- subset(logcounts,
                                                          grepl(
                                                            paste0(DGESubset_cellOxidative_salmon4$gene, collapse = "|"),
                                                            rownames(logcounts),
                                                            ignore.case = TRUE
                                                          )
)
##jpeg("glmQLF_cellOxidative_salmon4_heatmap.jpg")
heatmap(logcountsSubset_cellOxidative_salmon4, main= "Heatmap of GO:0034599 DGE in salmon4", margins = c(8, 1))
##dev.off()

# Heatmap of GO:0006298 DGE in salmon4 & Treatment Effect
# subset counts table by DE gene set
DGESubset_MMR_salmon4 <- resultsTable_salmon4[resultsTable_salmon4$gene %in% MMR,]
logcountsSubset_MMR_salmon4 <- subset(logcounts,
                                                grepl(
                                                  paste0(DGESubset_MMR_salmon4$gene, collapse = "|"),
                                                  rownames(logcounts),
                                                  ignore.case = TRUE
                                                )
)
##jpeg("glmQLF_MMR_salmon4_heatmap.jpg")
heatmap(logcountsSubset_MMR_salmon4, main= "Heatmap of GO:0006298 DGE in salmon4", margins = c(8, 1))
##dev.off()

