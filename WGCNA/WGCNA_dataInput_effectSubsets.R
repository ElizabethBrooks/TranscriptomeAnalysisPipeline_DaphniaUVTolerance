#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("WGCNA")

#Load the WGCNA and edgeR packages
library(WGCNA);
library("edgeR")

#The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

#Set working directory
#workingDir = args[1];
workingDir="/home/mae/Documents/RNASeq_Workshop_ND/WGCNA_PA42_v4.1/effectSubsets"
setwd(workingDir); 


#Import gene count data
#countsTable <- read.csv(file=args[2], row.names="gene")[ ,args[3]:args[4]]
countsTable <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/geneCounts_mergedHisat2_PA42_v4.1/cleaned.csv", row.names="gene", header=TRUE)[ ,1:24]

#Import grouping factor
targets <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_WGCNA_Olympics.csv", row.names="sample")

#Setup a design matrix
group <- factor(paste(targets$treatment,targets$tolerance,sep="."))
#Create DGE list object
list <- DGEList(counts=countsTable,group=group)
colnames(list) <- targets$sample
#Retain genes only if it is expressed at a minimum level
keep <- filterByExpr(list)
list <- list[keep, , keep.lib.sizes=FALSE]
#Use TMM normalization to eliminate composition biases between libraries
list <- calcNormFactors(list)
#Write normalized counts to file
countsTableNorm <- cpm(list, normalized.lib.sizes=TRUE)

#Import DEGs
geneCountsInter <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/geneCounts_mergedHisat2_PA42_v4.1/glmQLFAnalysis_FDR0.10/glmQLF_2WayANOVA_interaction_topTags_filtered.csv")
geneCountsTreat <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/geneCounts_mergedHisat2_PA42_v4.1/glmQLFAnalysis_FDR0.10/glmQLF_2WayANOVA_UVvsVIS_topTags_filtered.csv")
geneCountsTol <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/geneCounts_mergedHisat2_PA42_v4.1/glmQLFAnalysis_FDR0.10/glmQLF_2WayANOVA_TvsN_topTags_filtered.csv")
SETInterIn <- geneCountsInter[,1]
SETTreatIn <- geneCountsTreat[,1]
SETTolIn <- geneCountsTol[,1]

#Generate a subset of the gene counts
normListInter <- countsTableNorm[rownames(countsTableNorm) %in% SETInterIn,]
normListTreat <- countsTableNorm[rownames(countsTableNorm) %in% SETTreatIn,]
normListTol <- countsTableNorm[rownames(countsTableNorm) %in% SETTolIn,]

#Import annotation data
annotIn = read.csv(file = "/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/trinotate_annotation_report_PA42_v4.1_transcripts.csv", sep="\t");
annotUniprot = read.csv(file = "/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/trinotate_annotation_report_PA42_v4.1_transcripts_uniprot.csv", sep="\t");
annot <- cbind(annotIn,annotUniprot)

#Retrieve uniprot to enztrez ID mapping file
annotMap = read.csv(file = "/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/trinotate_annotation_report_PA42_v4.1_transcripts_uniprotEntrezMap.csv", sep="\t");

#Check if each uniprot ID has an entrez ID mapping
for(var in 1:nrow(annot))
{
  annot$entrezID[var] <- ifelse(annot[var,18] %in% annotMap[,1], annotMap[grep(annot[var,18], annotMap[,1]),2], ".")
}

#Retrieve the subset of genes with entrez IDs
annotTable <- subset(annot, entrezID!=".")

# The following is the % of genes in the subset without entrez annotation:
print("% interaction effect DEGs without entrez IDs: ")
geneAnnot = match(rownames(normListInter), annotTable[,1])
sum(is.na(geneAnnot))/nrow(normListInter)*100
print("% treatment effect DEGs without entrez IDs: ")
geneAnnot = match(rownames(normListTreat), annotTable[,1])
sum(is.na(geneAnnot))/nrow(normListTreat)*100
print("% tolerance effect DEGs without entrez IDs: ")
geneAnnot = match(rownames(normListTol), annotTable[,1])
sum(is.na(geneAnnot))/nrow(normListTol)*100


#Transpose the input interaction data
datExpr0Inter = as.data.frame(t(normListInter));
names(datExpr0Inter) = rownames(normListInter);
rownames(datExpr0Inter) = names(countsTable);
#Transpose the input treatment data
datExpr0Treat = as.data.frame(t(normListTreat));
names(datExpr0Treat) = rownames(normListTreat);
rownames(datExpr0Treat) = names(countsTable);
#Transpose the input tolerance data
datExpr0Tol = as.data.frame(t(normListTol));
names(datExpr0Tol) = rownames(normListTol);
rownames(datExpr0Tol) = names(countsTable);

#Check the genes across all samples
gsgInter = goodSamplesGenes(datExpr0Inter, verbose = 3);
gsgInter$allOK
gsgTreat = goodSamplesGenes(datExpr0Treat, verbose = 3);
gsgTreat$allOK
gsgTol = goodSamplesGenes(datExpr0Tol, verbose = 3);
gsgTol$allOK


#Remove genes that do not pass the check
#if (!gsg$allOK)
#{
  # Optionally, print the gene and sample names that were removed:
#  if (sum(!gsg$goodGenes)>0) 
#    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
#  if (sum(!gsg$goodSamples)>0) 
#    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
#  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
#}


#Cluster interaction samples to identify outliers
sampleTreeInter = hclust(dist(datExpr0Inter), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
#sizeGrWindow(12,9)
#pdf(file = "sampleClustering.pdf", width = 12, height = 9);
jpeg("sampleClusteringInter_cutLine.jpg", width = 480, height = 480)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTreeInter, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
dev.off()
#Cluster treatment samples to identify outliers
sampleTreeTreat = hclust(dist(datExpr0Treat), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
#sizeGrWindow(12,9)
#pdf(file = "sampleClustering.pdf", width = 12, height = 9);
jpeg("sampleClusteringTreat_cutLine.jpg", width = 480, height = 480)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTreeTreat, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
dev.off()
#Cluster tolerance samples to identify outliers
sampleTreeTol = hclust(dist(datExpr0Tol), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
#sizeGrWindow(12,9)
#pdf(file = "sampleClustering.pdf", width = 12, height = 9);
jpeg("sampleClusteringTol_cutLine.jpg", width = 480, height = 480)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTreeTol, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
dev.off()


# Determine cluster under the line
#clust = cutreeStatic(sampleTree, cutHeight = 8500, minSize = 10)
#table(clust)
# clust 1 contains the samples we want to keep.
#keepSamples = (clust==1)
#datExpr = datExpr0[keepSamples, ]
#nGenes = ncol(datExpr)
#nSamples = nrow(datExpr)

datExprInter = datExpr0Inter
datExprTreat = datExpr0Treat
datExprTol = datExpr0Tol


# remove columns that hold information we do not need.
#allTraits = read.csv(args[5]);
allTraits = read.csv("/home/mae/Documents/RNASeq_Workshop_ND/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_WGCNA_Olympics.csv");
#dim(allTraits)
#names(allTraits)


# Form a interaction data frame analogous to expression data that will hold the clinical traits.
inSamplesInter = rownames(datExprInter);
traitRowsInter = match(inSamplesInter, allTraits$sample);
datTraitsInter = allTraits[traitRowsInter, -1];
rownames(datTraitsInter) = allTraits[traitRowsInter, 1];
# Form a treatment data frame analogous to expression data that will hold the clinical traits.
inSamplesTreat = rownames(datExprTreat);
traitRowsTreat = match(inSamplesTreat, allTraits$sample);
datTraitsTreat = allTraits[traitRowsTreat, -1];
rownames(datTraitsTreat) = allTraits[traitRowsTreat, 1];
# Form a tolerance data frame analogous to expression data that will hold the clinical traits.
inSamplesTol = rownames(datExprTol);
traitRowsTol = match(inSamplesTol, allTraits$sample);
datTraitsTol = allTraits[traitRowsTol, -1];
rownames(datTraitsTol) = allTraits[traitRowsTol, 1];

collectGarbage();


# Re-cluster samples
#sampleTree2 = hclust(dist(datExpr), method = "average")
#sizeGrWindow(12,9)
#pdf(file = "sampleReClustering.pdf", width = 12, height = 9);
#jpeg("sampleReClustering.jpg", width = 480, height = 480)
#par(cex = 0.6);
#par(mar = c(0,4,2,0))
#plot(sampleTree2, main = "Sample clustering", sub="", xlab="", cex.lab = 1.5, 
#     cex.axis = 1.5, cex.main = 2)
#dev.off()


#Save the formatted data for input to the next stage of analysis
save(datExprInter, datTraitsInter, file = "PA42_v4.1_dataInputInter.RData")
save(datExprTreat, datTraitsTreat, file = "PA42_v4.1_dataInputTreat.RData")
save(datExprTol, datTraitsTol, file = "PA42_v4.1_dataInputTol.RData")
