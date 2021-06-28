#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("WGCNA")

#Load the WGCNA and edgeR packages
library(WGCNA);
library("edgeR")

#The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

#Set working directory
workingDir = args[1];
setwd(workingDir); 

#Import gene count data
countsTable <- read.csv(file=args[2], row.names="gene")[ ,args[3]:args[4]]
#head(countsTable)
#Import grouping factor
targets <- read.csv(file=args[5], row.names="sample")

#Setup a design matrix
group <- factor(paste(targets$treatment,targets$genotype,sep="."))
#cbind(targets,Group=group)
#Create DGE list object
list <- DGEList(counts=countsTable,group=group)
colnames(list) <- targets$sample

#Retain genes only if it is expressed at a minimum level
keep <- filterByExpr(list)
summary(keep)
list <- list[keep, , keep.lib.sizes=FALSE]

#Use TMM normalization to eliminate composition biases
# between libraries
list <- calcNormFactors(list)
#list$samples
#Write normalized counts to file
normList <- cpm(list, normalized.lib.sizes=TRUE)

#Take a quick look at what is in the data set:
dim(normList);
names(normList);

#Transpose the input data
datExpr0 = as.data.frame(t(normList));
names(datExpr0) = normList[,1];
rownames(datExpr0) = names(countsTable);

#Check the genes across all samples
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK


#Remove genes that do not pass the check
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}


#Cluster samples to identify outliers
sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
pdf(file = "sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)


# Plot a line to show the cut
jpeg("sampleClustering_cutLine.jpg")
abline(h = 8500, col = "red");
dev.off()

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 8500, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


# remove columns that hold information we do not need.
allTraits = read.csv(args[5]);
dim(allTraits)
names(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits.

femaleSamples = rownames(datExpr);
traitRows = match(femaleSamples, allTraits$sample);
datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 1];

collectGarbage();


# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
sizeGrWindow(12,9)
pdf(file = "sampleReClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree2, main = "Sample clustering", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)


#Save the formatted data for input to the next stage of analysis
save(datExpr, datTraits, file = "PA42_v4.1_WGCNA_inputData.RData")


