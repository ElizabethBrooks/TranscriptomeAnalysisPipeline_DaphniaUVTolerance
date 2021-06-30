#Set working directory
#workingDir = args[1];
workingDir="/home/mae/Documents/RNASeq_Workshop_ND/WGCNA_PA42_v4.1"
setwd(workingDir); 

# Load the WGCNA package
library(WGCNA)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

# Load the expression and trait data saved in the first part
lnames = load(file = "PA42_v4.1_entrezSubset_dataInput.RData");
#The variable lnames contains the names of loaded variables.
#lnames

# Load network data saved in the second part.
lnames = load(file = "PA42_v4.1_networkConstruction_auto_threshold7.RData");
#lnames


annot = read.csv(file = "geneAnnotations_entrezSubset.csv");
#dim(annot)
#names(annot)
probes = names(datExpr)
probes2annot = match(probes, annot[,1])
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.


# Get the corresponding entrez IDs
allLLIDs = annot$entrezID[probes2annot];
# $ Choose interesting modules
intModules = c("brown", "yellow", "red", "lightcyan")
for (module in intModules)
{
  # Select module probes
  modGenes = (moduleColors==module)
  # Get their entrez ID codes
  modLLIDs = allLLIDs[modGenes];
  # Write them into a file
  fileName = paste("LocusLinkIDs-", module, ".txt", sep="");
  write.table(as.data.frame(modLLIDs), file = fileName,
              row.names = FALSE, col.names = FALSE)
}
# As background in the enrichment analysis, we will use all probes in the analysis.
fileName = paste("LocusLinkIDs-all.txt", sep="");
write.table(as.data.frame(allLLIDs), file = fileName,
            row.names = FALSE, col.names = FALSE)


#Perform Go enrichment analysis
GOenr = GOenrichmentAnalysis(moduleColors, allLLIDs, organism = "fly", nBestP = 10);


tab = GOenr$bestPTerms[[4]]$enrichment


names(tab)


write.table(tab, file = "GOEnrichmentTable.csv", sep = ",", quote = TRUE, row.names = FALSE)


keepCols = c(1, 2, 5, 6, 7, 12, 13);
screenTab = tab[, keepCols];
# Round the numeric columns to 2 decimal places:
numCols = c(3, 4);
screenTab[, numCols] = signif(apply(screenTab[, numCols], 2, as.numeric), 2)
# Truncate the the term name to at most 40 characters
screenTab[, 7] = substring(screenTab[, 7], 1, 40)
# Shorten the column names:
colnames(screenTab) = c("module", "size", "p-val", "Bonf", "nInTerm", "ont", "term name");
rownames(screenTab) = NULL;
# Set the width of R's output. The reader should play with this number to obtain satisfactory output.
options(width=95)
# Finally, display the enrichment table:
screenTab

