#Set working directory
#workingDir = args[1];
workingDir="/home/mae/Documents/RNASeq_Workshop_ND/WGCNA_PA42_v4.1/effectSubsets"
setwd(workingDir); 

# Load the WGCNA package
library(WGCNA)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments. 
# See note above.
#enableWGCNAThreads()


# Load the data saved in the first part
lnames1 = load(file = "PA42_v4.1_dataInputInter.RData");

# Load network data saved in the second part.
lnames2 = load(file = "PA42_v4.1_networkConstructionInter_auto_threshold8_signed.RData");

#Import gene annotations
annot = read.csv(file = "geneAnnotations.csv");


# The following is the number of tolerance probes without annotation:
probes = names(datExprInter)
probes2annot = match(probes, annot[,1])
sum(is.na(probes2annot))
# Should return 0.


#Setup tolerance inputs
geneIDList = probes
uniprotIDList = annot$uniprotID[probes2annot]
entrezIDList = annot$entrezID[probes2annot]


# Define numbers of tolerance genes and samples
nGenes = ncol(datExprInter);
nSamples = nrow(datExprInter);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExprInter, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraitsInter, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


#sizeGrWindow(10,6)
jpeg("moduleTraitRelationshipsInter.jpg", width = 960, height = 960)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraitsInter),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


# Define variable tolerance containing the tolerance column of datTrait
tolerance = as.data.frame(datTraitsInter$tolerance);
names(tolerance) = "tolerance"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExprInter, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExprInter, tolerance, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(tolerance), sep="");
names(GSPvalue) = paste("p.GS.", names(tolerance), sep="");

# Create the starting data frame
geneInfo0 = data.frame(geneID = geneIDList,
                       uniprotID = uniprotIDList,
                       entrezID = entrezIDList,
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)

# Order modules by their significance for tolerance
modOrder = order(-abs(cor(MEs, tolerance, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.tolerance));
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file = "geneInfoInter_tolerance.csv")


# Define variable treatment containing the treatment column of datTrait
treatment = as.data.frame(datTraitsInter$treatment);
names(treatment) = "treatment"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExprInter, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExprInter, treatment, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(treatment), sep="");
names(GSPvalue) = paste("p.GS.", names(treatment), sep="");

# Create the starting data frame
geneInfo0 = data.frame(geneID = geneIDList,
                       uniprotID = uniprotIDList,
                       entrezID = entrezIDList,
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)


# Order modules by their significance for treatment
modOrder = order(-abs(cor(MEs, treatment, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.treatment));
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file = "geneInfoInter_treatment.csv")
