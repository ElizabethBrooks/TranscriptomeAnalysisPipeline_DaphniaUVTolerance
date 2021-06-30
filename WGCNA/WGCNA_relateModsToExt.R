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


#Setup inputs
geneIDList = probes
uniprotIDList = annot$uniprotID[probes2annot]
entrezIDList = annot$entrezID[probes2annot]


# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


#sizeGrWindow(10,6)
jpeg("moduleTraitRelationships.jpg", width = 960, height = 960)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
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
tolerance = as.data.frame(datTraits$tolerance);
names(tolerance) = "tolerance"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, tolerance, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(tolerance), sep="");
names(GSPvalue) = paste("p.GS.", names(tolerance), sep="");

#Plot MM vs GS for tolerance modules
module = "red"
column = match(module, modNames);
moduleGenes = moduleColors==module;

#sizeGrWindow(7, 7);
jpeg("MMvsGS_tolerance_red.jpg", width = 960, height = 480)
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for tolerance",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2)
dev.off()

module = "lightcyan"
column = match(module, modNames);
moduleGenes = moduleColors==module;

#sizeGrWindow(7, 7);
jpeg("MMvsGS_tolerance_lightcyan.jpg", width = 960, height = 480)
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for tolerance",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2)
dev.off()

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

write.csv(geneInfo, file = "geneInfo_tolerance.csv")


# Define variable treatment containing the treatment column of datTrait
treatment = as.data.frame(datTraits$treatment);
names(treatment) = "treatment"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, treatment, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(treatment), sep="");
names(GSPvalue) = paste("p.GS.", names(treatment), sep="");


#Plot MM vs GS for treatment modules
module = "yellow"
column = match(module, modNames);
moduleGenes = moduleColors==module;

#sizeGrWindow(7, 7);
jpeg("MMvsGS_treatment_yellow.jpg", width = 960, height = 480)
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for treatment",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2)
dev.off()

module = "brown"
column = match(module, modNames);
moduleGenes = moduleColors==module;

#sizeGrWindow(7, 7);
jpeg("MMvsGS_treatment_brown.jpg", width = 960, height = 480)
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for treatment",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2)
dev.off()

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

write.csv(geneInfo, file = "geneInfo_treatment.csv")