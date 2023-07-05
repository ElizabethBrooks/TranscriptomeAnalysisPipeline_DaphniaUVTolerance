#if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
#BiocManager::install("gridExtra")

#Load libraries
library(ggplot2)
library(gridExtra)
library(rcartocolor)

# Plotting Palettes
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
# https://github.com/Nowosad/rcartocolor
plotColors <- carto_pal(12, "Safe")
plotColorSubset <- c(plotColors[5], plotColors[6])

# turn off scientific notation
options(scipen = 999)

#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)

# retrieve working directory
#workingDir <- args[1]
workingDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Tolerance/GOAnalysis_OLYM_30"

# set working directory
setwd(workingDir)

# retrieve subset tag
#set <- args[2]
set <- "OLYM"

# set the minimum module size
#minModSize <- args[3]
minModSize <- "30"

# retrieve WGCNA directory
#inDir <- args[4]
inDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Tolerance"

# set the full subset tag name
tag <- paste(set, minModSize, sep="_")

# Load network data saved in the second part
importFile <- paste(tag, "networkConstruction-stepByStep.RData", sep="-")
importFile <- paste(inDir, importFile, sep="/")
lnames2 = load(file = importFile)

# create list of module colors mapped to numbers
num_mods <- length(unique(moduleColors))
color_list <- unique(moduleColors)

# create data frame to hold plotting data for each module
module_BP_results = data.frame(matrix(ncol = 4, nrow = 0))
colnames(module_BP_results) <- c("Color","Term","Significant","weightFisher")
module_MF_results = data.frame(matrix(ncol = 4, nrow = 0))
colnames(module_MF_results) <- c("Color","Term","Significant","weightFisher")
module_CC_results = data.frame(matrix(ncol = 4, nrow = 0))
colnames(module_CC_results) <- c("Color","Term","Significant","weightFisher")

# loop through each module color
for(j in 1:num_mods){
    # read in data on significant GO term (BP, MF, and CC) 
    BP_GO_sig <- read.csv(paste(color_list[j], 'BP_sigGO_terms.csv', sep="_"), row.names = NULL)[1:5, c("Term","Significant","weightFisher")]
    MF_GO_sig <- read.csv(paste(color_list[j], 'MF_sigGO_terms.csv', sep="_"), row.names = NULL)[1:5, c("Term","Significant","weightFisher")]
    CC_GO_sig <- read.csv(paste(color_list[j], 'CC_sigGO_terms.csv', sep="_"), row.names = NULL)[1:5, c("Term","Significant","weightFisher")]

    # add a column labeling the color to each GO term set
    BP_GO_sig <- cbind("Color" = color_list[j], BP_GO_sig)
    MF_GO_sig <- cbind("Color" = color_list[j], MF_GO_sig)
    CC_GO_sig <- cbind("Color" = color_list[j], CC_GO_sig)

    # create tables of all GO terms for each level
    #list_all_BP_GO_included <- unique(c(BP_GO_sig$GO.ID))
    #list_all_MF_GO_included <- unique(c(MF_GO_sig$GO.ID))
    #list_all_CC_GO_included <- unique(c(CC_GO_sig$GO.ID))

    # combine all tables
    module_BP_results <- rbind(module_BP_results, BP_GO_sig)
    module_MF_results <- rbind(module_MF_results, MF_GO_sig)
    module_CC_results <- rbind(module_CC_results, CC_GO_sig)
    
    # add effect tags
    plot_table$effect <- plot_table
}

# add GO level tags
module_BP_results <- cbind('Level' = 'BP', module_BP_results)
module_MF_results <- cbind('Level' = 'MF', module_MF_results)
module_CC_results <- cbind('Level' = 'CC', module_CC_results)

# combine GO level results for plotting
plot_table <- rbind(module_BP_results, module_MF_results, module_CC_results)

# remove NAs
plot_table <- na.omit(plot_table)

# remove < signs
plot_table$weightFisher <- gsub("<","",as.character(plot_table$weightFisher))

# import positively selected gene set p-values
positiveTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/selectionTests/fisherTest_positiveSelection_modules.csv", row.names="color")

# loop through each module color
plot_table$enrichment <- "NA"
plot_table$depletion <- "NA"
for(k in 1:num_mods){
  # add enrichment and depletion values
  plot_table[plot_table$Color == row.names(positiveTable)[k],]$enrichment <- positiveTable$enrichment[k]
  plot_table[plot_table$Color == row.names(positiveTable)[k],]$depletion <- positiveTable$depletion[k]
}

# add signifigance tags
plot_table$selection <- "None"
plot_table[plot_table$enrichment <= 0.05,]$selection <- "Enriched"
plot_table[plot_table$depletion <= 0.05,]$selection <- "Depleted"

# import effect ANOVA p-values
effectTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Tolerance/ANOVA_OLYM_30/aov_summary_pValues.csv", row.names="module")

# clean up row names
row.names(effectTable) <- gsub("ME","",as.character(row.names(effectTable)))

# loop through each module color
plot_table$treatment <- "NA"
plot_table$tolerance <- "NA"
plot_table$interaction <- "NA"
for(k in 1:num_mods){
  # add enrichment and depletion values
  plot_table[plot_table$Color == row.names(effectTable)[k],]$treatment <- effectTable$treatment[k]
  plot_table[plot_table$Color == row.names(effectTable)[k],]$tolerance <- effectTable$tolerance[k]
  plot_table[plot_table$Color == row.names(effectTable)[k],]$interaction <- effectTable$interaction[k]
}

# add signifigance tags
plot_table$effect <- "None"
plot_table[plot_table$treatment <= 0.05,]$effect <- "Treatment"
plot_table[plot_table$tolerance <= 0.05,]$effect <- "Tolerance"
plot_table[plot_table$interaction <= 0.05,]$effect <- "Interaction"

# create dot plot of significant GO terms
x_axis_order <- factor(plot_table$Color, levels = color_list)
facetLevel <- factor(plot_table$Level, levels = c('BP', 'MF', 'CC'))
facetSelection <- factor(plot_table$selection, levels = c('None', 'Enriched', 'Depleted'))
facetEffect <- factor(plot_table$effect, levels = c('None', 'Treatment', 'Tolerance', 'Interaction'))

dotplot <- ggplot(data = plot_table, aes(x = x_axis_order, y = Term, size = Significant, color = as.numeric(weightFisher))) + 
  facet_grid(rows = facet, space = 'free_y', scales = 'free') +
  geom_point() +
  scale_color_gradientn(colors = heat.colors(10), limits=c(0, 0.05)) + 
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90)) +
  xlab('Color') +
  ylab('GO Term') + 
  scale_x_discrete(labels=c("navajowhite2"=expression(bold("navajowhite2")), "bisque4"=expression(bold("bisque4")), parse=TRUE)) +
  scale_x_discrete(labels=c("navajowhite2"=expression(italics("navajowhite2")), "bisque4"=expression(italics("bisque4")), "magenta"=expression(italics("magenta")), parse=TRUE)) +
  labs(color = 'FDR Adjusted p-Value', size = 'Gene rank')

# view plot
dotplot

# save the plot to a PDF file
ggsave('dotplotModule_sigGO.pdf', plot = dotplot, device = 'pdf')


