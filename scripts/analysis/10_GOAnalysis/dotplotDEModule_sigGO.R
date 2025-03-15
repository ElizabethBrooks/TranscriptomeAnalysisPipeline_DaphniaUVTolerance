#if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
#BiocManager::install("gridExtra")

#Load libraries
library(ggplot2)
library(gridExtra)
library(rcartocolor)
library(dplyr)

# Plotting Palettes
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
# https://github.com/Nowosad/rcartocolor
plotColors <- carto_pal(12, "Safe")
plotColorSubset <- c(plotColors[5], plotColors[6])

# turn off scientific notation
options(scipen = 999)


# DE data
# retrieve working directory
workingDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Genotypes/GOAnalysis"

# set working directory
setwd(workingDir)

# read in data on significant GO terms (BP, MF, and CC) for each effect
treatment_BP_GO_terms <- read.csv('treatment_BP_sigGO_terms.csv')
#treatment_MF_GO_terms <- read.csv('treatment_MF_sigGO_terms.csv')
#treatment_CC_GO_terms <- read.csv('treatment_CC_sigGO_terms.csv')

## treatment effect
# filter for top 5 significant terms
treatment_BP_GO_sig <- treatment_BP_GO_terms#[1:5, ]
#treatment_MF_GO_sig <- treatment_MF_GO_terms#[1:5, ]
#treatment_CC_GO_sig <- treatment_CC_GO_terms#[1:5, ]

# add a column labeling the effect to each GO term set
treatment_BP_plot_table <- cbind("Set" = 'Treatment', treatment_BP_GO_sig)
#treatment_MF_plot_table <- cbind("Set" = 'Treatment', treatment_MF_GO_sig)
#treatment_CC_plot_table <- cbind("Set" = 'Treatment', treatment_CC_GO_sig)

# add GO level tags
treatment_BP_plot_table <- cbind('Level' = 'BP', treatment_BP_plot_table)
#treatment_MF_plot_table <- cbind('Level' = 'MF', treatment_MF_plot_table)
#treatment_CC_plot_table <- cbind('Level' = 'CC', treatment_CC_plot_table)

# plotting
#all_plot_table <- rbind(treatment_BP_plot_table, treatment_MF_plot_table, treatment_CC_plot_table)

# select columns
plot_tableDE <- select(treatment_BP_plot_table, Level, Set, GO.ID, Term, Significant, weightFisher)


# module data
# retrieve working directory
workingDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Genotypes/GOAnalysis_OLYM_30"

# set working directory
setwd(workingDir)

# retrieve subset tag
set <- "OLYM"

# set the minimum module size
minModSize <- "30"

# retrieve WGCNA directory
inDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Genotypes"

# set the full subset tag name
tag <- paste(set, minModSize, sep="_")

# Load network data saved in the second part
importFile <- paste(tag, "networkConstruction-stepByStep.RData", sep="-")
importFile <- paste(inDir, importFile, sep="/")
lnames2 = load(file = importFile)

# create list of module colors mapped to numbers
num_mods <- length(unique(moduleColors)) + 1
color_list <- c(unique(moduleColors), "None")

# create data frame to hold plotting data for each module
module_BP_results = data.frame(matrix(ncol = 4, nrow = 0))
colnames(module_BP_results) <- c("Set","Term","Significant","weightFisher")
#module_MF_results = data.frame(matrix(ncol = 4, nrow = 0))
#colnames(module_MF_results) <- c("Set","Term","Significant","weightFisher")
#module_CC_results = data.frame(matrix(ncol = 4, nrow = 0))
#colnames(module_CC_results) <- c("Set","Term","Significant","weightFisher")

# loop through each module color
for(j in 1:num_mods){
  # read in data on significant GO terms (BP, MF, and CC) 
  BP_GO_sig <- read.csv(paste(color_list[j], 'BP_sigGO_terms.csv', sep="_"), row.names = NULL)#[1:5,]
  #MF_GO_sig <- read.csv(paste(color_list[j], 'MF_sigGO_terms.csv', sep="_"), row.names = NULL)#[1:5,]
  #CC_GO_sig <- read.csv(paste(color_list[j], 'CC_sigGO_terms.csv', sep="_"), row.names = NULL)#[1:5,]
  
  # add a column labeling the color to each GO term set
  BP_GO_sig <- cbind("Set" = color_list[j], BP_GO_sig[,1], BP_GO_sig[,2], BP_GO_sig[,4], BP_GO_sig[,6])
  #MF_GO_sig <- cbind("Set" = color_list[j], MF_GO_sig[,1], MF_GO_sig[,2], MF_GO_sig[,4], MF_GO_sig[,6])
  #CC_GO_sig <- cbind("Set" = color_list[j], CC_GO_sig[,1], CC_GO_sig[,2], CC_GO_sig[,4], CC_GO_sig[,6])
  
  # combine all tables
  module_BP_results <- rbind(module_BP_results, BP_GO_sig)
  #module_MF_results <- rbind(module_MF_results, MF_GO_sig)
  #module_CC_results <- rbind(module_CC_results, CC_GO_sig)
}

# re-name columns
colnames(module_BP_results) <- c("Set","GO.ID","Term","Significant","weightFisher")
#colnames(module_MF_results) <- c("Set","GO.ID","Term","Significant","weightFisher")
#colnames(module_CC_results) <- c("Set","GO.ID","Term","Significant","weightFisher")

# add GO level tags
module_BP_results <- cbind('Level' = 'BP', module_BP_results)
#module_MF_results <- cbind('Level' = 'MF', module_MF_results)
#module_CC_results <- cbind('Level' = 'CC', module_CC_results)

# combine GO level results for plotting
#plotTableModule <- rbind(module_BP_results, module_MF_results, module_CC_results)

# remove NAs
plotTableModule <- na.omit(module_BP_results)

# remove < signs
plotTableModule$weightFisher <- gsub("<","",as.character(plotTableModule$weightFisher))

# select columns
plot_tableModule <- select(plotTableModule, Level, Set, GO.ID, Term, Significant, weightFisher)

# select the salmon4 and skyblue modules
plot_tableModule <- plot_tableModule[plot_tableModule$Set %in% c("salmon4","skyblue"),]

# combine tables
plot_table <- rbind(plot_tableDE, plot_tableModule)

# subset plot table
plot_tableSubset <- plot_table[plot_table$weightFisher < 0.026,]

# setup DE facet groups
x_axis_order <- factor(plot_tableSubset$Set, levels = c('Treatment', 'salmon4', 'skyblue'))
#facet <- factor(plot_tableSubset$Level, levels = c('BP', 'CC', 'MF'))

# create dot plot of significant GO terms
dotplot <- ggplot(data = plot_tableSubset, aes(x = x_axis_order, y = GO.ID, size = as.numeric(Significant), color = as.numeric(weightFisher))) + 
  #facet_grid(rows = facet, space = 'free_y', scales = 'free') +
  geom_point() +
  #scale_color_gradientn(colors = heat.colors(10), limits=c(0, 0.05)) + 
  scale_color_gradientn(colors = plotColorSubset) +
  theme_bw()+
  #theme(axis.text.x = element_text(angle = 90, color = plot_tableSubset$selection)) +
  theme(axis.text.x = element_text(angle = 90)) +
  #geom_text(position = position_dodge(width = 1), aes(x=effect, y=0)) +
  xlab('Color') +
  ylab('GO Term') + 
  labs(color = 'P-Value', size = 'Gene Rank')

# view plot
dotplot

# save the plot to a PDF file
#ggsave('dotplotDEModule_sigGO.pdf', plot = dotplot, device = 'pdf')
