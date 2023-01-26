#if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
#BiocManager::install("gridExtra")

#Load libraries
library(ggplot2)
library(gridExtra)

# turn off scientific notation
options(scipen = 999)

#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)

# retrieve working directory
workingDir <- args[1]
#workingDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Tolerance/GOAnalysis_OLYM_30"

# set working directory
setwd(workingDir)

# retrieve subset tag
set <- args[2]
#set <- "OLYM"

# set the minimum module size
minModSize <- args[3]
#minModSize <- "30"

# retrieve WGCNA directory
inDir <- args[4]
#inDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Tolerance"

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

#Loop through each module color
for(j in 1:num_mods){
    # read in data on significant GO term (BP, MF, and CC) 
    BP_GO_sig <- read.csv(paste(color_list[j], 'BP_sigGO_terms.csv', sep="_"), row.names = NULL)[1, c("Term","Significant","weightFisher")]
    MF_GO_sig <- read.csv(paste(color_list[j], 'MF_sigGO_terms.csv', sep="_"), row.names = NULL)[1, c("Term","Significant","weightFisher")]
    CC_GO_sig <- read.csv(paste(color_list[j], 'CC_sigGO_terms.csv', sep="_"), row.names = NULL)[1, c("Term","Significant","weightFisher")]

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

# create dot plot of significant GO terms
x_axis_order <- factor(plot_table$Color, levels = color_list)
facet <- factor(plot_table$Level, levels = c('BP', 'MF', 'CC'))

dotplot <- ggplot(data = plot_table, aes(x = x_axis_order, y = Term, size = Significant, color = as.numeric(weightFisher))) + 
  facet_grid(rows = facet, space = 'free_y', scales = 'free') +
  geom_point() +
  scale_color_gradientn(colors = heat.colors(10), limits=c(0, 0.05)) + 
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90)) +
  xlab('Color') +
  ylab('GO Term') + 
  labs(color = 'FDR Adjusted p-Value', size = 'Gene rank')

# view plot
dotplot

# save the plot to a PDF file
ggsave('dotplotModule_mostSigGO.pdf', plot = dotplot, device = 'pdf')


