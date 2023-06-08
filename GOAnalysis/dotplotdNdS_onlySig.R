#if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
#BiocManager::install("gridExtra")

#Load libraries
library(ggplot2)
library(gridExtra)

#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)

# retrieve working directory
#workingDir <- args[1]
workingDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/SelectionTests/GOAnalysis"

# set working directory
setwd(workingDir)

# read in data on significant GO terms (BP, MF, and CC) for each selection 
positive_BP_GO_terms <- read.csv('positive_BP_sigGO_terms.csv', row.names = 1)
positive_MF_GO_terms <- read.csv('positive_MF_sigGO_terms.csv', row.names = 1)
positive_CC_GO_terms <- read.csv('positive_CC_sigGO_terms.csv', row.names = 1)

negative_BP_GO_terms <- read.csv('negative_BP_sigGO_terms.csv', row.names = 1)
negative_MF_GO_terms <- read.csv('negative_MF_sigGO_terms.csv', row.names = 1)
negative_CC_GO_terms <- read.csv('negative_CC_sigGO_terms.csv', row.names = 1)

## positive selection
# filter for top 5 significant terms
positive_BP_GO_sig <- positive_BP_GO_terms[1:5, ]
positive_MF_GO_sig <- positive_MF_GO_terms[1:5, ]
positive_CC_GO_sig <- positive_CC_GO_terms[1:5, ]

# add a column labeling the selection to each GO term set
positive_BP_plot_table <- cbind("Selection" = 'positive', positive_BP_GO_sig)
positive_MF_plot_table <- cbind("Selection" = 'positive', positive_MF_GO_sig)
positive_CC_plot_table <- cbind("Selection" = 'positive', positive_CC_GO_sig)

## negative selection
# filter for top 5 significant terms
negative_BP_GO_sig <- negative_BP_GO_terms[1:5, ]
negative_MF_GO_sig <- negative_MF_GO_terms[1:5, ]
negative_CC_GO_sig <- negative_CC_GO_terms[1:5, ]

# add a column labeling the selection to each GO term set
negative_BP_plot_table <- cbind("Selection" = 'negative', negative_BP_GO_sig)
negative_MF_plot_table <- cbind("Selection" = 'negative', negative_MF_GO_sig)
negative_CC_plot_table <- cbind("Selection" = 'negative', negative_CC_GO_sig)

# create tables of all GO terms for each level
#list_all_BP_GO_included <- unique(c(positive_BP_GO_sig$GO.ID, negative_BP_GO_sig$GO.ID))
#list_all_MF_GO_included <- unique(c(positive_MF_GO_sig$GO.ID, negative_MF_GO_sig$GO.ID))
#list_all_CC_GO_included <- unique(c(positive_CC_GO_sig$GO.ID, negative_CC_GO_sig$GO.ID))

# combine all tables
# BP
all_BP_plot_table <- rbind(positive_BP_plot_table, negative_BP_plot_table)
all_BP_plot_table <- cbind('GO_cat' = 'BP', all_BP_plot_table)
# MF
all_MF_plot_table <- rbind(positive_MF_plot_table, negative_MF_plot_table)
all_MF_plot_table <- cbind('GO_cat' = 'MF', all_MF_plot_table)
# CC
all_CC_plot_table <- rbind(positive_CC_plot_table, negative_CC_plot_table)
all_CC_plot_table <- cbind('GO_cat' = 'CC', all_CC_plot_table)
# plotting
all_plot_table <- rbind(all_BP_plot_table, all_MF_plot_table, all_CC_plot_table)

# remove NAs
plot_table <- na.omit(all_plot_table)

# create dot plot of significant GO terms
x_axis_order <- factor(plot_table$Selection, levels = c('positive', 'negative'))
facet <- factor(plot_table$GO_cat, levels = c('BP', 'MF', 'CC'))

dotplot <- ggplot(data = plot_table, aes(x = x_axis_order, y = Term, size = Significant, color = as.numeric(weightFisher))) + 
  facet_grid(rows = facet, space = 'free_y', scales = 'free') +
  geom_point() +
  scale_color_gradientn(colors = heat.colors(10), limits=c(0, 0.05)) + 
  #scale_color_gradientn(colors = wes_palette("Zissou1", type = "continuous")) +
  theme_bw() +
  xlab('Selection') +
  ylab('GO Term') + 
  labs(color = 'FDR Adjusted p-Value', size = 'Gene rank')

# view plot
dotplot

# save the plot to a PDF file
ggsave('dotplot_sigGO.pdf', plot = dotplot, device = 'pdf')

# subset positively selected genes
plot_tableSubset <- plot_table[plot_table$Selection == "positive",]

# create dot plot of significant GO terms
x_axis_order <- factor(plot_tableSubset$Selection, levels = c('positive'))
facet <- factor(plot_tableSubset$GO_cat, levels = c('BP', 'MF', 'CC'))

dotplotSubset <- ggplot(data = plot_tableSubset, aes(x = x_axis_order, y = Term, size = Significant, color = as.numeric(weightFisher))) + 
  facet_grid(rows = facet, space = 'free_y', scales = 'free') +
  geom_point() +
  scale_color_gradientn(colors = heat.colors(10), limits=c(0, 0.05)) + 
  #scale_color_gradientn(colors = wes_palette("Zissou1", type = "continuous")) +
  theme_bw() +
  xlab('Selection') +
  ylab('GO Term') + 
  labs(color = 'FDR Adjusted p-Value', size = 'Gene rank')

# view plot
dotplotSubset

# save the plot to a PDF file
ggsave('dotplotPositive_sigGO.pdf', plot = dotplotSubset, device = 'pdf')

