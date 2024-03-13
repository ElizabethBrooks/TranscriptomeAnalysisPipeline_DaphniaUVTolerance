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

#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)

# retrieve working directory
workingDir <- args[1]
#workingDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Genotypes/GOAnalysis_ks"

# set working directory
setwd(workingDir)

# retrieve positive selection enrichment tests
positiveTable <- args[2]
#positiveTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/selectionTests/Genotypes/fisherTest_positiveSelection_DEGs.csv", row.names="effect")

# read in data on significant GO terms (BP, MF, and CC) for each effect 
interaction_BP_GO_terms <- read.csv('interaction_BP_sigGO_terms.csv', row.names = 1)
interaction_MF_GO_terms <- read.csv('interaction_MF_sigGO_terms.csv', row.names = 1)
interaction_CC_GO_terms <- read.csv('interaction_CC_sigGO_terms.csv', row.names = 1)

treatment_BP_GO_terms <- read.csv('treatment_BP_sigGO_terms.csv', row.names = 1)
treatment_MF_GO_terms <- read.csv('treatment_MF_sigGO_terms.csv', row.names = 1)
treatment_CC_GO_terms <- read.csv('treatment_CC_sigGO_terms.csv', row.names = 1)

tolerance_BP_GO_terms <- read.csv('tolerance_BP_sigGO_terms.csv', row.names = 1)
tolerance_MF_GO_terms <- read.csv('tolerance_MF_sigGO_terms.csv', row.names = 1)
tolerance_CC_GO_terms <- read.csv('tolerance_CC_sigGO_terms.csv', row.names = 1)

## interaction effect
# filter for top 5 significant terms
interaction_BP_GO_sig <- interaction_BP_GO_terms[1:5, ]
interaction_MF_GO_sig <- interaction_MF_GO_terms[1:5, ]
interaction_CC_GO_sig <- interaction_CC_GO_terms[1:5, ]

# add a column labeling the effect to each GO term set
interaction_BP_plot_table <- cbind("Effect" = 'Interaction', interaction_BP_GO_sig)
interaction_MF_plot_table <- cbind("Effect" = 'Interaction', interaction_MF_GO_sig)
interaction_CC_plot_table <- cbind("Effect" = 'Interaction', interaction_CC_GO_sig)

## treatment effect
# filter for top 5 significant terms
treatment_BP_GO_sig <- treatment_BP_GO_terms[1:5, ]
treatment_MF_GO_sig <- treatment_MF_GO_terms[1:5, ]
treatment_CC_GO_sig <- treatment_CC_GO_terms[1:5, ]

# add a column labeling the effect to each GO term set
treatment_BP_plot_table <- cbind("Effect" = 'Treatment', treatment_BP_GO_sig)
treatment_MF_plot_table <- cbind("Effect" = 'Treatment', treatment_MF_GO_sig)
treatment_CC_plot_table <- cbind("Effect" = 'Treatment', treatment_CC_GO_sig)

## tolerance effect
# filter for top 5 significant terms
tolerance_BP_GO_sig <- tolerance_BP_GO_terms[1:5, ]
tolerance_MF_GO_sig <- tolerance_MF_GO_terms[1:5, ]
tolerance_CC_GO_sig <- tolerance_CC_GO_terms[1:5, ]

# add a column labeling the effect to each GO term set
tolerance_BP_plot_table <- cbind("Effect" = 'Tolerance', tolerance_BP_GO_sig)
tolerance_MF_plot_table <- cbind("Effect" = 'Tolerance', tolerance_MF_GO_sig)
tolerance_CC_plot_table <- cbind("Effect" = 'Tolerance', tolerance_CC_GO_sig)

# create tables of all GO terms for each level
#list_all_BP_GO_included <- unique(c(interaction_BP_GO_sig$GO.ID, treatment_BP_GO_sig$GO.ID, tolerance_BP_GO_sig$GO.ID))
#list_all_MF_GO_included <- unique(c(interaction_MF_GO_sig$GO.ID, treatment_MF_GO_sig$GO.ID, tolerance_MF_GO_sig$GO.ID))
#list_all_CC_GO_included <- unique(c(interaction_CC_GO_sig$GO.ID, treatment_CC_GO_sig$GO.ID, tolerance_CC_GO_sig$GO.ID))

# combine all tables
# BP
all_BP_plot_table <- rbind(interaction_BP_plot_table, treatment_BP_plot_table, tolerance_BP_plot_table)
all_BP_plot_table <- cbind('GO_cat' = 'BP', all_BP_plot_table)
# MF
all_MF_plot_table <- rbind(interaction_MF_plot_table, treatment_MF_plot_table, tolerance_MF_plot_table)
all_MF_plot_table <- cbind('GO_cat' = 'MF', all_MF_plot_table)
# CC
all_CC_plot_table <- rbind(interaction_CC_plot_table, treatment_CC_plot_table, tolerance_CC_plot_table)
all_CC_plot_table <- cbind('GO_cat' = 'CC', all_CC_plot_table)
# plotting
all_plot_table <- rbind(all_BP_plot_table, all_MF_plot_table, all_CC_plot_table)

# remove NAs
plot_table <- na.omit(all_plot_table)

# create dot plot of significant GO terms
x_axis_order <- factor(plot_table$Effect, levels = c('Interaction', 'Treatment', 'Tolerance'))
facet <- factor(plot_table$GO_cat, levels = c('BP', 'CC', 'MF'))

dotplot <- ggplot(data = plot_table, aes(x = x_axis_order, y = Term, size = Significant, color = as.numeric(weightFisher))) + 
  facet_grid(rows = facet, space = 'free_y', scales = 'free') +
  geom_point() +
  #scale_color_gradientn(colors = heat.colors(10), limits=c(0, 0.05)) + 
  scale_color_gradientn(colors = plotColorSubset) +
  #scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_bw() +
  xlab('Effect') +
  ylab('GO Term') + 
  #scale_x_discrete(labels=c("Interaction"=expression(italic("Interaction")), parse=TRUE)) +
  labs(color = 'P-Value', size = 'Gene Rank')

# view plot
dotplot

# save the plot to a PDF file
ggsave('dotplot_sigGO.pdf', plot = dotplot, device = 'pdf')


