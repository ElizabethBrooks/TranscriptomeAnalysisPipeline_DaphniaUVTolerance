#if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
#BiocManager::install("gridExtra")

#Load libraries
library(ggplot2)
library(gridExtra)

#Retrieve input file name of gene counts
#args = commandArgs(trailingOnly=TRUE)

# retrieve working directory
#workingDir <- args[1]
workingDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Tolerance/GOAnalysis"

# set working directory
setwd(workingDir)

# read in data on significant GO terms (BP, MF, and CC) for each effect 
interaction_BP_GO_terms <- read.csv('interaction_BP_top5GO_terms.csv', row.names = 1)
interaction_MF_GO_terms <- read.csv('interaction_MF_top5GO_terms.csv', row.names = 1)
interaction_CC_GO_terms <- read.csv('interaction_CC_top5GO_terms.csv', row.names = 1)

treatment_BP_GO_terms <- read.csv('treatment_BP_top5GO_terms.csv', row.names = 1)
treatment_MF_GO_terms <- read.csv('treatment_MF_top5GO_terms.csv', row.names = 1)
treatment_CC_GO_terms <- read.csv('treatment_CC_top5GO_terms.csv', row.names = 1)

tolerance_BP_GO_terms <- read.csv('tolerance_BP_top5GO_terms.csv', row.names = 1)
tolerance_MF_GO_terms <- read.csv('tolerance_MF_top5GO_terms.csv', row.names = 1)
tolerance_CC_GO_terms <- read.csv('tolerance_CC_top5GO_terms.csv', row.names = 1)

## interaction effect
# filter for top 5 significant terms
interaction_BP_GO_top5 <- interaction_BP_GO_terms[1:5, ]
interaction_MF_GO_top5 <- interaction_MF_GO_terms[1:5, ]
interaction_CC_GO_top5 <- interaction_CC_GO_terms[1:5, ]

# add a column labeling the effect to each GO term set
interaction_BP_plot_table <- cbind("Effect" = 'interaction', interaction_BP_GO_top5)
interaction_MF_plot_table <- cbind("Effect" = 'interaction', interaction_MF_GO_top5)
interaction_CC_plot_table <- cbind("Effect" = 'interaction', interaction_CC_GO_top5)

## treatment effect
# filter for top 5 significant terms
treatment_BP_GO_top5 <- treatment_BP_GO_terms[1:5, ]
treatment_MF_GO_top5 <- treatment_MF_GO_terms[1:5, ]
treatment_CC_GO_top5 <- treatment_CC_GO_terms[1:5, ]

# add a column labeling the effect to each GO term set
treatment_BP_plot_table <- cbind("Effect" = 'treatment', treatment_BP_GO_top5)
treatment_MF_plot_table <- cbind("Effect" = 'treatment', treatment_MF_GO_top5)
treatment_CC_plot_table <- cbind("Effect" = 'treatment', treatment_CC_GO_top5)

## tolerance effect
# filter for top 5 significant terms
tolerance_BP_GO_top5 <- tolerance_BP_GO_terms[1:5, ]
tolerance_MF_GO_top5 <- tolerance_MF_GO_terms[1:5, ]
tolerance_CC_GO_top5 <- tolerance_CC_GO_terms[1:5, ]

# add a column labeling the effect to each GO term set
tolerance_BP_plot_table <- cbind("Effect" = 'tolerance', tolerance_BP_GO_top5)
tolerance_MF_plot_table <- cbind("Effect" = 'tolerance', tolerance_MF_GO_top5)
tolerance_CC_plot_table <- cbind("Effect" = 'tolerance', tolerance_CC_GO_top5)

# create tables of all GO terms for each level
list_all_BP_GO_included <- unique(c(interaction_BP_GO_top5$GO.ID, treatment_BP_GO_top5$GO.ID, tolerance_BP_GO_top5$GO.ID))
list_all_MF_GO_included <- unique(c(interaction_MF_GO_top5$GO.ID, treatment_MF_GO_top5$GO.ID, tolerance_MF_GO_top5$GO.ID))
list_all_CC_GO_included <- unique(c(interaction_CC_GO_top5$GO.ID, treatment_CC_GO_top5$GO.ID, tolerance_CC_GO_top5$GO.ID))

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
x_axis_order <- factor(plot_table$Effect, levels = c('interaction', 'treatment', 'tolerance'))
facet <- factor(plot_table$GO_cat, levels = c('BP', 'MF', 'CC'))

dotplot <- ggplot(data = plot_table, aes(x = x_axis_order, y = Term, size = Significant, color = weightFisher)) + 
  facet_grid(rows = facet, space = 'free_y', scales = 'free') +
  geom_point() +
  scale_color_gradientn(colors = heat.colors(10), limits=c(0, 0.05)) + 
  theme_bw() +
  xlab('Effect') +
  ylab('GO Term') + 
  labs(color = 'FDR Adjusted p-Value', size = 'Gene rank')

# view plot
dotplot

# save the plot to a PDF file
ggsave('GOenrich_Dotplot_CustomAnnotation.pdf', plot = dotplot, device = 'pdf')


