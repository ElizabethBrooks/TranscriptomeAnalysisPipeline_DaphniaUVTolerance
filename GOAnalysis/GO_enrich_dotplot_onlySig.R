#if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
#BiocManager::install('______')

#Load libraries
library(ggplot2)
library(gridExtra)

#Read in data on GO terms (BP, MF, and CC) for each genotype 
Y05_BP_GO_terms <- read.csv('Y05_all_genes/Y05_all_genes_BP_GO_terms.csv', row.names = 1)
Y05_MF_GO_terms <- read.csv('Y05_all_genes/Y05_all_genes_MF_GO_terms.csv', row.names = 1)
Y05_CC_GO_terms <- read.csv('Y05_all_genes/Y05_all_genes_CC_GO_terms.csv', row.names = 1)

Y023_BP_GO_terms <- read.csv('Y023_all_genes/Y023_all_genes_BP_GO_terms.csv', row.names = 1)
Y023_MF_GO_terms <- read.csv('Y023_all_genes/Y023_all_genes_MF_GO_terms.csv', row.names = 1)
Y023_CC_GO_terms <- read.csv('Y023_all_genes/Y023_all_genes_CC_GO_terms.csv', row.names = 1)

E05_BP_GO_terms <- read.csv('E05_all_genes/E05_all_genes_BP_GO_terms.csv', row.names = 1)
E05_MF_GO_terms <- read.csv('E05_all_genes/E05_all_genes_MF_GO_terms.csv', row.names = 1)
E05_CC_GO_terms <- read.csv('E05_all_genes/E05_all_genes_CC_GO_terms.csv', row.names = 1)

R2_BP_GO_terms <- read.csv('R2_all_genes/R2_all_genes_BP_GO_terms.csv', row.names = 1)
R2_MF_GO_terms <- read.csv('R2_all_genes/R2_all_genes_MF_GO_terms.csv', row.names = 1)
R2_CC_GO_terms <- read.csv('R2_all_genes/R2_all_genes_CC_GO_terms.csv', row.names = 1)

PA_BP_GO_terms <- read.csv('PA_all_genes/PA_all_genes_BP_GO_terms.csv', row.names = 1)
PA_MF_GO_terms <- read.csv('PA_all_genes/PA_all_genes_MF_GO_terms.csv', row.names = 1)
PA_CC_GO_terms <- read.csv('PA_all_genes/PA_all_genes_CC_GO_terms.csv', row.names = 1)

Sierra_BP_GO_terms <- read.csv('Sierra_all_genes/Sierra_all_genes_BP_GO_terms.csv', row.names = 1)
Sierra_MF_GO_terms <- read.csv('Sierra_all_genes/Sierra_all_genes_MF_GO_terms.csv', row.names = 1)
Sierra_CC_GO_terms <- read.csv('Sierra_all_genes/Sierra_all_genes_CC_GO_terms.csv', row.names = 1)

#--------------------------------------------------------------------------------------------------

#ALL BP WORK
    #Filter for top 5 significant in each (get the 1st 5 rows...then take only the significant ones).
    #They should already be in decreasing significance order.
Y05_BP_GO_top5 <- Y05_BP_GO_terms[1:5, ]
Y023_BP_GO_top5 <- Y023_BP_GO_terms[1:5, ]
E05_BP_GO_top5 <- E05_BP_GO_terms[1:5, ]
R2_BP_GO_top5 <- R2_BP_GO_terms[1:5, ]
PA_BP_GO_top5 <- PA_BP_GO_terms[1:5, ]
Sierra_BP_GO_top5 <- Sierra_BP_GO_terms[1:5, ]

Y05_BP_GO_sig <- Y05_BP_GO_top5[which(Y05_BP_GO_top5$p_adjusted <= 0.05), ]
Y023_BP_GO_sig <- Y023_BP_GO_top5[which(Y023_BP_GO_top5$p_adjusted <= 0.05), ]
E05_BP_GO_sig <- E05_BP_GO_top5[which(E05_BP_GO_top5$p_adjusted <= 0.05), ]
R2_BP_GO_sig <- R2_BP_GO_top5[which(R2_BP_GO_top5$p_adjusted <= 0.05), ]
PA_BP_GO_sig <- PA_BP_GO_top5[which(PA_BP_GO_top5$p_adjusted <= 0.05), ]
Sierra_BP_GO_sig <- Sierra_BP_GO_top5[which(Sierra_BP_GO_top5$p_adjusted <= 0.05), ]

    #List of all 30 BP for each genotype (duplicates removed)
list_all_BP_GO_included <- unique(c(Y05_BP_GO_sig$GO.ID, Y023_BP_GO_sig$GO.ID, E05_BP_GO_sig$GO.ID, 
                                    R2_BP_GO_sig$GO.ID, PA_BP_GO_sig$GO.ID, Sierra_BP_GO_sig$GO.ID))


    #Add a column labeling the genotype to each 
Y05_BP_plot_table <- cbind("Genotype" = 'Y05', Y05_BP_GO_sig)
Y023_BP_plot_table <- cbind("Genotype" = 'Y023', Y023_BP_GO_sig)
E05_BP_plot_table <- cbind("Genotype" = 'E05', E05_BP_GO_sig)
#R2_BP_plot_table <- cbind("Genotype" = 'R2', R2_BP_GO_sig)
PA_BP_plot_table <- cbind("Genotype" = 'PA', PA_BP_GO_sig)
Sierra_BP_plot_table <- cbind("Genotype" = 'Sierra', Sierra_BP_GO_sig)

#--------------------------------------------------------------------------------------------------

#ALL MF WORK
    #Filter for top 5 significant in each (get the 1st 5 rows...then take only the significant ones).
    #They should already be in decreasing significance order.
Y05_MF_GO_top5 <- Y05_MF_GO_terms[1:5, ]
Y023_MF_GO_top5 <- Y023_MF_GO_terms[1:5, ]
E05_MF_GO_top5 <- E05_MF_GO_terms[1:5, ]
R2_MF_GO_top5 <- R2_MF_GO_terms[1:5, ]
PA_MF_GO_top5 <- PA_MF_GO_terms[1:5, ]
Sierra_MF_GO_top5 <- Sierra_MF_GO_terms[1:5, ]


Y05_MF_GO_sig <- Y05_MF_GO_top5[which(Y05_MF_GO_top5$p_adjusted <= 0.05), ]
Y023_MF_GO_sig <- Y023_MF_GO_top5[which(Y023_MF_GO_top5$p_adjusted <= 0.05), ]
E05_MF_GO_sig <- E05_MF_GO_top5[which(E05_MF_GO_top5$p_adjusted <= 0.05), ]
R2_MF_GO_sig <- R2_MF_GO_top5[which(R2_MF_GO_top5$p_adjusted <= 0.05), ]
PA_MF_GO_sig <- PA_MF_GO_top5[which(PA_MF_GO_top5$p_adjusted <= 0.05), ]
Sierra_MF_GO_sig <- Sierra_MF_GO_top5[which(Sierra_MF_GO_top5$p_adjusted <= 0.05), ]

    #List of all 30 MF for each genotype (duplicates removed)
list_all_MF_GO_included <- unique(c(Y05_MF_GO_sig$GO.ID, Y023_MF_GO_sig$GO.ID, E05_MF_GO_sig$GO.ID, 
                                    R2_MF_GO_sig$GO.ID, PA_MF_GO_sig$GO.ID, Sierra_MF_GO_sig$GO.ID))


    #Add a column labeling the genotype to each 
Y05_MF_plot_table <- cbind("Genotype" = 'Y05', Y05_MF_GO_sig)
Y023_MF_plot_table <- cbind("Genotype" = 'Y023', Y023_MF_GO_sig)
E05_MF_plot_table <- cbind("Genotype" = 'E05', E05_MF_GO_sig)
R2_MF_plot_table <- cbind("Genotype" = 'R2', R2_MF_GO_sig)
PA_MF_plot_table <- cbind("Genotype" = 'PA', PA_MF_GO_sig)
Sierra_MF_plot_table <- cbind("Genotype" = 'Sierra', Sierra_MF_GO_sig)

#--------------------------------------------------------------------------------------------------

#ALL CC WORK
    #Filter for top 5 significant in each (get the 1st 5 rows...then take only the significant ones).
    #They should already be in decreasing significance order.
Y05_CC_GO_top5 <- Y05_CC_GO_terms[1:5, ]
Y023_CC_GO_top5 <- Y023_CC_GO_terms[1:5, ]
E05_CC_GO_top5 <- E05_CC_GO_terms[1:5, ]
R2_CC_GO_top5 <- R2_CC_GO_terms[1:5, ]
PA_CC_GO_top5 <- PA_CC_GO_terms[1:5, ]
Sierra_CC_GO_top5 <- Sierra_CC_GO_terms[1:5, ]


Y05_CC_GO_sig <- Y05_CC_GO_top5[which(Y05_CC_GO_top5$p_adjusted <= 0.05), ]
Y023_CC_GO_sig <- Y023_CC_GO_top5[which(Y023_CC_GO_top5$p_adjusted <= 0.05), ]
E05_CC_GO_sig <- E05_CC_GO_top5[which(E05_CC_GO_top5$p_adjusted <= 0.05), ]
R2_CC_GO_sig <- R2_CC_GO_top5[which(R2_CC_GO_top5$p_adjusted <= 0.05), ]
PA_CC_GO_sig <- PA_CC_GO_top5[which(PA_CC_GO_top5$p_adjusted <= 0.05), ]
Sierra_CC_GO_sig <- Sierra_CC_GO_top5[which(Sierra_CC_GO_top5$p_adjusted <= 0.05), ]

    #List of all 30 CC for each genotype (duplicates removed)
list_all_CC_GO_included <- unique(c(Y05_CC_GO_sig$GO.ID, Y023_CC_GO_sig$GO.ID, E05_CC_GO_sig$GO.ID, 
                                    R2_CC_GO_sig$GO.ID, PA_CC_GO_sig$GO.ID, Sierra_CC_GO_sig$GO.ID))


    #Add a column labeling the genotype to each 
Y05_CC_plot_table <- cbind("Genotype" = 'Y05', Y05_CC_GO_sig)
Y023_CC_plot_table <- cbind("Genotype" = 'Y023', Y023_CC_GO_sig)
E05_CC_plot_table <- cbind("Genotype" = 'E05', E05_CC_GO_sig)
#R2_CC_plot_table <- cbind("Genotype" = 'R2', R2_CC_GO_sig)
PA_CC_plot_table <- cbind("Genotype" = 'PA', PA_CC_GO_sig)
Sierra_CC_plot_table <- cbind("Genotype" = 'Sierra', Sierra_CC_GO_sig)

#--------------------------------------------------------------------------------------------------

#Combine all tables into 1
all_genotype_BP_plot_table <- rbind(Y05_BP_plot_table, Y023_BP_plot_table, E05_BP_plot_table, PA_BP_plot_table, 
                                    Sierra_BP_plot_table)
all_genotype_BP_plot_table <- cbind('GO_cat' = 'BP', all_genotype_BP_plot_table)

all_genotype_MF_plot_table <- rbind(Y05_MF_plot_table, Y023_MF_plot_table, E05_MF_plot_table, R2_MF_plot_table, PA_MF_plot_table, 
                                    Sierra_MF_plot_table)
all_genotype_MF_plot_table <- cbind('GO_cat' = 'MF', all_genotype_MF_plot_table)

all_genotype_CC_plot_table <- rbind(Y05_CC_plot_table, Y023_CC_plot_table, E05_CC_plot_table, PA_CC_plot_table, 
                                    Sierra_CC_plot_table)
all_genotype_CC_plot_table <- cbind('GO_cat' = 'CC', all_genotype_CC_plot_table)

all_genotype_plot_table <- rbind(all_genotype_BP_plot_table, all_genotype_MF_plot_table, all_genotype_CC_plot_table)

#--------------------------------------------------------------------------------------------------

#Make plot
x_axis_order <- factor(all_genotype_plot_table$Genotype, levels = c('Y05', 'Y023', 'E05', 'R2', 'PA', 'Sierra'))
facet <- factor(all_genotype_plot_table$GO_cat, levels = c('BP', 'MF', 'CC'))

p <- ggplot(data = all_genotype_plot_table, aes(x = x_axis_order, y = Term, size = Significant, color = p_adjusted)) 
p <- p + facet_grid(rows = facet, space = 'free_y', scales = 'free')
p <- p + geom_point() + scale_color_gradientn(colors = heat.colors(10), limits=c(0, 0.05)) + theme_bw()
p <- p + xlab('Genotype') + ylab('GO Term') 
p <- p + labs(color = 'Adjusted p-value', size = 'Gene rank')
#--------------------------------------------------------------------------------------------------

#Save to file
final_plot <- p
ggsave('GOenrich_Dotplot_CustomAnnotation_AllGenes.pdf', plot = final_plot, device = 'pdf')


