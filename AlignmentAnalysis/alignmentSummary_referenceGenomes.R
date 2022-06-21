#Import libraries
library(ggplot2)
library(stringr)

#Read in alignment summaries
align_KAP4 <- read.csv("/Users/bamflappy/PfrenderLab/dMelUV/AlignmentAnalysis/alignmentSummarized_run1_hisat2_run1_KAP4.csv")
align_PA42_v3.0 <- read.csv("/Users/bamflappy/PfrenderLab/dMelUV/AlignmentAnalysis/alignmentSummarized_run2_hisat2_run2_PA42_v3.0.csv")
align_PA42_v4.1 <- read.csv("/Users/bamflappy/PfrenderLab/dMelUV/AlignmentAnalysis/alignmentSummarized_run4_hisat2_run4_PA42_v4.1.csv")

#Collect alignment stats
(align_stats <- data.frame(
  sample = c(str_remove(align_KAP4$sample, "140327_I481_FCC3P1PACXX_L[234]_"), str_remove(align_PA42_v4.1$sample, "140327_I481_FCC3P1PACXX_L[234]_"), str_remove(align_PA42_v3.0$sample, "140327_I481_FCC3P1PACXX_L[234]_")),
  ref = c(rep("KAP4", each=nrow(align_KAP4)), rep("PA42_v4.1", each=nrow(align_PA42_v4.1)), rep("PA42_v3.0", each=nrow(align_PA42_v3.0))),
  overall = sapply(c(str_remove(align_KAP4$overall, "%"), str_remove(align_PA42_v4.1$overall, "%"), str_remove(align_PA42_v3.0$overall, "%")), as.numeric),
  concordant = sapply(c(str_remove(align_KAP4$concordant, "%"), str_remove(align_PA42_v4.1$concordant, "%"), str_remove(align_PA42_v3.0$concordant, "%")), as.numeric)
))

#Check data types in DF
print(summary(align_stats))

#Barplot of overall alignment stats
align_plot_overall <- ggplot(align_stats, aes(x=sample, y=overall, fill=ref)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5))+
  labs(title="Alignment Percentages to D. Pulex Reference Genomes", x="Daphnia Sample", y="Overall Alignment")
ggsave("/Users/bamflappy/PfrenderLab/dMelUV/AlignmentAnalysis/align_stats_ref_overall.png", align_plot_overall)

#Barplot of concordant alignment stats
align_plot_conc <- ggplot(align_stats, aes(x=sample, y=concordant, fill=ref)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5))+
  labs(title="Alignment Percentages to D. Pulex Reference Genomes", x="Daphnia Sample", y="Uniquely Mapped")
ggsave("/Users/bamflappy/PfrenderLab/dMelUV/AlignmentAnalysis/align_stats_ref_conc.png", align_plot_conc)




####BAD EXAMPLE####
#Collect alignment stats
(align_stats_bad <- data.frame(
  sample = c(str_remove(align_KAP4$sample, "140327_I481_FCC3P1PACXX_L[234]_"), str_remove(align_PA42_v4.1$sample, "140327_I481_FCC3P1PACXX_L[234]_"), str_remove(align_PA42_v3.0$sample, "140327_I481_FCC3P1PACXX_L[234]_")),
  ref = c(rep("KAP4", each=nrow(align_KAP4)), rep("PA42_v4.1", each=nrow(align_PA42_v4.1)), rep("PA42_v3.0", each=nrow(align_PA42_v3.0))),
  overall = c(str_remove(align_KAP4$overall, "%"), str_remove(align_PA42_v4.1$overall, "%"), str_remove(align_PA42_v3.0$overall, "%")),
  concordant = c(str_remove(align_KAP4$concordant, "%"), str_remove(align_PA42_v4.1$concordant, "%"), str_remove(align_PA42_v3.0$concordant, "%"))
))

#Check data types in DF
print(summary(align_stats_bad))

#Barplot of alignment stats
align_plot_bad <- ggplot(align_stats, aes(x=sample, y=overall, fill=ref)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title="Alignment Percentages to D. Pulex References", x="D. melanica Sample", y="Overall Alignment")

