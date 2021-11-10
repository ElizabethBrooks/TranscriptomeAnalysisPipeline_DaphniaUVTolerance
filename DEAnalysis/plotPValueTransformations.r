# Load libraries
library(ggplot2)
library("dplyr")

#Set working directory
#workingDir = args[1];
workingDir="/Users/bamflappy/PfrenderLab/dMelUV/DEA_PA42_v4.1/glmQLFAnalysis/"
setwd(workingDir); 

#Turn off scientific notation
options(scipen = 999)

#Import gene counts
inputCounts <- read.csv(file="glmQLF_2WayANOVA_interaction_topTags.csv")

#Count and visualize the distribution of FDR adjusted p-values
table(inputCounts$FDR < 0.05)

#Transformation 1 - reverse ordered 
inputCounts$Inverted <- 1-inputCounts$FDR
table(inputCounts$Inverted > 0.95)

#Transformation 2 - squared
inputCounts$Squared <- inputCounts$FDR^2
table(inputCounts$Squared < 0.0025)

#Transformation 3 - directional nominal p-value
inputCounts$Nominal <- -log10(inputCounts$FDR) * sign(inputCounts$logFC)
table(inputCounts$Nominal > 1.30103) + table(inputCounts$Nominal < -1.30103)

#Filter to significant genes
inputCounts.keep <- inputCounts$FDR <= 0.05
inputCounts.filt <- inputCounts[inputCounts.keep,]

#Plot p-values
jpeg("DmelUV_pValues_interactionEffect.jpeg")
ggplot(inputCounts, aes(x=seq_along(PValue), y=PValue)) + 
  geom_point() +
  geom_hline(yintercept=0.05, linetype="dashed", color = "red") +
  xlab("Rank") +
  ggtitle("P-Values by Rank") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

jpeg("DmelUV_pValues_interactionEffectSig.jpeg")
ggplot(inputCounts.filt, aes(x=seq_along(PValue), y=PValue)) + 
  geom_point() +
  xlab("Rank") +
  ggtitle("Significant P-Values by Rank") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#Plot FDR adjust p-values
jpeg("DmelUV_FDRAdjustedPValues_interactionEffect.jpeg")
ggplot(inputCounts, aes(x=seq_along(FDR), y=FDR)) + 
  geom_point() +
  geom_hline(yintercept=0.05, linetype="dashed", color = "red") +
  xlab("Rank") +
  ggtitle("FDR Adjusted P-Values by Rank") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

jpeg("DmelUV_FDRAdjustedPValues_interactionEffectSig.jpeg")
ggplot(inputCounts.filt, aes(x=seq_along(FDR), y=FDR)) + 
  geom_point() +
  xlab("Rank") +
  ggtitle("Significant FDR Adjusted P-Values by Rank") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#qplot(inputCounts$FDR, geom="histogram", binwidth = 0.001) +
#  geom_vline(xintercept=0.05, linetype="dashed", color = "red")

#Plot transformation 1
jpeg("DmelUV_FDRAdjustedPValues_reverseOrdered_interactionEffect.jpeg")
ggplot(inputCounts, aes(x=seq_along(Inverted), y=Inverted)) + 
  geom_point() +
  geom_hline(yintercept=0.95, linetype="dashed", color = "red") +
  xlab("Rank") +
  ggtitle("Reverse Ordered FDR Adjusted P-Values by Rank") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

jpeg("DmelUV_FDRAdjustedPValues_reverseOrdered_interactionEffectSig.jpeg")
ggplot(inputCounts.filt, aes(x=seq_along(Inverted), y=Inverted)) + 
  geom_point() +
  xlab("Rank") +
  ggtitle("Reverse Ordered Significant FDR Adjusted P-Values by Rank") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#qplot(inputCounts$Inverted, geom="histogram", binwidth = 0.001) +
#  geom_vline(xintercept=0.95, linetype="dashed", color = "red")

#Plot transformation 2
jpeg("DmelUV_FDRAdjustedPValues_squareTransform_interactionEffect.jpeg")
ggplot(inputCounts, aes(x=seq_along(Squared), y=Squared)) + 
  geom_point() +
  geom_hline(yintercept=0.0025, linetype="dashed", color = "red") +
  xlab("Rank") +
  ggtitle("Squared Transform of FDR Adjusted P-Values by Rank") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

jpeg("DmelUV_FDRAdjustedPValues_squaredTransform_interactionEffectSig.jpeg")
ggplot(inputCounts.filt, aes(x=seq_along(Squared), y=Squared)) + 
  geom_point() +
  xlab("Rank") +
  ggtitle("Squared Transform of Significant FDR Adjusted P-Values by Rank") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#qplot(inputCounts$Squared, geom="histogram", binwidth = 0.001) +
#  geom_vline(xintercept=0.0025, linetype="dashed", color = "red")

#Plot transformation 3
jpeg("DmelUV_FDRAdjustedPValues_directionalNominalTransform_interactionEffect.jpeg")
ggplot(inputCounts, aes(x=seq_along(Nominal), y=Nominal)) + 
  geom_point() +
  geom_hline(yintercept=1.30103, linetype="dashed", color = "red") +
  geom_hline(yintercept=-1.30103, linetype="dashed", color = "red") +
  xlab("Rank") +
  ggtitle("Directional Nominal Transform of FDR Adjusted P-Values by Rank") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

jpeg("DmelUV_FDRAdjustedPValues_directionalNominalTransform_interactionEffectSig.jpeg")
ggplot(inputCounts.filt, aes(x=seq_along(Nominal), y=Nominal)) + 
  geom_point() +
  xlab("Rank") +
  ggtitle("Directional Nominal Transform of Significant FDR Adjusted P-Values by Rank") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#qplot(inputCounts$Nominal, geom="histogram", binwidth = 0.1) +
#  geom_vline(xintercept=1.30103, linetype="dashed", color = "red") +
#  geom_vline(xintercept=-1.30103, linetype="dashed", color = "red")

