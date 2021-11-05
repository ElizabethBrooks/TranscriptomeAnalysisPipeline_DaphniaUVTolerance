# Load libraries
library(ggplot2)
library("dplyr")

#Set working directory
#workingDir = args[1];
workingDir="/Users/bamflappy/PfrenderLab/dMelUV/DEA_PA42_v4.1/glmQLFAnalysis_FDR0.10/"
setwd(workingDir); 

#Turn off scientific notation
options(scipen = 999)

#Import gene counts
inputCounts <- read.csv(file="glmQLF_2WayANOVA_interaction_topTags.csv")

#Count and visualize the distribution of FDR adjusted p-values
table(inputCounts$FDR < 0.01)

#ggplot(inputCounts, aes(x=seq_along(FDR), y=FDR)) + 
#  geom_point() +
#  geom_hline(yintercept=0.01, linetype="dashed", color = "red")

qplot(inputCounts$FDR, geom="histogram") +
  geom_vline(xintercept=0.01, linetype="dashed", color = "red")

#Transformation 1 - invert 
inputCounts$Inverted <- 1-inputCounts$FDR

#ggplot(inputCounts, aes(x=seq_along(Inverted), y=Inverted)) + 
#  geom_point() +
#  geom_hline(yintercept=0.99, linetype="dashed", color = "red")

qplot(inputCounts$Inverted, geom="histogram") +
  geom_vline(xintercept=0.99, linetype="dashed", color = "red")

#Transformation 2 - squared
inputCounts$Squared <- inputCounts$FDR^2
table(inputCounts$Squared < 0.0001)

#ggplot(inputCounts, aes(x=seq_along(Squared), y=Squared)) + 
#  geom_point() +
#  geom_hline(yintercept=0.0001, linetype="dashed", color = "red")

qplot(inputCounts$Squared, geom="histogram") +
  geom_vline(xintercept=0.0001, linetype="dashed", color = "red")

#Transformation 3 - directional nominal p-value
inputCounts$Nominal <- -log10(inputCounts$FDR) * sign(inputCounts$logFC)
table(inputCounts$Nominal > 2) + table(inputCounts$Nominal < -2)

#ggplot(inputCounts, aes(x=seq_along(Nominal), y=Nominal)) + 
#  geom_point() +
#  geom_hline(yintercept=2, linetype="dashed", color = "red") +
#  geom_hline(yintercept=-2, linetype="dashed", color = "red")

qplot(inputCounts$Nominal, geom="histogram") +
  geom_vline(xintercept=2, linetype="dashed", color = "red") +
  geom_vline(xintercept=-2, linetype="dashed", color = "red")
