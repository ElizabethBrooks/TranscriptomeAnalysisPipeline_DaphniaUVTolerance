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
inputCounts <- read.csv(file="glmQLF_2WayANOVA_interaction_topTags_filtered.csv")

#Count and visualize the distribution of FDR adjusted p-values
table(inputCounts$FDR < 0.05)

#ggplot(inputCounts, aes(x=seq_along(FDR), y=FDR)) + 
#  geom_point() +
#  geom_hline(yintercept=0.05, linetype="dashed", color = "red")

#ggplot(inputCounts, aes(FDR)) +
#  geom_bar() +
#  scale_x_binned() +
#  geom_vline(xintercept=0.05, linetype="dashed", color = "red")

qplot(inputCounts$FDR, geom="histogram", binwidth = 0.001) +
  geom_vline(xintercept=0.05, linetype="dashed", color = "red")

#Transformation 1 - invert 
inputCounts$Inverted <- 1-inputCounts$FDR
table(inputCounts$Inverted > 0.95)

#ggplot(inputCounts, aes(x=seq_along(Inverted), y=Inverted)) + 
#  geom_point() +
#  geom_hline(yintercept=0.95, linetype="dashed", color = "red")

#ggplot(inputCounts, aes(Inverted)) +
#  geom_bar() +
#  scale_x_binned() +
#  geom_vline(xintercept=0.95, linetype="dashed", color = "red")

qplot(inputCounts$Inverted, geom="histogram", binwidth = 0.001) +
  geom_vline(xintercept=0.95, linetype="dashed", color = "red")

#Transformation 2 - squared
inputCounts$Squared <- inputCounts$FDR^2
table(inputCounts$Squared < 0.0025)

#ggplot(inputCounts, aes(x=seq_along(Squared), y=Squared)) + 
#  geom_point() +
#  geom_hline(yintercept=0.0025, linetype="dashed", color = "red")

#ggplot(inputCounts, aes(Squared)) +
#  geom_bar() +
#  scale_x_binned() +
#  geom_vline(xintercept=0.0025, linetype="dashed", color = "red")

qplot(inputCounts$Squared, geom="histogram", binwidth = 0.0001) +
  geom_vline(xintercept=0.0025, linetype="dashed", color = "red")

#Transformation 3 - directional nominal p-value
inputCounts$Nominal <- -log10(inputCounts$FDR) * sign(inputCounts$logFC)
table(inputCounts$Nominal > 1.30103) + table(inputCounts$Nominal < -1.30103)

#ggplot(inputCounts, aes(x=seq_along(Nominal), y=Nominal)) + 
#  geom_point() +
#  geom_hline(yintercept=1.30103, linetype="dashed", color = "red") +
#  geom_hline(yintercept=-1.30103, linetype="dashed", color = "red")

#ggplot(inputCounts, aes(Nominal)) +
#  geom_bar(binwidth = 0.5) +
#  scale_x_binned() +
#  geom_vline(xintercept=1.30103, linetype="dashed", color = "red") +
#  geom_vline(xintercept=-1.30103, linetype="dashed", color = "red")

qplot(inputCounts$Nominal, geom="histogram", binwidth = 0.1) +
  geom_vline(xintercept=1.30103, linetype="dashed", color = "red") +
  geom_vline(xintercept=-1.30103, linetype="dashed", color = "red")

