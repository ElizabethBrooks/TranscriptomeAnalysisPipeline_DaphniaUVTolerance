#!/usr/bin/env Rscript
#R script to generate scree plot for binned PCA of gene count matrices
#Library {tidyverse} used to generate scree plots
library(tidyverse)
#Retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)
#Test if there is one input argument
if (length(args)!=1) {
  stop("One file name must be supplied.n", call.=FALSE)
}
#Retrieve gene count tables with gene IDs from input file
gCount0 = read.csv(args[1], sep=",", row.names=1)
#Generate principal componants of PCA performed with prcomp
gpca <- prcomp(gCount0[ , names(gCount0) != "method"])
#Create percentage contributions of components for scree plot
data.frame(sd = gpca$sdev) %>% 
  mutate(pct = 100 * (sd/sum(sd))) %>% 
  ggplot(aes(1:4, pct)) + geom_col()