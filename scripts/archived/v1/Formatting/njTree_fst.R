#!/usr/bin/env Rscript

# script to create a NJ tree from Fst values

# set the working directory
setwd("/Users/bamflappy/PfrenderLab/OLYM_dMelUV/OLYM_past")

# install packages
install.packages("ape")

# load libraries
library(ape)

# create vector of Fst values from the lower triangle of the data mtrix
x <- c(0.12, 0.22, 0.02, 0.13, 0.30, 0.19, 
       0.42, 0.18, 0.26, 0.51, 0.36,
       0.28, 0.14, 0.14, 0.14,
       0.19, 0.30, 0.16, 
       0.27, 0.16, 
       0.07)

# create a matrix of zeros
M <- matrix(0, 7, 7)

# add the Fst values to the lower triangle
M[lower.tri(M)] <- x

# transpose
M <- t(M)

# add the Fst values to the lower triangle
M[lower.tri(M)] <- x

# add the column and row names
dimnames(M) <- list(c("H","C (E05)","G","A (R2)","B","K (Y05)","I (Y023)"), c("H","C (E05)","G","A (R2)","B","K (Y05)","I (Y023)"))

# perform neighbor joining
tr <- nj(M)

# create and save the NJ tree plot
jpeg("olympicsFst_njTree.jpg", width = 700, height = 350)
plot(tr, "u")
dev.off()
