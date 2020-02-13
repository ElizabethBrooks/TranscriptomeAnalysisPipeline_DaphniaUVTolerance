#!/usr/bin/env Rscript
#Usage: Rscript 
#Usage Ex: Rscript 
#R script to generate grouped and colored bar plots

library("png") # for reading in PNGs

# example image
img <- readPNG(system.file("img", "Rlogo.png", package="png"))

# setup plot
par(mar=rep(0,4)) # no margins

# layout the plots into a matrix w/ 12 columns, by row
layout(matrix(1:120, ncol=12, byrow=TRUE))

# do the plotting
for(i in 1:120) {
    plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
    rasterImage(img,0,0,1,1)
}

# write to PDF
dev.print(pdf, "output.pdf")