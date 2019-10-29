#!/bin/bash
#Script to run Rscripts that generate kmeans PCA plots

#Plot merged data kmeans PCA
Rscript geneCounts_kMeans.r $1
#Rename produced plot
outFile=$(echo $1 | sed 's/\.csv/\.pdf/')
mv Rplots.pdf $outFile