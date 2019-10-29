#!/bin/bash
#Script to run Rscripts that generate binned kmeans PCA plots

#Plot merged data binned kmeans PCA
Rscript geneCounts_kMeans_binned.r $1
#Rename produced plot
outFile=$(echo $1 | sed 's/\.csv/\.pdf/')
mv Rplots.pdf $outFile