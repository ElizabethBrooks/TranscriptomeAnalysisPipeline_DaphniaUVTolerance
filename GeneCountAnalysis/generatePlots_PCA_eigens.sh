#!/bin/bash
#Script to run Rscripts that generate binned PCA plots

#Plot merged data binned PCA
Rscript geneCounts_PCA_eigens.r $1
#Rename produced plot
outFile=$(echo $1 | sed 's/\.csv/\.pdf/')
mv Rplots.pdf $outFile