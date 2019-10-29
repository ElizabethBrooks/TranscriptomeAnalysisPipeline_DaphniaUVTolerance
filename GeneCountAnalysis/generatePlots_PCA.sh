#!/bin/bash
#Script to run Rscripts that generate PCA plots

#Plot merged data PCA
Rscript geneCounts_PCA.r $1
#Rename produced plot
outFile=$(echo $1 | sed 's/\.csv/\.pdf/')
mv Rplots.pdf $outFile