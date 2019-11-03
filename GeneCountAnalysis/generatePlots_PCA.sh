#!/bin/bash
#Usage: bash generatePlots_PCA.sh mergedCounts_file_transposed.csv
#Usage Ex: bash generatePlots_PCA.sh mergedCounts_legacy_transposed.csv
#Script to run Rscripts that generate PCA plots

#Plot merged data PCA
Rscript geneCounts_PCA.r $1
#Rename produced plot
outFile=$(echo $1 | sed 's/\.csv//')
mv Rplots.pdf $outFile".pdf"