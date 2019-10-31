#!/bin/bash
#Usage: bash generatePlots_PCA.sh mergedCounts_file_transposed.csv SCALE
#Usage Ex: bash generatePlots_PCA.sh mergedCounts_legacy_transposed.csv FALSE
#Script to run Rscripts that generate PCA plots

#Plot merged data PCA
Rscript geneCounts_PCA.r $1 $2
#Rename produced plot
outFile=$(echo $1 | sed 's/\.csv//')
mv Rplots.pdf $outFile"_scale"$2".pdf"