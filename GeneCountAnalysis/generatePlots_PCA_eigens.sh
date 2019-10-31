#!/bin/bash
#Usage: bash generatePlots_PCA_eigens.sh mergedCounts_file_transposed.csv SCALE
#Usage Ex: bash generatePlots_PCA_eigens.sh mergedCounts_legacy_transposed.csv FALSE
#Script to run Rscripts that generate binned PCA plots

#Plot merged data binned PCA
Rscript geneCounts_PCA_eigens.r $1 $2
#Rename produced plot
outFile=$(echo $1 | sed 's/\.csv//')
mv Rplots.pdf $outFile"_scale"$2".pdf"