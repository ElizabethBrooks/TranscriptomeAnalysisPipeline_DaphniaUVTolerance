#!/bin/bash
#Usage: bash generatePlots_PCA_eigens_binned.sh mergedCounts_file_transposed.csv bin SCALE
#Usage Ex: bash generatePlots_PCA_eigens_binned.sh mergedCounts_legacy_transposed.csv treatment FALSE
#Script to run Rscripts that generate binned PCA plots

#Plot merged data binned PCA
Rscript geneCounts_PCA_eigens_binned.r $1 $2 $3
#Rename produced plot
outFile=$(echo $1 | sed 's/\.csv//')
mv Rplots.pdf $outFile"_bin"$2"_scale"$3".pdf"