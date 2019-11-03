#!/bin/bash
#Usage: bash generatePlots_PCA_binned.sh mergedCounts_file_transposed.csv bin
#Usage Ex: bash generatePlots_PCA_binned.sh mergedCounts_legacy_transposed.csv treatment
#Script to run Rscripts that generate binned PCA plots

#Plot merged data binned PCA
Rscript geneCounts_PCA_binned.r $1 $2
#Rename produced plot
outFile=$(echo $1 | sed 's/\.csv//')
mv Rplots.pdf $outFile"_bin"$2".pdf"