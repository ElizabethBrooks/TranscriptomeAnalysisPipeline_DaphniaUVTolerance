#!/bin/bash
#Usage: bash generatePlots_kMeans_binned.sh mergedCounts_file_transposed.csv bin k
#Usage Ex: bash generatePlots_kMeans_binned.sh mergedCounts_legacy_transposed.csv treatment 3
#Script to run Rscripts that generate binned kMeans plots

#Plot merged data binned kMeans clustering
Rscript geneCounts_kMeans_binned.r $1 $2 $3
#Rename produced plot
outFile=$(echo $1 | sed 's/\.csv//')
mv Rplots.pdf $outFile"_k"$2"_bin"$3".pdf"