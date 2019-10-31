#!/bin/bash
#Usage: bash generatePlots_kMeans.sh mergedCounts_file_transposed.csv k SCALE
#Usage Ex: bash generatePlots_kMeans.sh mergedCounts_legacy_transposed.csv 3 FALSE
#Script to run Rscripts that generate kMeans plots

#Plot merged data kMeans clustering
Rscript geneCounts_kMeans.r $1 $2 $3
#Rename produced plot
outFile=$(echo $1 | sed 's/\.csv//')
mv Rplots.pdf $outFile"_k"$2"_scale"$3".pdf"