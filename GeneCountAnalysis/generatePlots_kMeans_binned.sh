#!/bin/bash
#Usage: bash generatePlots_kMeans_binned.sh GeneCountAnalysisFilePath bin k
#Usage Ex: bash generatePlots_kMeans_binned.sh GeneCountAnalysis_subset_run1/geneCounts_merged_counted_htseqTophat2_run1_subset_transposed.csv treatment 3
#Script to run Rscripts that generate binned kMeans plots
#Retrieve gene count analysis outputs absolute path
outputsPath=$(grep "geneTableAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/geneTableAnalysis://g")
#Plot merged data binned kMeans clustering
Rscript geneCounts_kMeans_binned.r "$outputsPath"/"$1" $2 $3
#Rename produced plot
outFile=$(echo "$1" | sed 's/\.csv//')
mv Rplots.pdf "$outputsPath"/"$outFile"_k"$2"_bin_"$3".pdf