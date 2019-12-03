#!/bin/bash
#Usage: bash generatePlots_kMeans.sh GeneCountAnalysisFilePath k
#Usage Ex: bash generatePlots_kMeans.sh GeneCountAnalysis_subset_run1/geneCounts_merged_counted_htseqTophat2_run1_subset_transposed.csv 3
#Script to run Rscripts that generate kMeans plots
#Retrieve gene count analysis outputs absolute path
outputsPath=$(grep "geneTableAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/geneTableAnalysis://g")
#Plot merged data kMeans clustering
Rscript geneCounts_kMeans.r "$outputsPath"/"$1" $2
#Rename produced plot
outFile=$(echo "$1" | sed 's/\.csv//')
mv Rplots.pdf "$outputsPath"/"$outFile"_k"$2".pdf