#!/bin/bash
#Usage: bash generatePlots_PCA_eigens.sh GeneCountAnalysisFilePath
#Usage Ex: bash generatePlots_PCA_eigens.sh GeneCountAnalysis_subset_run1/mergedCounts_legacy_transposed.csv
#Script to run Rscripts that generate binned PCA plots
#Retrieve gene count analysis outputs absolute path
outputsPath=$(grep "geneTableAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/geneTableAnalysis://g")
#Plot merged data binned PCA
Rscript geneCounts_PCA_eigens.r "$outputsPath"/"$1"
#Rename produced plot
outFile=$(echo "$1" | sed 's/\.csv//')
mv Rplots.pdf "$outputsPath"/"$outFile".pdf