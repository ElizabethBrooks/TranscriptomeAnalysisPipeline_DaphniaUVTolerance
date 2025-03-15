#!/bin/bash
#Usage: bash generatePlots_PCA_eigens.sh GeneCountAnalysisFilePath
#Usage Ex: bash generatePlots_PCA_eigens.sh GeneCountAnalysis_subset_run1/mergedCounts_legacy_transposed.csv
#Script to run Rscripts that generate binned PCA plots
#Retrieve gene count analysis inputs absolute path
inputsPath=$(grep "geneTableAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/geneTableAnalysis://g")
#Retrieve gene count analysis outputs absolute path
outputsPath=$(grep "statistics:" ../InputData/outputPaths.txt | tr -d " " | sed "s/statistics://g")
#Plot merged data binned PCA
Rscript geneCounts_PCA_eigens.r "$inputsPath"/GeneCounts_Formatted/"$1"
#Rename and move produced plot
outFolder=$(dirname "$1")
mkdir "$outputsPath"/GeneCounts_Stats/"$outFolder"
outFile=$(echo "$1" | sed 's/\.csv//')
mv Rplots.pdf "$outputsPath"/GeneCounts_Stats/"$outFile".pdf