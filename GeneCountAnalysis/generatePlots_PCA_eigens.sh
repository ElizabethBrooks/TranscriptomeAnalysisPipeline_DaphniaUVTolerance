#!/bin/bash
#Usage: bash generatePlots_PCA_eigens.sh mergedCounts_file_transposed.csv
#Usage Ex: bash generatePlots_PCA_eigens.sh mergedCounts_legacy_transposed.csv
#Script to run Rscripts that generate binned PCA plots
#Retrieve gene count analysis outputs absolute path
outputsPath=$(grep "geneCountAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/geneCountAnalysis://g")
#Move to outputs directory
cd "$outputsPath"
#Create directory for gene count analysis
outputAnalysis=GeneCountAnalysis
mkdir "$outputAnalysis"
#Plot merged data binned PCA
Rscript geneCounts_PCA_eigens.r "$1"
#Rename produced plot
outFile=$(echo "$1" | sed 's/\.csv//')
mv Rplots.pdf "$outputAnalysis"/"$outFile".pdf