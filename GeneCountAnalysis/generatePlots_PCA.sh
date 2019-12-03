#!/bin/bash
#Usage: bash generatePlots_PCA.sh GeneCountAnalysisFilePath
#Usage Ex: bash generatePlots_PCA.sh GeneCountAnalysis_subset_run1/geneCounts_merged_counted_htseqTophat2_run1_subset_transposed.csv
#Script to run Rscripts that generate PCA plots
#Retrieve gene count analysis outputs absolute path
outputsPath=$(grep "geneTableAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/geneTableAnalysis://g")
#Plot merged data PCA
Rscript geneCounts_PCA.r "$outputsPath"/"$1"
#Rename produced plot
outFile=$(echo "$1" | sed 's/\.csv//')
mv Rplots.pdf "$outputsPath"/"$outFile".pdf