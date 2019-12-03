#!/bin/bash
#Usage: bash generatePlots_PCA_binned.sh GeneCountAnalysisFilePath bin
#Usage Ex: bash generatePlots_PCA_binned.sh GeneCountAnalysis_subset_run1/geneCounts_merged_counted_htseqTophat2_run1_subset_transposed.csv treatment
#Script to run Rscripts that generate binned PCA plots
#Retrieve gene count analysis outputs absolute path
outputsPath=$(grep "geneTableAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/geneTableAnalysis://g")
#Plot merged data binned PCA
Rscript geneCounts_PCA_binned.r "$outputsPath"/"$1" $2
#Rename produced plot
outFile=$(echo "$1" | sed 's/\.csv//')
mv Rplots.pdf "$outputsPath"/"$outFile"_bin_"$2".pdf