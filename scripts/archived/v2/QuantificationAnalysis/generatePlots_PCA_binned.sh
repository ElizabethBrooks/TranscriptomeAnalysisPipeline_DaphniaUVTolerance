#!/bin/bash
#Usage: bash generatePlots_PCA_binned.sh GeneCountAnalysisFilePath bin
#Usage Ex: bash generatePlots_PCA_binned.sh GeneCountAnalysis_subset_run1/geneCounts_merged_counted_htseqTophat2_run1_subset_final_merged_counts_annotatedTreatment_transposed.csv treatment
#Script to run Rscripts that generate binned PCA plots
#Retrieve gene count analysis inputs absolute path
inputsPath=$(grep "geneTableAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/geneTableAnalysis://g")
#Retrieve gene count analysis outputs absolute path
outputsPath=$(grep "statistics:" ../InputData/outputPaths.txt | tr -d " " | sed "s/statistics://g")
#Plot merged data binned PCA
Rscript geneCounts_PCA_binned.r "$inputsPath"/GeneCounts_Formatted/"$1" $2
#Rename and move produced plot
outFolder=$(dirname "$1")
mkdir "$outputsPath"/GeneCounts_Stats/"$outFolder"
outFile=$(echo "$1" | sed 's/\.csv//')
mv Rplots.pdf "$outputsPath"/GeneCounts_Stats/"$outFile"_"$2"Bin.pdf