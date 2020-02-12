#!/bin/bash
#Script to run Rscripts that perform DE analysis of gene count tables
#Usage: bash statistics_filtered_DEAnalysis_edgeR.sh countsFile startColPos endColPos
#Usage Ex: bash statistics_filtered_DEAnalysis_edgeR.sh GeneCountAnalysis_subset_run1/geneCounts_merged_counted_htseqTophat2_run1_subset_cleaned.csv 1 6

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve statistics outputs absolute path
outputsPath=$(grep "statistics:" ../InputData/outputPaths.txt | tr -d " " | sed "s/statistics://g")
outputCounts="$outputsPath"/geneCountStats_cols"$2"to"$3"
#Retrieve analysis inputs path
inputsPath=$(grep "geneTableAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/geneTableAnalysis://g")
outFile=$(basename "$inputsPath"/"$1" | sed 's/\.csv//g')
#Perform DE analysis using edgeR
Rscript DEStatistics_filtered_edgeR.r "$inputsPath"/"$1" $2 $3
#Rename and move produced filtered table
mv tmpOut.csv "$outputCounts"_"$outFile".csv
#Rename and move produced plot
mv Rplots.pdf "$outputCounts"_"$outFile".pdf