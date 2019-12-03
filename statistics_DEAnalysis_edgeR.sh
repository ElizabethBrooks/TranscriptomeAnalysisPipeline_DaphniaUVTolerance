#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N stats_edgeR_jobOutput
#Script to run Rscripts that perform DE analysis of gene count tables
#Usage: bash statistics_DEAnalysis_edgeR.sh countsFile startColPos endColPos
#Usage Ex: bash statistics_DEAnalysis_edgeR.sh GeneCountAnalysis_subset_run1/geneCounts_merged_counted_htseqTophat2_run1_subset_cleaned.csv 1 6
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve statistics outputs absolute path
outputsPath=$(grep "statistics:" InputData/outputPaths.txt | tr -d " " | sed "s/statistics://g")
#Retrieve analysis inputs path
inputsPath=$(grep "geneTableAnalysis:" InputData/outputPaths.txt | tr -d " " | sed "s/geneTableAnalysis://g")
#Perform DE analysis using edgeR
Rscript statistics_edgeR.r "$inputsPath"/"$1" $2 $3
#Move produce stats file
outFile=$(basename "$inputsPath"/"$1" | sed 's/\.csv//g')
mv stats_tmpOut.csv "$outputsPath"/alignmentStats_cols"$2"to"$3"_"$outFile"