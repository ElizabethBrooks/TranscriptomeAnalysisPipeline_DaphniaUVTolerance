#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N stats_edgeR_jobOutput
#Script to run Rscripts that perform DE analysis of gene count tables
#Usage: bash statistics_DEAnalysis_edgeR.r countsFile startColPos endColPos
#Usage Ex: bash statistics_DEAnalysis_edgeR.r ../GeneCounts_Merged/merged_counts_legacy_cleaned.csv 1 6
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve outputs absolute path
outputsFile="InputData/outputsPath.txt"
outputsPath=$(head -n 1 $outputsFile)
#Move to outputs directory
cd "$outputsPath"
#Perform DE analysis using edgeR
Rscript statistics_filtered_edgeR.r "$1" $2 $3
#Make directory for output stats files
#mkdir ../AlignmentStats_Analysis
#Rename and move produced plot
outFile=$(basename "$1")
outFile=$(echo "$outFile" | sed 's/\.csv//')
mv Rplots.pdf ../AlignmentStats_Analysis/alignmentStats_filtered_cols"$2"to"$3"_"$outFile".pdf