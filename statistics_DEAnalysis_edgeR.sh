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
#Retrieve statistics outputs absolute path
outputsPath=$(grep "statistics:" InputData/outputPaths.txt | tr -d " " | sed "s/statistics://g")
#Move to outputs directory
cd "$outputsPath"
#Perform DE analysis using edgeR
Rscript statistics_edgeR.r "$1" $2 $3
#Make directory for output stats files
#mkdir ../AlignmentStats_Analysis
#Move produce stats file
outFile=$(basename "$1")
mv stats_tmpOut.csv ../AlignmentStats_Analysis/alignmentStats_cols"$2"to"$3"_"$outFile"