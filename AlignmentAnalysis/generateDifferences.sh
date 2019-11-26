#!/bin/bash
#Usage: bash generateDifferences.sh percentsFile_legacy.csv percentsFile_newMethod.csv
#Usage Ex: bash generateDifferences.sh ../../AlignmentStats_Analysis/alignmentSummarized_legacy_subset_trimmed.csv ../../AlignmentStats_Analysis/alignmentSummarized_hisat2_trimmed.csv
#Script to run Rscripts that generate csv files of column differences
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve alignment analysis outputs absolute path
outputsPath=$(grep "alignmentAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/alignmentAnalysis://g")
#Move to outputs directory
cd "$outputsPath"
#Create directory for alignment analysis
outputAnalysis=AlignmentAnalysis
mkdir "$outputAnalysis"
#Determine what analysis method was used for the input first folder of data
if [[ "$1" == *"hisat2"*  || "$1" == *"tophat2"* ]]; then
	echo "ERROR: The "$1" file should be the legacy stats... exiting"
	exit 1
fi
#Determine what analysis method was used for the input second folder of data
if [[ "$2" == *"hisat2"*  ]]; then
	#Set analysis method for folder naming
	analysisMethod="Hisat2"
elif [[ "$2" == *"tophat2"* ]]; then
	#Set analysis method for folder naming
	analysisMethod="Tophat2"
else
	echo "ERROR: The "$2" file was not found... exiting"
	exit 1
fi
#Generate csv files of differences
Rscript generateMatrix_differences.r "$1" "$2"
#Clean up headers
sed -i "s/\"x\"/overallDifferences/g" alignmentSummarized_differences_overall.csv
sed -i "s/\"x\"/concordantDifferences/g" alignmentSummarized_differences_concordant.csv
#Move and rename produced csv files
outFile="$outputAnalysis"/alignmentSummarized_legacy"$analysisMethod"_differences
mv alignmentSummarized_differences_overall.csv "$outFile"_overall.csv
mv alignmentSummarized_differences_concordant.csv "$outFile"_concordant.csv