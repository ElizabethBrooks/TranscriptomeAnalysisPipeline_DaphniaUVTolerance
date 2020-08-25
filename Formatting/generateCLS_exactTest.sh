#!/bin/bash
#Usage: bash generateCLS_exactTest.sh sampleName
#Usage Ex: bash generateCLS_exactTest.sh R2
#Script to generate a categorical class file format (*.cls)

#Check for input argument of file name
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve statistics outputs absolute path
outputsPath=$(grep "geneCountAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/geneCountAnalysis://g")
#Initialize values
numSamples=6
numTreatments=2
#Set output file name
outFile="$outputsPath"/GeneCounts_Formatted/CLS/"$1"_exactTest.cls
#Set header line and write to file
headerLine="$numSamples $numTreatments 1"
echo $headerLine > "$outFile"
#Set class name line and write to file
classLine="# "$1"_VIS "$1"_UV"
echo $classLine >> "$outFile"
#Set class label line and write to file
labelLine="VIS VIS VIS UV UV UV"
echo $labelLine >> "$outFile"
#Print a script completion confirmation message
echo "$1 exact test CLS file has been generated!"