#!/bin/bash
#Usage: bash generateCLS_olympicsTreatment.sh
#Usage Ex: bash generateCLS_olympicsTreatment.sh
#Script to generate a categorical class file format (*.cls)

#Retrieve statistics outputs absolute path
outputsPath=$(grep "DEAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/DEAnalysis://g")
#Initialize values
numSamples=24
numPheno=2
#Set output file name
outFile="$outputsPath"/dMelUV_olympicsTreatment.cls
#Set header line and write to file
headerLine="$numSamples $numPheno 1"
echo $headerLine > "$outFile"
#Set class name line and write to file
classLine="# VIS UV"
echo $classLine >> "$outFile"
#Set class label line and write to file
labelLine="VIS VIS VIS UV UV UV VIS VIS VIS UV UV UV VIS VIS VIS UV UV UV VIS VIS VIS UV UV UV"
echo $labelLine >> "$outFile"
#Print a script completion confirmation message
echo "CLS phenotype file has been generated!"