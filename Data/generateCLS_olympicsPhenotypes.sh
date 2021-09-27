#!/bin/bash
#Usage: bash generateCLS_olympicsPhenotypes.sh
#Usage Ex: bash generateCLS_olympicsPhenotypes.sh
#Script to generate a categorical class file format (*.cls)

#Retrieve statistics outputs absolute path
outputsPath=$(grep "DEAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/DEAnalysis://g")
#Initialize values
numSamples=24
numPheno=2
#Set output file name
outFile="$outputsPath"/dMelUV_olympicsPhenotypes.cls
#Set header line and write to file
headerLine="$numSamples $numPheno 1"
echo $headerLine > "$outFile"
#Set class name line and write to file
classLine="# VIS NTol"
echo $classLine >> "$outFile"
#Set class label line and write to file
labelLine="Tol Tol Tol Tol Tol Tol Tol Tol Tol Tol Tol Tol NTol NTol NTol NTol NTol NTol NTol NTol NTol NTol NTol NTol"
echo $labelLine >> "$outFile"
#Print a script completion confirmation message
echo "CLS phenotype file has been generated!"