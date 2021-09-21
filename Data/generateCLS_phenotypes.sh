#!/bin/bash
#Usage: bash generateCLS_phenotypes.sh
#Usage Ex: bash generateCLS_phenotypes.sh
#Script to generate a categorical class file format (*.cls)

#Retrieve statistics outputs absolute path
outputsPath=$(grep "DEAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/DEAnalysis://g")
#Initialize values
numSamples=36
numPheno=2
#Set output file name
outFile="$outputsPath"/dMelUV_phenotypes.cls
#Set header line and write to file
headerLine="$numSamples $numPheno 1"
echo $headerLine > "$outFile"
#Set class name line and write to file
classLine="# Tol NTol"
echo $classLine >> "$outFile"
#Set class label line and write to file
labelLine="Tol Tol Tol Tol Tol Tol Tol Tol Tol Tol Tol NTol NTol NTol NTol NTol NTol NTol NTol NTol NTol NTol NTol NTol NTol NTol NTol NTol NTol Tol Tol Tol Tol Tol Tol"
echo $labelLine >> "$outFile"
#Print a script completion confirmation message
echo "CLS phenotype file has been generated!"