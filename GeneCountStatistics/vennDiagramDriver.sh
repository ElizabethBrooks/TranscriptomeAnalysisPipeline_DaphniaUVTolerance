#!/bin/bash
#Script to generate venn diagrams from edgeR stats
#Usage: bash 
#Usage Ex: bash 

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve statistics outputs absolute path
outputsPath=$(grep "statistics:" ../InputData/outputPaths.txt | tr -d " " | sed "s/statistics://g")
#Retrieve analysis inputs path
inputsPath=$(grep "geneTableAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/geneTableAnalysis://g")
outFile=$(basename "$inputsPath"/"$1" | sed 's/\.csv//g')
#Create directory for output files
outDir="$outputsPath"/"$outFile"
outputStats="$outDir"/geneCountStats_cols"$2"to"$3"_"$outFile"
mkdir "$outDir"
mkdir "$outputStats"

#TO DO
#Prepare for analysis
Rscript geneSets_vennDiagram.r "$1"
Rscript topTags_vennDiagram.r "$1"
#Move output JPGs