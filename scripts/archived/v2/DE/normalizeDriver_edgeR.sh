#!/bin/bash

# script to run Rscripts that perform DE analysis of gene count tables using glm in edgeR
# usage: bash normalizeDriver_edgeR.sh
# usage Ex: bash normalizeDriver_edgeR.sh Genotypes
# usage Ex: bash normalizeDriver_edgeR.sh Tolerance

# load module for R
#module load bio

# check for input arguments
if [ $# -eq 0 ]; then
   	echo "ERROR: No argument(s) supplied... exiting"
   	exit 1
fi

# retrieve analysis type
analysisType=$1

#Create directory for output files
outDir=$(grep "DEAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/DEAnalysis://g")

# name output directory for the analysis type
outDir=$outDir"/"$analysisType

# retrieve current directory
curDir=$(pwd)

# retrieve experimental design data path
cd $(dirname "../InputData/expDesign_Olympics"$analysisType".csv")
designPath=$(pwd)
designPath=$designPath"/expDesign_Olympics"$analysisType".csv"

# move back to origional directory
cd $curDir

#Retrieve analysis inputs path
inFile=$(grep "cleanedGeneCounts:" ../InputData/outputPaths.txt | tr -d " " | sed "s/cleanedGeneCounts://g")
inFile=$inFile"/cleaned.csv"

#Perform normalization of raw counts
Rscript "normalizeLogCounts"$analysisType"_edgeR.r" "$outDir" "$inFile" 1 24 "$designPath"

# clean produced tables
file=$outDir"/normalizedLogCounts.csv"

# fix header
tail -n+2 "$file" > $outDir"/tmpTail.csv"
head -1 "$file" | sed -e 's/^/gene,/' > $outDir"/tmpHeader.csv"

# update table
cat $outDir"/tmpHeader.csv" > "$file"
cat $outDir"/tmpTail.csv" >> "$file"

# clean up
rm $outDir"/tmpHeader.csv"
rm $outDir"/tmpTail.csv"
