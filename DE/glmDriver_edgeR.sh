#!/bin/bash
# script to run Rscripts that perform DE analysis of gene count tables using glm in edgeR
# usage: bash glmDriver_edgeR.sh analysisType
# default usage Ex: bash glmDriver_edgeR.sh Tolerance
# usage Ex: bash glmDriver_edgeR.sh Genotypes
## potential Usage: bash glmDriver_edgeR.sh analysisType referenceGenome

#Load module for R
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
mkdir $outDir

# name output directory for the analysis type
outDir=$outDir"/"$analysisType
mkdir $outDir
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outDir directory already exsists... please remove before proceeding."
	exit 1
fi

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

#Set FDR cut off
#fdrCut=0.10

# output status message
echo "Performing $analysisType DE analysis of $inFile"

#Perform DE analysis using glmLRT in edgeR and output analysis results to a txt file
Rscript glmQLF_Olympics"$analysisType"_edgeR.R "$outDir" "$inFile" 1 24 "$designPath" > "$outDir"/glmQLF_analysisResults.txt

# clean produced tables
for f in "$outDir"/*.csv; do
	file="$f"
	#sed -i 's/"//g' "$file"
	#Fix header
	tail -n+2 "$file" > tmpTail.csv
	head -1 "$file" | sed -e 's/^/gene,/' > tmpHeader.csv
	#Update table
	cat tmpHeader.csv > "$file"
	cat tmpTail.csv >> "$file"
	rm tmp*.csv
done
