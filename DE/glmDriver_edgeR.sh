#!/bin/bash
#Script to run Rscripts that perform DE analysis of gene count tables using glm in edgeR
#Usage: bash glmDriver_edgeR.sh referenceGenome analysisType
#Usage Ex: bash glmDriver_edgeR.sh KAP4 Tolerance
#Usage Ex: bash glmDriver_edgeR.sh KAP4 Genotypes

#Load module for R
#module load bio

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi

#Retrieve experimental design data path
designPath="../InputData/expDesign_Olympics"$2".csv"
#Retrieve analysis inputs path
inFile=$(grep "geneCounts:" ../InputData/inputPaths_"$1".txt | tr -d " " | sed "s/geneCounts://g")
echo "$inFile"
#Set FDR cut off
#fdrCut=0.10

#Create directory for output files
outDir="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/"$1"/DE"$2
mkdir $outDir
#Check if the folder already exists
#if [ $? -ne 0 ]; then
#	echo "The $outDir directory already exsists... please remove before proceeding."
#	exit 1
#fi

#Perform DE analysis using glmLRT in edgeR and output analysis results to a txt file
Rscript glmQLF_Olympics"$2"_edgeR.r "$inFile" 1 24 "$designPath" > "$outDir"/glmQLF_analysisResults.txt

#Move produced tables
for f in *.csv; do
	file="$f"
	sed -i '' 's/"//g' "$file"
	#Fix header
	tail -n+2 "$file" > tmpTail.csv
	head -1 "$file" | sed -e 's/^/gene,/' > tmpHeader.csv
	#Update table
	cat tmpHeader.csv > "$file"
	cat tmpTail.csv >> "$file"
	rm tmp*.csv
	#Move updated table
	mv "$file" "$outDir"
done
#Move produced plots
for p in *.jpg; do mv $p "$outDir"; done
