#!/bin/bash
#Script to run Rscripts that perform DE analysis of gene count tables using glm in edgeR
#Usage: bash normalizeDriver_edgeR.sh
#Usage Ex: bash normalizeDriver_edgeR.sh

#Load module for R
#module load bio

#Retrieve gene ontology data path
designPath="../InputData/expDesign_fullSet.csv"
#Retrieve analysis inputs path
inputsPath=$(grep "DEAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/DEAnalysis://g")
inFile="$inputsPath"/cleaned.csv

#Create directory for output files
outDir="$inputsPath"

#Perform normalization of raw counts
Rscript normalize_edgeR.r "$inFile" "$designPath"

#Move produced tables
for f in *.csv; do
	file="$f"
	sed -i 's/"//g' "$file"
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
