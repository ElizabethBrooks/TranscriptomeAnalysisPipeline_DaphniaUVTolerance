#!/bin/bash
#Script to run Rscripts that perform DE analysis of gene count tables using glm in edgeR
#Usage: bash glmDriver_pvalues_edgeR.sh

#Load module for R
#module load bio

#Retrieve experimental design data path
designPath="../InputData/expDesign_olympics.csv"
#Retrieve analysis inputs path
inputsPath=$(grep "DEAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/DEAnalysis://g")
inFile="$inputsPath"/cleaned.csv

#Create directory for output files
outDir="$inputsPath"/glm"$1"Analysis
mkdir $outDir

#Perform DE analysis using glmLRT in edgeR 
Rscript glmQLF_edgeR_pvalues.r "$inFile" 1 24 "$designPath"

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
#Move produced plots
for p in *.jpg; do mv $p "$outDir"; done
