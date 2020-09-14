#!/bin/bash
#Script to run Rscripts that perform DE analysis of gene count tables using glm in edgeR
#Usage: bash glmDriver_edgeR.sh analysisType countsFolder
#Usage Ex: bash glmDriver_edgeR.sh QLF genome_sortedName_samtoolsHisat2_run1_counted_htseq_run1
#Usage Ex: bash glmDriver_edgeR.sh LRT genome_sortedName_samtoolsHisat2_run1_counted_htseq_run1

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi

#Retrieve gene ontology data path
ontologyPath=$(grep "geneOntology:" ../InputData/inputPaths.txt | tr -d " " | sed "s/geneOntology://g")
#Retrieve analysis inputs path
inputsPath=$(grep "DEAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/DEAnalysis://g")
inFile="$inputsPath"/"$2"/cleaned.csv

#Create directory for output files
outDir="$inputsPath"/"$2"/DEAnalyzed
mkdir $outDir

#Determine analysis method
if [[ "$1" == "LRT" ]]; then
	#Perform DE analysis using glmLRT in edgeR and output analysis results to a txt file
	Rscript glmLRT_edgeR.r "$inFile" 1 24 "$ontologyPath" > "$outDir"/glmLRT_analysisResults.txt
elif [[ "$1" == "QLF" ]]; then
	#Perform DE analysis using glmLRT in edgeR and output analysis results to a txt file
	Rscript glmQLF_edgeR.r "$inFile" 1 24 "$ontologyPath" > "$outDir"/glmQLF_analysisResults.txt
else
	echo "Invalid analysis type entered... exiting!"
	exit 1
fi

#Rename and move produced tables
for f in *.csv; do mv $f "$outDir"; done
#Rename and move produced plots
for p in *.jpg; do mv $p "$outDir"; done
