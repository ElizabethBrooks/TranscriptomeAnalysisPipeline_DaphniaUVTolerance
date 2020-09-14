#!/bin/bash
#Script to run Rscripts that retieve annotations for DE analysis results
#Usage: bash geneAnnotations_edgeR.sh analysisResultsFile
#Usage Ex: bash geneAnnotations_edgeR.sh glmQLF_2WayANOVA_topTags_filtered.csv

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi

#Retrieve gene ontology data path
ontologyPath=$(grep "geneOntology:" ../InputData/inputPaths.txt | tr -d " " | sed "s/geneOntology://g")
#Retrieve analysis inputs path
outDir=$(grep "DEAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/DEAnalysis://g")

#Move to current outputs folder
cd "$outDir"

#Remove file extension
inputResults=$(basename "$1" | sed 's/\.csv//g')

#Retrieve all GO IDs
cut -d"," -f7 "$ontologyPath" | sed "s/\"//g" | cut -d"^" -f1 | sed -n '/\./!p' > sprotIDs_GO.txt

#Make table of GO data for the input DE analysis results
tail -n+2 "$1" | cut -f1 -d ',' | sed "s/\"//g" > tmpIDs.txt

#Rertieve GO annotations for DE analysis results
head -1 "$ontologyPath" > GO_"$inputResults".csv
while IFS= read -r line; do grep $line"," "$ontologyPath" >> GO_"$inputResults".csv ; done < tmpIDs.txt

#Retieve sprot IDs from selected annotations
cut -d"," -f7 GO_"$inputResults".csv | sed "s/\"//g" | cut -d"^" -f1 | sed -n '/\./!p' > sprotIDs_"$inputResults".txt

#Clean up
rm tmpIDs.txt