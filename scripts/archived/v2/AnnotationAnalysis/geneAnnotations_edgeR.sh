#!/bin/bash
#Script to run Rscripts that retieve annotations for DE analysis results
#Usage: bash geneAnnotations_edgeR.sh analysisType analysisResultsFile
#Usage Ex: bash geneAnnotations_edgeR.sh QLF glmQLF_2WayANOVA_TvsN_topTags_filtered.csv
#Usage Ex: bash geneAnnotations_edgeR.sh QLF glmQLF_2WayANOVA_UVvsVIS_topTags_filtered.csv
#Usage Ex: bash geneAnnotations_edgeR.sh GSEA daphniaOlympics_leading_edge_matrix_for_results.1.gmx

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi

#Retrieve gene ontology data path
ontologyPath=$(grep "geneOntology:" ../InputData/inputPaths.txt | tr -d " " | sed "s/geneOntology://g")

#Determine input type
if [[ "$1" == "QLF" || "$1" == "LRT" ]]; then
	#Retrieve DE analysis inputs path
	inputsPath=$(grep "DEAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/DEAnalysis://g")
	outDir="$inputsPath"/glm"$1"Analysis
	#Remove file extension
	inputResults=$(echo "$2" | sed 's/\.csv//g')
elif [[ "$1" == "GSEA" ]]; then
	#Retrieve GSEA inputs path
	outDir=$(grep "geneSets:" ../InputData/inputPaths.txt | tr -d " " | sed "s/geneSets://g")
	#Remove file extension
	inputResults=$(echo "$2" | sed 's/\.gmx//g')
else
	echo "Invalid analysis type entered... exiting!"
	exit 1
fi

#Move to current outputs folder
cd "$outDir"

#Retrieve all GO IDs
cut -d"," -f7 "$ontologyPath" | sed "s/\"//g" | cut -d"^" -f1 | sed -n '/\./!p' > sprotIDs_GO.txt

#Make table of GO data for the input DE analysis results
tail -n+2 "$2" | cut -f1 -d ',' | sed "s/\"//g" > tmpIDs.txt

#Rertieve GO annotations for DE analysis results
head -1 "$ontologyPath" > GO_"$inputResults".csv
while IFS= read -r line; do grep $line"," "$ontologyPath" >> GO_"$inputResults".csv ; done < tmpIDs.txt

#Retieve sprot IDs from selected annotations
cut -d"," -f7 GO_"$inputResults".csv | sed "s/\"//g" | cut -d"^" -f1 | sed -n '/\./!p' > sprotIDs_"$inputResults".txt

#Clean up
rm tmpIDs.txt