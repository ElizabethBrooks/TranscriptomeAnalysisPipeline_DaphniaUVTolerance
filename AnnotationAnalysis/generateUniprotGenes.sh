#!/bin/bash
#Script to run Rscripts that retieve annotations for DE analysis results
#Usage: bash generateUniprotGenes.sh analysisResultsFile
#Usage Ex: bash generateUniprotGenes.sh geneset_dnaRepair.gmx

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi

#Retrieve gene ontology data path
ontologyPath=$(grep "geneOntology:" ../InputData/inputPaths.txt | tr -d " " | sed "s/geneOntology://g")
#Retrieve inputs path
inputPath=$(grep "geneSets:" ../InputData/inputPaths.txt | tr -d " " | sed "s/geneSets://g")
inputSet="$inputPath"/"$1"
#Retrieve analysis outputs path
outDir=$(grep "GSEA:" ../InputData/outputPaths.txt | tr -d " " | sed "s/GSEA://g")

#Move to current outputs folder
cd "$outDir"

#Remove file extension
inputResults=$(basename "$inputSet" | sed 's/\.gmx//g')

#Make table of GO data for the input DE analysis results
tail -n+3 "$inputSet" > tmpIDs.txt

#Rertieve GO annotations for DE analysis results
while IFS= read -r line; do grep $line "$ontologyPath" >> tmp_"$inputResults".txt ; done < tmpIDs.txt

#Retieve sprot IDs from selected annotations
cut -d"," -f1 tmp_"$inputResults".txt > geneIDs_"$inputResults".gmx

#Remove duplicate gene ID matches
awk '!visited[$0]++' geneIDs_"$inputResults".gmx > geneIDs_noDups_"$inputResults".gmx

#Clean up
rm tmp*.txt