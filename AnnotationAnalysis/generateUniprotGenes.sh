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
inFile=$(basename "$1")
inputSet="$inputPath"/"$inFile"
#Retrieve analysis outputs path
outDir="$inputPath"

#Move to current outputs folder
cd "$outDir"

#Remove file extension
inputResults=$(basename "$inputSet" | sed 's/\.gmx//g')

#Make table of GO data for the input DE analysis results
tail -n+3 "$inputSet" > tmpIDs.txt

#Rertieve GO annotations for DE analysis results
head -1 "$ontologyPath" > tmpGO.txt
while IFS= read -r line; do grep $line "$ontologyPath" >> tmpGO.txt ; done < tmpIDs.txt

#Retieve sprot IDs from selected annotations
cut -d"," -f1 tmpGO.txt > tmpGeneIDs.gmx

#Remove duplicate gene ID matches
head -2 "$inputSet" > geneIDs_noDups_"$inputResults".gmx
tail -n+2 tmpGeneIDs.gmx | awk '!visited[$0]++' >> geneIDs_noDups_"$inputResults".gmx

#Clean up
rm tmp*