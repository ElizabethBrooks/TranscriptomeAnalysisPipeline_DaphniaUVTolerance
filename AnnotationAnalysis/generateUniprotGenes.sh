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
head -1 "$ontologyPath" > tmpList.txt
while IFS= read -r line
do
	if grep -q $line "$ontologyPath"; then #Match
		go=$(grep $line "$ontologyPath")
		echo $go >> tmpGO.txt
		echo $go >> tmpList.txt
	else #Unique
		echo $line",NA" >> tmpList.txt
	fi
done < tmpIDs.txt

#Retieve gene IDs from selected annotations
cut -d"," -f1 tmpGO.txt > tmpGOGeneIDs.gmx
cut -d"," -f1 tmpList.txt >> geneIDs_fullList_"$inputResults".txt

#Remove duplicate gene ID matches
head -2 "$inputSet" > geneIDs_noDups_"$inputResults".gmx
tail -n+2 tmpGOGeneIDs.gmx | awk '!visited[$0]++' >> geneIDs_noDups_"$inputResults".gmx

#Clean up
rm tmp*