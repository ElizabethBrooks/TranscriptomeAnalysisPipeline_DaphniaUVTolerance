#!/bin/bash
#Usage: bash retrieveGeneCounts.sh -addFlag GeneCountAnalysisFilePath geneID
#Usage Ex: bash retrieveGeneCounts.sh -no GeneCountAnalysis_subset_run1/geneCounts_merged_counted_htseqTophat2_run1_subset_cleaned.csv gene10997
#Script to retrieve gene counts for a specififed gene ID
# and add to a csv file if selected
#Retrieve gene count analysis inputs absolute path
inputsPath=$(grep "geneTableAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/geneTableAnalysis://g")
#Retrieve gene count analysis outputs absolute path
outputsPath=$(grep "statistics:" ../InputData/outputPaths.txt | tr -d " " | sed "s/statistics://g")
#Retrieve gene counts for the input ID 
counts=$(grep "$3" "$inputsPath"/GeneCounts_Formatted/"$2")
#Write counts to file if selected
if [[ "$1" == "-yes" || "$1" == "-y" || "$1" == "-Yes" || "$1" == "-Y" || "$1" == "-YES" ]]; then
	outFolder=$(dirname "$2")
	mkdir "$outputsPath"/GeneCounts_Stats/"$outFolder"
	echo "$counts" >> "$outputsPath"/GeneCounts_Stats/"$outFolder"/selectedCounts.csv
elif [[ "$1" == "-no" || "$1" == "-n" || "$1" == "-No" || "$1" == "-N" || "$1" == "-NO" ]]; then
	echo "$counts"
else
	echo "$counts"
fi