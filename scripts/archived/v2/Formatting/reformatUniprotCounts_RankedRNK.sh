#!/bin/bash
#Script to re-format merged count tables to a ranked list file format (*.rnk)
#Usage: bash reformatUniprotCounts_RankedGCT.sh countsFile
#Usage Ex: bash reformatUniprotCounts_RankedRNK.sh /Users/bamflappy/PfrenderLab/dMelUV/WGCNA_PA42_v4.1/filteredCountFiles/normalizedCountsInter_PA42_v4.1_uniprot_pAdjust_cleaned.csv

#Check for input argument of file name
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi

#Retrieve input filename and path
inCountsFile="$1"
inputsPath=$(dirname $1)
countsFile=$(basename "$inCountsFile" | sed 's/\.csv//g')
#Set output file name
outFile="$inputsPath"/"$countsFile"_noDups.rnk
outFileDups="$inputsPath"/"$countsFile"_dups.rnk
#Keep only the gene ID and weight columns
tail -n+2 $inCountsFile | cut -d"," -f1,2 > tmp.csv
#Add header to output files
head -1 $inCountsFile | cut -d"," -f1,2 > $outFileDups
head -1 $inCountsFile | cut -d"," -f1,2 > $outFile
#Filter the ranked list to remove duplicate IDs
#Keep the first ID occurance since the list is sorted by ascending adjusted p-value
while read -r line; do
	gTag=$(echo "$line" | cut -d"," -f1)
	gTag=$gTag","
	if(grep -q "$gTag" $outFile); then
		echo "$line" >> $outFileDups
	else
		echo "$line" >> $outFile
	fi
done < tmp.csv
#Print a script completion confirmation message
echo "$countsFile has been reformatted!"
#Clean up
rm tmp.csv