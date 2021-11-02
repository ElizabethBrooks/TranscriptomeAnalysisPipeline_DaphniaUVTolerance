#!/bin/bash
#Script to re-format merged count tables to a ranked list file format (*.rnk)
#Usage: bash reformatUniprotCounts_AdjustedWeightRNK.sh countsFile
#Usage Ex: bash reformatUniprotCounts_AdjustedWeightRNK.sh /Users/bamflappy/PfrenderLab/dMelUV/WGCNA_PA42_v4.1/filteredCountFiles/normalizedCountsInter_PA42_v4.1_uniprot_pAdjust_cleaned_noDups.rnk

#Check for input argument of file name
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi

#Retrieve input filename and path
inCountsFile="$1"
inputsPath=$(dirname $1)
countsFile=$(basename "$inCountsFile" | sed 's/\.rnk//g')
#Set output file name
outFile="$inputsPath"/"$countsFile"_reWeighted.rnk
#Keep only the gene ID and weight columns
tail -n+2 $inCountsFile | cut -d"," -f1,2 > tmp.csv
#Add header to output file
head -1 $inCountsFile | cut -d"," -f1,2 > $outFile
#Filter the ranked list to remove duplicate IDs
#Keep the first ID occurance since the list is sorted by ascending adjusted p-value
while read -r line; do
	gTag=$(echo "$line" | cut -d"," -f1)
	pTag=$(echo "$line" | cut -d"," -f2)
	var="1"
	pTag=$(echo "$var - $pTag" | bc)
	echo "$gTag","$pTag" >> $outFile
done < tmp.csv
#Clean up decimals
sed -i '.bak' "s/,\./,0\./g" $outFile
#Print a script completion confirmation message
echo "$countsFile has been reformatted!"
#Clean up
rm tmp.csv
rm *.bak