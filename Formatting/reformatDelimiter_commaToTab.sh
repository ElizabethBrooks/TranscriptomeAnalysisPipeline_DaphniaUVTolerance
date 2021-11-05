#!/bin/bash
#Script to re-format merged count tables to a ranked list file format (*.rnk)
#Usage: bash reformatDelimiter_commaToTab.sh countsFile
#Usage Ex: bash reformatDelimiter_commaToTab.sh /Users/bamflappy/PfrenderLab/dMelUV/WGCNA_PA42_v4.1/filteredCountFiles/normalizedCountsInter_PA42_v4.1_uniprot_pAdjust_cleaned_noDups_reRanked.rnk

#Check for input argument of file name
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi

#Retrieve input filename and path
inCountsFile="$1"
inputsPath=$(dirname $1)
countsFile=$(basename "$inCountsFile" | sed 's/\.rnk//g' | sed 's/\.csv//g')
#Set output file name
outFile="$inputsPath"/"$countsFile"_tabDelimited.rnk
#Replace commas with tabs
cat $inCountsFile > $outFile
sed -i '.bak' "s/,/\t/g" $outFile
#Print a script completion confirmation message
echo "$countsFile has been reformatted!"
#Clean up
rm "$inputsPath"/*.bak