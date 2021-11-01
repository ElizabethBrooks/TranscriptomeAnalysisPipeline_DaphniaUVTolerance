#!/bin/bash
#Script to re-format merged count tables to Gene Cluster Text file format (*.gct)
#Usage: bash reformatUniprotCounts_RankedGCT.sh countsFile
#Usage Ex: bash reformatUniprotCounts_RankedGCT.sh /Users/bamflappy/PfrenderLab/dMelUV/WGCNA_PA42_v4.1/filteredCountFiles/normalizedCountsInter_PA42_v4.1_uniprot_pAdjust_cleaned.csv

#Check for input argument of file name
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Initialize values
numRows=0
numCols=0
#Retrieve input filename and path
inCountsFile="$1"
inputsPath=$(dirname $1)
countsFile=$(basename "$inCountsFile" | sed 's/\.csv//g')
#Remove the column of adjusted p-values
cut -d"," -f1 $inCountsFile > tmpCols1.csv
cut -d"," -f3- $inCountsFile > tmpCols2.csv
paste -d"," tmpCols1.csv tmpCols2.csv > tmpinFile.csv
#Change delimiter and clean up input csv file
cat tmpinFile.csv | sed 's/"//g' | sed 's/,/\t/g' > tmpinCountsFile.txt
#Retrieve number of rows
numRows=$(($(wc -l < tmpinCountsFile.txt | sed 's/ //g')-1))
#Retrieve number of samples
numCols=$(($(head -n1 tmpinCountsFile.txt | awk '{print NF}')-1))
#Set output file name
outFile="$inputsPath"/"$countsFile"_reformatted.gct
#Output headers for GCT formatting
echo "#1.2" > tmpHeader.txt
echo -e "$numRows \t $numCols" >> tmpHeader.txt
#Update header
sed -i '.bak' 's/sprot\tgene/Name\tDescription/g' tmpinCountsFile.txt
#Append header to reformatted counts table
cat tmpHeader.txt tmpinCountsFile.txt > "$outFile"
#Print a script completion confirmation message
echo "$countsFile has been reformatted!"
#Clean up
rm tmp*
rm *.bak