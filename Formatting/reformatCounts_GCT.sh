#!/bin/bash
#Script to re-format merged count tables to Gene Cluster Text file format (*.gct)
#Usage: bash reformatCounts_GCT.sh countsFolder countsFile
#Usage Ex: bash reformatCounts_GCT.sh genome_sortedName_samtoolsHisat2_run2_counted_htseq_run1Analysis cleaned.csv

#Check for input argument of file name
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve analysis inputs path
inputsPath=$(grep "DEAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/DEAnalysis://g")
#Initialize values
numRows=0
numCols=0
#Retrieve input filename
inFile="$inputsPath"/"$1"/"$2"
countsFile=$(basename "$inFile" | sed 's/\.csv//g')
#Change delimiter for csv files
sed 's/,/\t/g' "$inFile" > tmpInFile.txt
#Retrieve number of rows
wc -l tmpInFile.txt > tmpNumRows.txt
numRows=$(cut -d ' ' -f1 tmpNumRows.txt)
#Retrieve number of samples
numCols=$(($(head -n1 tmpInFile.txt | awk '{print NF}')-1))
#Set output file name
outFile="$inFile"/"$countsFile"_reformatted.gct
#Output headers for GCT formatting
echo "#1.2" > tmpHeader.gct
echo -e "$numRows \t $numCols" >> tmpHeader.gct
#Create temporary file with added empty second column for the 'description' field
cut -f1 tmpInFile.txt > tmpData1.gct
sed -e "s/$/\tNA/" -i tmpData1.gct
cut -f2- tmpInFile.txt > tmpData2.gct
paste tmpData1.gct tmpData2.gct > tmpData3.gct
sed -i "s/gene\tNA/Name\tDescription/g" tmpData3.gct
#Append header to reformatted counts table
cat tmpHeader.gct tmpData3.gct > "$outFile"
#Print a script completion confirmation message
echo "Gene counts file has been reformatted!"
#Clean up
rm tmp*.txt
rm tmp*.gct