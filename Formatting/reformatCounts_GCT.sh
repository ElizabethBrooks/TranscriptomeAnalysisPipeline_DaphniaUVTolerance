#!/bin/bash
#Script to re-format merged count tables to Gene Cluster Text file format (*.gct)
#Usage: bash reformatCounts_GCT.sh countsFile
#Usage Ex: bash reformatCounts_GCT.sh ~/PfrenderLab/PA42_v4.1/geneCounts_cleaned_PA42_v4.1.csv
#Usage Ex: bash reformatCounts_GCT.sh ~/PfrenderLab/WGCNA_PA42_v4.1/normalizedCountsInter_PA42_v4.1.csv

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
#Change delimiter and clean up input csv file
cat "$inCountsFile" | sed 's/"",/gene,/g' | sed 's/"//g' | sed 's/,/\t/g' > tmpinCountsFile.txt
#Retrieve number of rows
numRows=$(($(wc -l < tmpinCountsFile.txt | sed 's/ //g')-1))
#Retrieve number of samples
numCols=$(($(head -n1 tmpinCountsFile.txt | awk '{print NF}')-1))
#Set output file name
outFile="$inputsPath"/"$countsFile"_reformatted.gct
#Output headers for GCT formatting
echo "#1.2" > tmpHeader.gct
echo -e "$numRows \t $numCols" >> tmpHeader.gct
#Create temporary file with added empty second column for the 'description' field
cut -f1 tmpinCountsFile.txt > tmpData1.gct
sed -i '.bak' -e "s/$/\tNA/" tmpData1.gct 
cut -f2- tmpinCountsFile.txt > tmpData2.gct
paste tmpData1.gct tmpData2.gct > tmpData3.gct
sed -i '.bak' 's/gene\tNA/Name\tDescription/g' tmpData3.gct
#Append header to reformatted counts table
cat tmpHeader.gct tmpData3.gct > "$outFile"
#Print a script completion confirmation message
echo "$countsFile has been reformatted!"
#Clean up
rm tmp*.txt
rm tmp*.gct
rm *.bak