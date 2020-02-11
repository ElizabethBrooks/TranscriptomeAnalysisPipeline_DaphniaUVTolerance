#!/bin/bash
#Script to re-format merged count tables to Gene Cluster Text file format (*.gct)
#Usage: bash reformatCounts_DESeq2.sh countsFile
#Usage Ex: bash reformatCounts_DESeq2.sh GeneCountsMerged/GeneCountsAnalyzed_countedCoordinate_htseqHisat2_run1_fullset_run1/geneCounts_merged_countedCoordinate_htseqHisat2_run1_fullset.txt

#Check for input argument of file name
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve statistics outputs absolute path
outputsPath=$(grep "statistics:" ../InputData/outputPaths.txt | tr -d " " | sed "s/geneTableAnalysis://g")
#Retrieve analysis inputs path
inputsPath=$(grep "geneTableAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/geneTableAnalysis://g")
#Retrieve input filename
countsFile=$(basename "$1" | sed 's/\.txt//g')
#Set output file name
outFile="$1_reformatted.gct"
#Retrieve number of rows
numRows=$(wc -l "$1")
#Retrieve number of samples
numCols=$(head -n1 "$1" | awk '{print NF}')
#Output headers for GCT formatting
echo "#1.2" >> tmpHeader.gct
printf '$numRows \t $numCols' >> tmpHeader.gct
#Create temporary file with added empty second column for the 'description' field
cut -f1 "$1" > tmpData.gct
sed -e 's/$/\t/' -i tmpData.gct
cut -f2- "$1" >> tmpData.gct
#Append header to reformatted counts table
cat tmpHeader.gct tmpData.gct > "$outputsPath/$outFile"
#Clean up
rm tmpHeader.gct
rm tmpData.gct