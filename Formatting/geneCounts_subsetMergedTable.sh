#!/bin/bash
#Usage: bash geneCounts_subsetMergedTable.sh countedGenesFile sampleName
#Usage Ex: bash geneCounts_subsetMergedTable.sh geneCounts_merged_countedCoordinate_htseqHisat2_run1_fullset_reformatted.gct R2
#Script to select a specified subset of the gene count data

#Check for input argument of file name
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve statistics outputs absolute path
outputsPath=$(grep "geneTableAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/geneTableAnalysis://g")
#Retrieve analysis inputs path
inputsPath=$(grep "geneTableAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/geneTableAnalysis://g")
#Initialize values
colNum=0
#Retrieve input filename
inFile="$inputsPath/$1"
#Set output file name
outFile="$outputsPath"/subset"$2"_"$1"
#Set sample subset
subsetStart="$2"_VIS_Pool1
#Check file type
if [ ${inFile: -4} == ".gct" ]; then #GCT formatted
	#Retrieve column number of selected sample subset
	head -3 "$inFile" > tmpHeader.txt
	tail tmpHeader.txt | tr "\t" "\n" | grep -n "$subsetStart" > tmpColNum.txt
	colNum=$(cut -d ':' -f1 tmpColNum.txt)
	#Retrieve the selected subset (default of 3 replicates for 2 treatments), 
	# including the first column with gene IDs and second with description
	cut -f 1,2,$colNum,$(($colNum+1)),$(($colNum+2)),$(($colNum+3)),$(($colNum+4)),$(($colNum+5)) "$inFile" > "$outFile"
else #TXT and CSV formatted
	#Retrieve column number of selected sample subset
	head -1 "$inFile" | tr "\t" "\n" | grep -n "$subsetStart" > tmpColNum.txt
	colNum=$(cut -d ':' -f1 tmpColNum.txt)
	#Retrieve the selected subset (default of 3 replicates for 2 treatments), 
	# including the first column with gene IDs
	cut -f 1,$colNum,$(($colNum+1)),$(($colNum+2)),$(($colNum+3)),$(($colNum+4)),$(($colNum+5)) "$inFile" > "$outFile"
fi
#Print a script completion confirmation message
echo "Selected $2 subset has been generated!"
#Clean up
rm tmp*.txt