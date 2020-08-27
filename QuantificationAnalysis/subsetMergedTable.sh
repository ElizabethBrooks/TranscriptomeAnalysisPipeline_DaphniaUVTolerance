#!/bin/bash
#Usage: bash subsetMergedTable.sh countedGenesFile sampleName
#Usage Ex: bash subsetMergedTable.sh geneCounts_merged_countedCoordinate_htseqHisat2_run1_fullset_reformatted.gct R2
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
numSamples=6
numRows=0
numCols=0
#Retrieve input filename and current number of columns
inFile="$inputsPath/GeneCounts_Formatted/$1"
numCols=$(($(head -n1 "$inFile" | awk '{print NF}')-1))
#Set output file name
countsFile=$(basename "$inFile" | sed 's/\.txt//g' | sed 's/\.csv//g' | sed 's/\.gct//g')
outFile="$outputsPath"/GeneCounts_Formatted/GeneCounts_Merged/"$countsFile"/subset"$2"_"$1"
#Set sample subset
subsetStart="$2"_VIS_Pool1
#Check file type
if [ ${inFile: -4} == ".gct" ]; then #GCT formatted
	#Retrieve column number of selected sample subset
	head -3 "$inFile" > tmpHeader.txt
	tail tmpHeader.txt | tr "\t" "\n" | grep -n "$subsetStart" > tmpColNum.txt
	colNum=$(($(cut -d ':' -f1 tmpColNum.txt)-3))
	#Retrieve the selected subset (default of 3 replicates for 2 treatments), 
	# including the first column with gene IDs and second with description
	cut -f 1,2,$colNum,$(($colNum+1)),$(($colNum+2)),$(($colNum+3)),$(($colNum+4)),$(($colNum+5)) "$inFile" > "$outFile"
	#Reset sample count for the second line of GCT formatted file
	head -1 tmpHeader.txt > tmpHeadLine.txt
	head -2 "$outFile" > tmpHeader.txt
	tail -1 tmpHeader.txt | sed "/\t$numCols/\t$numSamples/g" > tmpSampleLine.txt
	#Retrieve remaining data
	wc -l "$inFile" > tmpNumRows.txt
	numRows=$(($(cut -d ' ' -f1 tmpNumRows.txt)-2))
	tail -$numRows "$outFile" > tmpData.txt
	#Write separated lines to a single file
	cat tmpHeadLine.txt tmpSampleLine.txt tmpData.txt > "$outFile"
else #TXT and CSV formatted
	#Retrieve column number of selected sample subset
	head -1 "$inFile" | tr "\t" "\n" | grep -n "$subsetStart" > tmpColNum.txt
	colNum=$(($(cut -d ':' -f1 tmpColNum.txt)-1))
	#Retrieve the selected subset (default of 3 replicates for 2 treatments), 
	# including the first column with gene IDs
	cut -f 1,$colNum,$(($colNum+1)),$(($colNum+2)),$(($colNum+3)),$(($colNum+4)),$(($colNum+5)) "$inFile" > "$outFile"
fi
#Print a script completion confirmation message
echo "Selected $2 subset has been generated!"
#Clean up
rm tmp*.txt