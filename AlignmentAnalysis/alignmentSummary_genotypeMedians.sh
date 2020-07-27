#!/bin/bash
#Script to prepare alignment summary matrices
#Usage: bash alignmentSummary_genotypeMedians.sh analysisTarget runNum
#Usage Ex: bash alignmentSummary_genotypeMedians.sh tophat2 3 numGenotypes
#Alternate usage Ex: bash alignmentSummary_genotypeMedians.sh hisat2_run1 1 6

#Retrieve alignment analysis outputs absolute path
outputsPath=$(grep "alignmentAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/alignmentAnalysis://g")
#Set outputs directory and file
outDir="$outputsPath"/AlignmentsAnalyzed
outFile="$outDir"/alignmentSummarized_"$1"_run"$2"_medians.csv
outFileTmp="$outDir"/alignmentSummarized_"$1"_run"$2"_tmp.csv
#Set input path
inFile="$outDir"/alignmentSummarized_"$1"_run"$2"_formatted.csv
#Write header to output file
head -1 $inFile | sed 's/sample/genotype/g' > $outFile
#Insert new fields
startCol=$(cat $outFile | cut -d "," -f1)
endCols=$(cat $outFile | cut -d "," -f4-7)
echo $startCol",overallMedian,concordantMedian,overallSd,concordantSd,"$endCols > $outFile
#Determine the number of lines in the input file
numLines=0
numSamples=0
numLines=$(wc -l $inFile | cut -d " " -f1)
numLines=$(($numLines-1))
numSamples=$(($numLines/$3))
#Loop through sample subsets
for i in $(seq 1 $3); do
	#Calculate current window range
	start=$((($i-1)*$numSamples+2))
	end=$(($start+$numSamples-1))
	#Set current genotype tmp file
	currTmp="$outDir"/tmp"$start".csv
	#Retrieve current genotype data subset
	genotypeChunk=$(cat $inFile | awk "NR >= $start && NR <= $end { print }")
	echo $genotypeChunk > $currTmp
	#Fix formatting
	sed -i 's/% /%\n/g' $currTmp
	sed -i 's/run. /run.\n/g' $currTmp
	#Retrieve columns and calculate median values
	overallMedian=$(cut -d "," -f2 $currTmp | sed 's/%//g' | sort -n | awk -f ../util/median.awk)
	concordantMedian=$(cut -d "," -f3 $currTmp | sed 's/%//g' | sort -n | awk -f ../util/median.awk)
	#Retrieve columns to calculate standard deviations
	overallSd=$(cut -d "," -f2 $currTmp | sed 's/%//g' | sort -n | awk -f ../util/sd.awk)
	concordantSd=$(cut -d "," -f3 $currTmp | sed 's/%//g' | sort -n | awk -f ../util/sd.awk)
	#Output current genotype info
	startCol=$(head -1 $currTmp | cut -d "," -f1 | cut -d "_" -f1)
	endCols=$(tail -1 $currTmp | cut -d "," -f4-7)
	echo $startCol","$overallMedian","$concordantMedian","$overallSd","$concordantSd","$endCols >> $outFile
done
#Clean up
rm "$outDir"/tmp*.csv
