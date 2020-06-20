#!/bin/bash
#Script to perform t-test analysis of all samples in an input set
#Usage: bash search_driver.sh assembledFolder sampleList
#Usage Ex: bash search_driver.sh search trimmed_run1 swissprot Y05 Y023_5 E05 R2 PA Sierra
#Usage Ex: bash search_driver.sh reciprocal trimmed_run1 Y05 Y023_5 E05 R2 PA Sierra
#Usage Ex: bash search_driver.sh merge trimmed_run1 Y05 Y023_5 E05 R2 PA Sierra
#Usage Ex: bash search_driver.sh plot trimmed_run1

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#set output summary file path
outPath=$(grep "assemblyAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblyAnalysis://g")
#Set merged summary file name and header
if [[ "$1" == merge ]]; then
	#Set output file name
	outFile="$outPath"/"$2""_blastp_summary.txt"
	#Add header to output summary file
	echo "query,db,queryHits,reciprocalHits,bestHits" > "$outFile"
fi
#Initialize variables
counter=0
#Loop through all input sets of treatments and perform t-test analsysis
for i in "$@"; do
	#Determine what type of data folder was input
	if [[ "$2" == trimmed* ]]; then
		inputFolder=$(echo "$2""$i"_assemblyTrinity)
	elif [[ "$2" == sorted* ]]; then
		inputFolder=$(echo "$2""$i"_assemblyGenomeTrinity)
	else
		echo "ERROR: Input folder for analysis is not a valid option... exiting!"
		exit 1
	fi
	#Check input option
	if [[ "$1" == reciprocal ]]; then #Skip first two arguments 
		if [ $counter -ge 2 ]; then
			#Usage: qsub reciprocalSearch_blastp.sh transcriptomeFastaFolder
			echo "Searching $i transcriptome with blastp for $2..."
			qsub reciprocalSearch_blastp.sh "$inputFolder"
		fi
	elif [[ "$1" == merge ]]; then #Skip first two arguments
		if [ $counter -ge 2 ]; then
			#Usage: bash mergeSearches_blastp.sh transcriptomeFastaFolder
			echo "Merging $i blastp result for $2..."
			bash mergeSearches_blastp.sh "$inputFolder" "$i" >> "$outFile"
		fi
	elif [[ "$1" == search ]]; then #Skip first three arguments
		if [ $counter -ge 3 ]; then
			#Usage: qsub search_blastp.sh transcriptomeFastaFolder proteinDB
			echo "Searching $i transcriptome with blastp for $3..."
			qsub search_blastp.sh "$inputFolder" "$3"
		fi
	fi
	counter=$(($counter+1))
done
#Check if plotting was selected
if [ "$1" == plot ]; then
	#Load necessary modules
	module load R
	#Retrieve output file name
	outFile="$outPath"/"$2""_blastp_summary.txt"
	#Usage: Rscript blastpStats_barPlots.r title blastpSummaryFile
	echo "Plotting blastp results for $2..."
	Rscript blastpStats_barPlots.r "$2" "$outFile"
	echo "Blastp results for $2 have been plotted!"
fi
