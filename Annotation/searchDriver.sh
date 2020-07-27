#!/bin/bash
#Script to perform sequence searches using a selected program for an input transcript data set
#Usage: bash search_driver.sh method PA42Target assembledFolder sampleList
#Usage Ex: bash search_driver.sh blastp PA42_proteins trimmed_run1 swissprot Y05 Y023_5 E05 R2 PA Sierra
#Usage Ex: bash search_driver.sh hmmscan PA42_proteins trimmed_run1 Y05 Y023_5 E05 R2 PA Sierra
#Usage Ex: bash search_driver.sh reciprocal PA42_proteins trimmed_run1 Y05 Y023_5 E05 R2 PA Sierra PA42_cds PA42_transcripts PA42_proteins

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#set output summary file path
outPath=$(grep "reciprocalSearch:" ../InputData/outputPaths.txt | tr -d " " | sed "s/reciprocalSearch://g")
#Initialize variables
counter=0
#Loop through all input sets of treatments and perform t-test analsysis
for i in "$@"; do
	#Determine what type of data folder was input
	if [[ "$3" == trimmed* ]]; then
		if [[ "$i" == PA42* ]]; then
			inputFolder="$i"
		else
			inputFolder=$(echo "$3""$i"_assemblyTrinity)
		fi
	elif [[ "$3" == sorted* ]]; then
		if [[ "$i" == PA42* ]]; then
			inputFolder="$i"
		else
			inputFolder=$(echo "$3""$i"_assemblyGenomeTrinity)
		fi
	else
		echo "ERROR: Input folder for analysis is not a valid option... exiting!"
		exit 1
	fi
	#Check input option
	if [[ "$1" == reciprocal ]]; then #Skip first three arguments 
		if [ $counter -ge 3 ]; then
			#Usage: qsub reciprocalSearch_blastp.sh transcriptomeFastaFolder
			echo "Searching $i transcriptome with blastp for $3 and $2..."
			qsub reciprocalSearch_blastp.sh "$inputFolder" "$2"
		fi
	elif [[ "$1" == hmmscan ]]; then #Skip first three arguments 
		if [ $counter -ge 3 ]; then
			#Usage: qsub search_hmmscan.sh transcriptomeFastaFolder
			echo "Searching $i transcriptome with hmmscan for $3 and $2..."
			qsub search_hmmscan.sh "$inputFolder" "$2"
		fi
	elif [[ "$1" == blastp ]]; then #Skip first three arguments
		if [ $counter -ge 3 ]; then
			#Usage: qsub search_blastp.sh transcriptomeFastaFolder proteinDB
			echo "Searching $i transcriptome with blastp using $2..."
			qsub search_blastp.sh "$inputFolder" "$2"
		fi
	fi
	counter=$(($counter+1))
done
