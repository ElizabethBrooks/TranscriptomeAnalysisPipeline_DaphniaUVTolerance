#!/bin/bash
#Script to perform t-test analysis of all samples in an input set
#Usage: bash searchDriver.sh assembledFolder sampleList
#Usage Ex: bash searchDriver.sh search trimmed_run1 Y05 Y023_5 E05 R2 PA Sierra
#Usage Ex: bash searchDriver.sh merge trimmed_run1 Y05 Y023_5 E05 R2 PA Sierra
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Initialize variables
counter=0
#Loop through all input sets of treatments and perform t-test analsysis
for i in "$@"; do
	#Determine what type of data folder was input
	if [[ "$2" == trimmed* ]]; then
		inputFolder=$(echo "$1""$i"_assemblyTrinity)
	elif [[ "$2" == sorted* ]]; then
		inputFolder=$(echo "$1""$i"_assemblyGenomeTrinity)
	else
		echo "ERROR: Input folder for analysis is not a valid option... exiting!"
		exit 1
	fi
	#Skip first two arguments
	if [ $counter -ge 2 ]; then
		#Check input option
		if [[ "$1" == search ]]; then
			#Usage: qsub search_blastp.sh transcriptomeFastaFolder
			qsub search_blastp.sh "$inputFolder"
		elif [[ "$1" == merge ]]; then
			#Usage: bash mergeSearches_blastp.sh transcriptomeFastaFolder
			bash mergeSearches_blastp.sh "$inputFolder"
		fi
	fi
	counter=$(($counter+1))
done