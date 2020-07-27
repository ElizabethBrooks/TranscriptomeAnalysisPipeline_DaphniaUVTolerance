#!/bin/bash
#Script to perform annotation of input transcript data sets
#Usage: bash annotation_driver.sh assembledFolder sampleList
#Usage Ex: bash annotation_driver.sh trimmed_run1 Y05 Y023_5 E05 R2 PA Sierra

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
	if [[ "$1" == trimmed* ]]; then
		inputFolder=$(echo "$1""$i"_assemblyTrinity)
	elif [[ "$1" == sorted* ]]; then
		inputFolder=$(echo "$1""$i"_assemblyGenomeTrinity)
	else
		echo "ERROR: Input folder for analysis is not a valid option... exiting!"
		exit 1
	fi
	#Skip first two arguments
	if [ $counter -ge 1 ]; then
		#Usage: qsub annotation_trinotate.sh transcriptomeFasta
		echo "Annotating $i transcriptome with trinotate..."
		qsub annotation_trinotate.sh "$inputFolder"
	fi
	counter=$(($counter+1))
done
