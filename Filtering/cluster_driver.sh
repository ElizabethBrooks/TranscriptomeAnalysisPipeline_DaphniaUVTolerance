#!/bin/bash
#Script to perform sequence searches using a selected program for an input transcript data set
#Usage: bash cluster_driver.sh assembledFolder sampleList
#Usage Ex: bash cluster_driver.sh reciprocalNucleotide 0.90 trimmed_run1 Y05 Y023_5 E05 R2 PA Sierra
#Usage Ex: bash cluster_driver.sh protein 0.95 trimmed_run1 Y05 Y023_5 E05 R2 PA Sierra

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
	if [[ "$3" == trimmed* ]]; then
		inputFolder=$(echo "$3""$i"_assemblyTrinity)
	elif [[ "$3" == sorted* ]]; then
		inputFolder=$(echo "$3""$i"_assemblyGenomeTrinity)
	else
		echo "ERROR: Input folder for analysis is not a valid option... exiting!"
		exit 1
	fi
	#Check input option
	if [[ "$1" == reciprocal* ]]; then #Skip first three arguments 
		if [[ "$1" == protein ]]; then #Protein clustering
			if [ $counter -ge 3 ]; then
				#Usage: qsub reciprocalClusterProtein_cdhit.sh transcriptomeFasta clusterPercent
				echo "Reciprocal protein clustering $i transcriptome with cdhit at $2 percent identity..."
				qsub reciprocalClusterProtein_cdhit.sh "$inputFolder" "$2"
			fi
		elif [[ "$1" == nucleotide ]]; then #Nucleotide clustering
			if [ $counter -ge 3 ]; then
				#Usage: qsub reciprocalClusterNucleotide_cdhit.sh transcriptomeFasta clusterPercent
				echo "Reciprocal nucleotide clustering $i transcriptome with cdhit at $2 percent identity..."
				qsub reciprocalClusterNucleotide_cdhit.sh "$inputFolder" "$2"
			fi
		fi
	else
		if [[ "$1" == protein ]]; then #Protein clustering
			if [ $counter -ge 3 ]; then
				#Usage: qsub clusterProtein_cdhit.sh transcriptomeFasta clusterPercent
				echo "Protein clustering $i transcriptome with cdhit at $2 percent identity..."
				qsub clusterProtein_cdhit.sh "$inputFolder" "$2"
			fi
		elif [[ "$1" == nucleotide ]]; then #Nucleotide clustering
			if [ $counter -ge 3 ]; then
				#Usage: qsub clusterNucleotide_cdhit.sh transcriptomeFasta clusterPercent
				echo "Nucleotide clustering $i transcriptome with cdhit at $2 percent identity..."
				qsub clusterNucleotide_cdhit.sh "$inputFolder" "$2"
			fi
		fi
	fi
	counter=$(($counter+1))
done
