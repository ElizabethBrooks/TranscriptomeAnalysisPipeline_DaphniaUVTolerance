#!/bin/bash
#Script to perform merge multifasta files and retain only
#the specified unique data (by sequence, ID, or both)
#Usage: bash clusterFasta_cdHit.sh mergeBy sortedFolder genotypes
#Usage Ex: bash clusterFasta_cdHit.sh sequence sortedCoordinate_samtoolsHisat2_run1 Y05 Y023_5 E05 R2 PA Sierra
#Usage Ex: bash clusterFasta_cdHit.sh sequence sortedCoordinate_samtoolsHisat2_run1 Y05 Y023_5 E05 R2 PA Sierra
#Default usage Ex: bash clusterFasta_cdHit.sh sequence assemblyTrinity_all

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder inputs supplied... exiting"
   	exit 1
fi
#Initialize variables
counter=1
fastaList=""
#Retrieve cd-hit software package path
softsPath=$(grep "packageCdhit:" ../InputData/inputPaths.txt | tr -d " " | sed "s/packageCdhit://g")
#Retrieve fasta output absolute path
outputsPath=$(grep "multiFASTA:" ../InputData/outputPaths.txt | tr -d " " | sed "s/multiFASTA://g")
#Perform stringent clustering of sequences
echo "Beginning fasta file clustering..."
#Determine input folder source
#Determine assembly target
if [[ "$2" == sorted* ]]; then
	#Retrieve fasta file path
	inputsPath=$(grep "assembling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assembling://g")
	#Create output directory
	outputFolder="$outputsPath/$2""_assemblyGenomeTrinity_multiFasta"
	#Set merged fasta file name
	multiFastaFile="$outputFolder/assemblyGenomeTrinity_multiFasta.fasta"
	summaryFile="$outputFolder/$2""_assemblyGenomeTrinity_multiFasta_summary.txt"
	#Retrieve selected fasta files
	#Loop through all input genotypes and merge fasta files
	for i in "$@"; do
		#Skip first two arguments
		if [[ $counter -ge 3 ]]; then
			#Set current fasta file path
			fastaFile="$inputsPath/$2$i""_assemblyGenomeTrinity/Trinity.fasta "
			./"$softsPath"/cd-hit-est -o "$outputFolder"/cdhit -c 0.98 -i "$fastaFile" -p 1 -d 0 -b 3 -T 10
		fi
		counter=$(($counter+1))
	done
elif [[ "$2" == trimmed* ]]; then
	#Retrieve fasta file path
	inputsPath=$(grep "assembling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assembling://g")
	#Create output directory
	outputFolder="$outputsPath/$2""_assemblyTrinity_multiFasta"
	#Set merged fasta file name
	multiFastaFile="$outputFolder/assemblyTrinity_multiFasta.fasta"
	summaryFile="$outputFolder/$2""_assemblyTrinity_multiFasta_summary.txt"
	#Retrieve selected fasta files
	#Loop through all input genotypes and merge fasta files
	for i in "$@"; do
		#Skip first two arguments
		if [[ $counter -ge 3 ]]; then
			#Set current fasta file path
			fastaFile="$inputsPath/$2$i""_assemblyTrinity/Trinity.fasta "
			./"$softsPath"/cd-hit-est -o "$outputFolder"/cdhit -c 0.98 -i "$fastaFile" -p 1 -d 0 -b 3 -T 10
		fi
		counter=$(($counter+1))
	done
else #Default accept a list of full file paths
	#Create output directory
	outputFolder="$outputsPath/$2_multiFasta"
	#Set merged fasta file name
	multiFastaFile="$outputFolder/$2""_multiFasta.fasta"
	summaryFile="$outputFolder/$2""_multiFasta_summary.txt"
	#Retrieve selected fasta files
	#Loop through all input genotypes and merge fasta files
	for i in "$@"; do
		#Set current fasta file path
		if [[ $counter -ge 3 ]]; then
			./"$softsPath"/cd-hit-est -o "$outputFolder"/cdhit -c 0.98 -i "$i" -p 1 -d 0 -b 3 -T 10
		fi
		counter=$(($counter+1))
	done
fi
#Check if the folder already exists
mkdir "$outputFolder"
if [ $? -ne 0 ]; then
	echo "The $outputFolder directory already exsists... please remove before proceeding."
	exit 1
fi