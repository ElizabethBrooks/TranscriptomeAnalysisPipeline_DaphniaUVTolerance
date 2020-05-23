#!/bin/bash
#Script to perform merge mergedfasta files and retain only
#the specified unique data (by sequence, ID, or both)
#Usage: bash fastaStats_driver.sh mergeBy sortedFolder genotypes
#Usage Ex: bash fastaStats_driver.sh sequence sortedCoordinate_samtoolsHisat2_run1 Y05 Y023_5 E05 R2 PA Sierra
#Usage Ex: bash fastaStats_driver.sh sequence sortedCoordinate_samtoolsTophat2_run1 Y05 Y023_5 E05 R2 PA Sierra
#Usage Ex: bash fastaStats_driver.sh sequence trimmed_run1 Y05 Y023_5 E05 R2 PA Sierra
#Default usage Ex: bash fastaStats_driver.sh sequence assemblyTrinity_all

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder inputs supplied... exiting"
   	exit 1
fi

#Initialize variables
counter=1
fastaList=""

#Retrieve fasta output absolute path
outputsPath=$(grep "mergedFASTA:" ../InputData/outputPaths.txt | tr -d " " | sed "s/mergedFASTA://g")

#Determine assembly target
if [[ "$2" == sorted* ]]; then
	#Retrieve fasta file path
	inputsPath=$(grep "assembling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assembling://g")
	#Create output directory
	outputFolder="$outputsPath/$2""_assemblyGenomeTrinity_mergedFasta"
	#Set merged fasta file name
	mergedFastaFile="$outputFolder/assemblyGenomeTrinity_mergedFasta.fasta"
	summaryFile="$outputFolder/$2""_assemblyGenomeTrinity_mergedFasta_summary.txt"
	#Retrieve selected fasta files
	#Loop through all input genotypes and merge fasta files
	for i in "$@"; do
		#Skip first two arguments
		if [[ $counter -eq $# ]]; then
			#Add fasta file to list
			fastaFile="$inputsPath/$2$i""_assemblyGenomeTrinity/Trinity.fasta"
			fastaList="$fastaList$fastaFile"
		elif [[ $counter -ge 3 ]]; then
			#Add fasta file to list
			fastaFile="$inputsPath/$2$i""_assemblyGenomeTrinity/Trinity.fasta "
			fastaList="$fastaList$fastaFile "
		fi
		counter=$(($counter+1))
	done
elif [[ "$2" == trimmed* ]]; then
	#Retrieve fasta file path
	inputsPath=$(grep "assembling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assembling://g")
	#Retrieve selected fasta files
	#Loop through all input genotypes and merge fasta files
	for i in "$@"; do
		#Skip first two arguments
		if [[ $counter -eq $# ]]; then
			#Add fasta file to list
			fastaFile="$inputsPath/$2$i""_assemblyTrinity/Trinity.fasta"
			fastaList="$fastaList$fastaFile"
		elif [[ $counter -ge 3 ]]; then
			#Add fasta file to list
			fastaFile="$inputsPath/$2$i""_assemblyTrinity/Trinity.fasta "
			fastaList="$fastaList$fastaFile "
		fi
		counter=$(($counter+1))
	done
	#Create output directory
	outputFolder="$outputsPath/$2""_assemblyTrinity_mergedFasta"
	#Set merged fasta file name
	mergedFastaFile="$outputFolder/assemblyTrinity_mergedFasta.fasta"
	summaryFile="$outputFolder/$2""_assemblyTrinity_mergedFasta_summary.txt"
else #Default accept a list of full file paths
	#Retrieve selected fasta files
	#Loop through all input genotypes and merge fasta files
	for i in "$@"; do
		#Skip first two arguments
		if [[ $counter -eq $# ]]; then
			fastaList="$fastaList$i"
		elif [[ $counter -ge 3 ]]; then
			fastaList="$fastaList$i "
		fi
		counter=$(($counter+1))
	done
	#Create output directory
	outputFolder="$outputsPath/$2_mergedFasta"
	#Set merged fasta file name
	mergedFastaFile="$outputFolder/$2""_mergedFasta.fasta"
	summaryFile="$outputFolder/$2""_mergedFasta_summary.txt"
fi

#Write list of fasta files to the summary file
echo "File list:" > $summaryFile
printf '%s\n' "${fastaList[@]}" >> $summaryFile

#Write fasta stats to the summary file
bash fastaStats.sh $fastaList $mergedFastaFile >> $summaryFile

#Write fasta stats to the summary file
summaryFileCSV=$(echo $summaryFile | sed 's/\.txt/\.csv/g')
bash fastaStats_csvFormatted.sh $fastaList $mergedFastaFile >> $summaryFileCSV

#Plot fasta stats from summary file
Rscript fastaStats_barPlot.r $summaryFileCSV