#!/bin/bash
#Script to perform merge multifasta files and retain only
#the specified unique data (by sequence, ID, or both)
#Usage: bash uniqueMerge_fasta.sh mergeBy sortedFolder genotypes
#Usage Ex: bash uniqueMerge_fasta.sh sequence sortedCoordinate_samtoolsHisat2_run1 Y05 Y023_5 E05 R2 PA Sierra
#Usage Ex: bash uniqueMerge_fasta.sh sequence sortedCoordinate_samtoolsHisat2_run1 Y05 Y023_5 E05 R2 PA Sierra
#Default usage Ex: bash uniqueMerge_fasta.sh sequence assemblyTrinity_all

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder inputs supplied... exiting"
   	exit 1
fi
#Initialize variables
counter=1
fastaList=""
genotypeList=""
#Retrieve fasta output absolute path
outputsPath=$(grep "multiFASTA:" ../InputData/outputPaths.txt | tr -d " " | sed "s/multiFASTA://g")
#Determine input folder type
if [[ "$2" == *assembly* ]]; then
	#Retrieve fasta file path
	inputsPath=$(grep "assembling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assembling://g")
	#Retrieve selected fasta files
	#Loop through all input genotypes and merge fasta files
	for i in "$@"; do
		#Skip first two arguments
		if [[ $counter -eq $# ]]; then
			fastaList="$fastaList$inputsPath/$2$i""_assemblyGenomeTrinity/Trinity.fasta"
			genotypeList="$genotypeList$i"
		elif [[ $counter -ge 3 ]]; then
			fastaList="$fastaList$inputsPath/$2$i""_assemblyGenomeTrinity/Trinity.fasta "
			genotypeList="$genotypeList$i""_"
		fi
		counter=$(($counter+1))
	done
	#Determine assembly target
	if [[ "$2" == sorted* ]]; then
		#Create output directory
		outputFolder="$outputsPath/$2$genotypeList""_assemblyGenomeTrinity_multiFasta"
		#Set merged fasta file name
		multiFastaFile="$outputFolder/assemblyGenomeTrinity_multiFasta.fasta"
		summaryFile="$outputFolder/$2$genotypeList""_assemblyGenomeTrinity_multiFasta_summary.txt"
	else
		#Create output directory
		outputFolder="$outputsPath/$2$genotypeList""_assemblyTrinity_multiFasta"
		#Set merged fasta file name
		multiFastaFile="$outputFolder/assemblyTrinity_multiFasta.fasta"
		summaryFile="$outputFolder/$2$genotypeList""_assemblyTrinity_multiFasta_summary.txt"
	fi
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
	outputFolder="$outputsPath/$2_multiFasta"
	#Set merged fasta file name
	multiFastaFile="$outputFolder/$2""_multiFasta.fasta"
	summaryFile="$outputFolder/$2""_multiFasta_summary.txt"
fi
#Check if the folder already exists
mkdir "$outputFolder"
if [ $? -ne 0 ]; then
	echo "The $outputFolder directory already exsists... please remove before proceeding."
	exit 1
fi
#Move to outputs directory
cd "$outputFolder"
#Merge the set of fasta files
echo "Beginning fasta file merging..."
#Determine which method to merge fasta files by
if [[ "$1" == sequence ]]; then
	#Sequence identical merge
	awk 'BEGIN{RS=">"; FS="\n"; ORS=""}
		(FNR==1){next}
		{ name=$1; seq=$0; gsub(/(^[^\n]*|)\n/,"",seq) }
		!(seen[seq]++){ print ">" $0 }' $fastaList > $multiFastaFile
elif [[ "$1" == sequenceName ]]; then
	#First part of sequence name identical merge
	awk 'BEGIN{RS=">"; FS="\n"; ORS=""}
		(FNR==1){next}
		{ name=$1; seq=$0; gsub(/(^[^\n]*|)\n/,"",seq) }
		{ key=substr(name,1,index(s,"|")) }
		!(seen[key]++){ print ">" $0 }' $fastaList > $multiFastaFile
elif [[ "$1" == sequenceAndName ]]; then
	#Sequence name and sequence identical merge
	awk 'BEGIN{RS=">"; FS="\n"; ORS=""}
		(FNR==1){next}
		{ name=$1; seq=$0; gsub(/(^[^\n]*|)\n/,"",seq) }
		!(seen[name,seq]++){ print ">" $0 }' $fastaList > $multiFastaFile
else
	echo "Selected merge format for fasta files not valid... exiting!"
	exit 1
fi
echo "Fasta file merging complete!"
#Write list of fasta files to the summary file
echo "Fasta list: $fastaList" > $summaryFile
#Write fasta stats to the summary file
bash fastaStats.sh $fastaList $multiFastaFile >> $summaryFile