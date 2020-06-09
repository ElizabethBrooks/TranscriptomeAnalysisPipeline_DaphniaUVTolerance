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
#Retrieve fasta file path
inputsPath=$(grep "assembling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assembling://g")

#Determine assembly target
if [[ "$2" == sorted* ]]; then
	#Create output directory
	outputFolder="$outputsPath/$2""_assemblyGenomeTrinity_merged"$1"Fasta"
	#Set output file names
	mergedFastaFile="$outputFolder/Trinity.fasta"
	summaryFile="$outputFolder/$2""_assemblyGenomeTrinity_merged"$1"Fasta_summary.txt"
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
	#Create output directory
	outputFolder="$outputsPath/$2""_assemblyTrinity_merged"$1"Fasta"
	#Set output file names
	mergedFastaFile="$outputFolder/Trinity.fasta"
	summaryFile="$outputFolder/$2""_assemblyTrinity_merged"$1"Fasta_summary.txt"
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
else #Default accept a list of full file paths
	#Create output directory
	outputFolder="$outputsPath/$2_merged"$1"Fasta"
	#Set output file names
	mergedFastaFile="$outputFolder/Trinity.fasta"
	summaryFile="$outputFolder/$2""_merged"$1"Fasta_summary.txt"
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
fi
#Check if the folder already exists
mkdir "$outputFolder"
if [ $? -ne 0 ]; then
	echo "The $outputFolder directory already exsists... please remove before proceeding."
	#exit 1
fi

#Merge the set of fasta files
echo "Beginning fasta file merging..."
#bash fastaMerge.sh $1 $mergedFastaFile $fastaList

#Write fasta stats to the summary file
echo "Beginning file statistics summarizing..."
#bash fastaStats.sh $mergedFastaFile $fastaList > $summaryFile

#Write fasta stats to the csv formatted summary file
echo "Beginning file statistics formatting..."
summaryFileCSV=$(echo "$summaryFile" | sed 's/\.txt/\.csv/g')
bash fastaStats_csvFormatted.sh $mergedFastaFile $fastaList > $summaryFileCSV
#Fix file tags
genotypeList=$(echo "${@:3}")
counter=1
for i in $genotypeList; do
	genotypeTag="file"$counter
	echo "GENOTYPE: $genotypeTag, $i"
	sed -i "s/$genotypeTag/$i/g" $summaryFileCSV
	counter=$(($counter+1))
done
#Re-set header
sed -i 's/file/genotype/g' $summaryFileCSV

#Plot fasta stats from summary file
echo "Beginning file statistics plotting..."
#Rscript fastaStats_barPlots.r $1 $summaryFileCSV

