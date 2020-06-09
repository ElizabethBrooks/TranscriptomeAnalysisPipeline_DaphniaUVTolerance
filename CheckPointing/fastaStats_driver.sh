#!/bin/bash
#Script to perform merge mergedfasta files and retain only
#the specified unique data (by sequence, ID, or both)
#Usage: bash fastaStats_driver.sh mergeBy sortedFolder genotypes
#Usage Ex: bash fastaStats_driver.sh sequenceAssembled sortedCoordinate_samtoolsHisat2_run1 Y05 Y023_5 E05 R2 PA Sierra
#Usage Ex: bash fastaStats_driver.sh sequenceAssembled sortedCoordinate_samtoolsTophat2_run1 Y05 Y023_5 E05 R2 PA Sierra
#Usage Ex: bash fastaStats_driver.sh sequenceDecoded trimmed_run1 Y05 Y023_5 E05 R2 PA Sierra
#Default usage Ex: bash fastaStats_driver.sh sequence assemblyTrinity_all

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder inputs supplied... exiting"
   	exit 1
fi

#Initialize variables
fastaList=""

#Retrieve fasta output absolute path
outputsPath=$(grep "mergedFASTA:" ../InputData/outputPaths.txt | tr -d " " | sed "s/mergedFASTA://g")
#Retrieve fasta file path
inputsPath=$(grep "assembling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assembling://g")

#Determine assembly target
if [[ "$2" == sorted* ]]; then
	#Create output directory
	outputFolder="$outputsPath/$2""_assemblyGenomeTrinity_"$1"_mergedFasta"
	#Set output file names
	mergedFastaFile="$outputFolder/merged_Trinity.fasta"
	summaryFile="$outputFolder/$2""_assemblyGenomeTrinity_"$1"_mergedFasta_summary.txt"
	#Retrieve selected fasta files
	#Loop through all input genotypes and add fasta files to a list
	for i in "${@:3}"; do #Skip first two arguments
		#Determine which fastas were input
		if [[ "$1" == *Assembled ]]; then
			#Add fasta file to list
			fastaFile="$inputsPath/$2$i""_assemblyGenomeTrinity/Trinity.fasta "
			fastaList="$fastaList$fastaFile"
			plotTitle="Trinity"
		elif [[ "$1" == *Decoded ]]; then
			#Add fasta file to list
			fastaFile="$inputsPath/$2$i""_assemblyGenomeTrinity/decoded_transdecoder/Trinity.fasta.transdecoder.pep "
			fastaList="$fastaList$fastaFile"
			plotTitle="Transdecoder"
		else
			echo "Invalid fasta input... exiting!"
			exit 1
		fi
	done
elif [[ "$2" == trimmed* ]]; then
	#Create output directory
	outputFolder="$outputsPath/$2""_assemblyTrinity_"$1"_mergedFasta"
	#Set output file names
	mergedFastaFile="$outputFolder/merged_Trinity.fasta"
	summaryFile="$outputFolder/$2""_assemblyTrinity_"$1"_mergedFasta_summary.txt"
	#Retrieve selected fasta files
	#Loop through all input genotypes and add fasta files to a list
	for i in "${@:3}"; do #Skip first two arguments
		#Determine which fastas were input
		if [[ "$1" == *Assembled ]]; then
			#Add fasta file to list
			fastaFile="$inputsPath/$2$i""_assemblyTrinity/Trinity.fasta "
			fastaList="$fastaList$fastaFile"
			plotTitle="Trinity"
		elif [[ "$1" == *Decoded ]]; then
			#Add fasta file to list
			fastaFile="$inputsPath/$2$i""_assemblyTrinity/decoded_transdecoder/Trinity.fasta.transdecoder.pep "
			fastaList="$fastaList$fastaFile"
			plotTitle="Transdecoder"
		else
			echo "Invalid fasta input... exiting!"
			exit 1
		fi
	done
else #Default accept a list of full file paths
	#Create output directory
	outputFolder="$outputsPath"/"$2"_"$1"_mergedFasta
	#Set output file names
	mergedFastaFile="$outputFolder"/merged_Trinity.fasta
	summaryFile="$outputFolder"/"$2""$1"_mergedFasta_summary.txt
	plotTitle="Input"
	#Retrieve selected fasta files
	#Loop through all input genotypes and add fasta files to a list
	for i in "${@:3}"; do #Skip first two arguments
		fastaList="$fastaList$i "
	done
fi
#Check if the folder already exists
mkdir "$outputFolder"
if [ $? -ne 0 ]; then
	echo "The $outputFolder directory already exsists... please remove before proceeding."
	exit 1
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
#bash fastaStats_csvFormatted.sh $mergedFastaFile $fastaList > $summaryFileCSV
#Fix file tags
genotypeList=$(echo "${@:3}")
genotypeList="Merged "$genotypeList
counter=1
for i in $genotypeList; do
	genotypeTag="file"$counter
	replaceTag=$genotypeTag","$i
	#sed -i "s/$genotypeTag/$replaceTag/g" $summaryFileCSV
	counter=$(($counter+1))
done
#sed -i 's/fileTotal,/fileTotal,Total,/g' $summaryFileCSV
#sed -i 's/fileDuplicates,/fileDuplicates,Duplicates,/g' $summaryFileCSV
#Re-set header
#sed -i 's/file,/file,genotype,/g' $summaryFileCSV

#Plot fasta stats from summary file
echo "Beginning file statistics plotting..."
egrep -v "Merged|Total|Duplicates" $summaryFileCSV > "$outputFolder"/tmp.csv
Rscript fastaStats_barPlots.r "$2" "$plotTitle" "$outputFolder"/tmp.csv
#Clean up
rm "$outputFolder"/tmp.csv
