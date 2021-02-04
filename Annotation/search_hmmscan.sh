#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N search_hmmscan_jobOutput
#$ -pe smp 8
#Script to use hmmscan to use hmmer to identify protein domains
#Usage: qsub search_hmmscan.sh transcriptomeFasta
#Usage Ex: qsub search_hmmscan.sh trimmed_run1E05_assemblyTrinity/clusteredNucleotides_cdhit_0.98
#Alternate usage Ex: qsub search_hmmscan.sh PA42_v4.1_proteins
#Alternate usage Ex: qsub search_hmmscan.sh PA42_v4.1_cds
#Alternate usage Ex: qsub search_hmmscan.sh PA42_v4.1_transcripts

#Load necessary modules for ND CRC servers
module load bio
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Retreive pfam database storage path
dbFile=$(grep "trinotatePfam:" ../InputData/databasePaths.txt | tr -d " " | sed "s/trinotatePfam://g")
#Determine input query transcriptome for blastp
if [[ "$1" == *assemblyTrinity* ]]; then
	#Retrieve reads input absolute path
	assemblyPath=$(grep "assemblingFree:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingFree://g")
	inputsPath="$assemblyPath"/"$1"/decoded_transdecoder
	#Set outputs absolute path
	outputFolder="$assemblyPath"/"$1"/searched_hmmscan
	#Determine input file type
	if [[ "$1" == */clusteredNucleotide* ]]; then
		inputsPath="$inputsPath"/cdhitEst.transdecoder.pep
	elif [[ "$1" == */clusteredProtein* ]]; then
		inputsPath="$inputsPath"/cdhit.transdecoder.pep
	else 
		inputsPath="$inputsPath"/Trinity.fasta.transdecoder.pep
	fi
elif [[ "$1" == *assemblyGenome* ]]; then
	#Retrieve reads input absolute path
	assemblyPath=$(grep "assemblingGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingGenome://g")
	inputsPath="$assemblyPath"/"$1"/decoded_transdecoder
	#Set outputs absolute path
	outputFolder="$assemblyPath"/"$1"/searched_hmmscan
	#Determine input file type
	if [[ "$1" == */clusteredNucleotide* ]]; then
		inputsPath="$inputsPath"/cdhitEst.transdecoder.pep
	elif [[ "$1" == */clusteredProtein* ]]; then
		inputsPath="$inputsPath"/cdhit.transdecoder.pep
	else 
		inputsPath="$inputsPath"/Trinity.fasta.transdecoder.pep
	fi
elif [[ "$1" == *proteins ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "proteinSequences:" ../InputData/inputPaths.txt | tr -d " " | sed "s/proteinSequences://g")
	#Set outputs absolute path
	outputPath=$(dirname "$inputsPath")
	outputFolder="$outputPath"/searched_hmmscan
elif [[ "$1" == *cds ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "codingSequences:" ../InputData/inputPaths.txt | tr -d " " | sed "s/codingSequences://g")
	outputPath=$(dirname "$inputsPath")
	inputsPath=$(echo "$outputPath"/decoded_transdecoder/*transdecoder.pep)
	#Set outputs absolute path
	outputFolder="$outputPath"/searched_hmmscan	
elif [[ "$1" == *transcripts ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "transcriptSequences:" ../InputData/inputPaths.txt | tr -d " " | sed "s/transcriptSequences://g")
	outputPath=$(dirname "$inputsPath")
	inputsPath=$(echo "$outputPath"/decoded_transdecoder/*transdecoder.pep)
	#Set outputs absolute path
	outputFolder="$outputPath"/searched_hmmscan	
else
	#Error message
	echo "Invalid fasta entered (assembled transcriptome expected)... exiting!"
	exit 1
fi
#Make output directory
mkdir "$outputFolder"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputFolder directory already exsists... please remove before proceeding."
	exit 1
fi
#Move to output folder
cd "$outputFolder"
#Name output file of inputs
inputOutFile="$outputFolder"/searched_hmmscan_summary.txt
#Use hmmscan to search Pfam domain entries
echo "Beginning hmmscan pfam database search..."
hmmscan --cpu 8 --domtblout TrinotatePFAM.out "$dbFile" "$inputsPath" > pfam.log
echo "Finished hmmscan pfam database search!"
#Output run commands to summary file
echo "hmmscan --cpu 8 --domtblout TrinotatePFAM.out $dbFile $inputsPath > pfam.log" >> "$inputOutFile"
