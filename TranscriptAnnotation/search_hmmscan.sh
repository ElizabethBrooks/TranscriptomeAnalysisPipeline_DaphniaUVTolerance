#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N search_hmmscan_jobOutput
#$ -pe smp 8
#Script to use hmmscan to use hmmer to identify protein domains
#Usage: qsub search_hmmscan.sh transcriptomeFasta
#Usage Ex: qsub search_hmmscan.sh trimmed_run1E05_assemblyTrinity
#Alternate usage Ex: qsub search_hmmscan.sh PA42
#Alternate usage Ex: qsub search_hmmscan.sh PA42_cds
#Alternate usage Ex: qsub search_hmmscan.sh PA42_transcripts

#Load necessary modules for ND CRC servers
module load bio/hmmer
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Retreive pfam database storage path
dbFile=$(grep "trinotatePfamDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/trinotatePfamDB://g")
#Determine input query transcriptome for blastp
if [[ "$1" == *assembly* ]]; then
	#Retrieve reads input absolute path
	assemblyPath=$(grep "assembling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assembling://g")
	inputsPath="$assemblyPath"/"$1"/decoded_transdecoder/Trinity.fasta.transdecoder.pep
	#Set outputs absolute path
	outputFolder="$assemblyPath"/"$1"/searched_hmmscan
elif [[ "$1" == PA42 ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "proteinSequencesDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/proteinSequencesDB://g")
	#Set outputs absolute path
	outputPath=$(dirname "$inputsPath")
	outputFolder="$outputPath"/searched_hmmscan
elif [[ "$1" == PA42_cds ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "codingSequencesDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/codingSequencesDB://g")
	inputsPath="$inputsPath"/decoded_transdecoder/PA42.3.0.cds_new.fasta.transdecoder.pep
	#Set outputs absolute path
	outputPath=$(dirname "$inputsPath")
	outputFolder="$outputPath"/searched_hmmscan	
elif [[ "$1" == PA42_transcripts ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "transcriptSequencesDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/transcriptSequencesDB://g")
	inputsPath="$inputsPath"/decoded_transdecoder/PA42.3.0.transcripts_new.fasta.transdecoder.pep
	#Set outputs absolute path
	outputPath=$(dirname "$inputsPath")
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
inputOutFile="$outputFolder"/searched_hmmscan_"$1"_summary.txt
#Use blastp to search a database
# and output with outfmt6 header:
#qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
echo "Beginning hmmscan pfam database search..."
hmmscan --cpu 8 --domtblout TrinotatePFAM.out "$dbFile" "$inputsPath" > pfam.log
echo "Finished hmmscan pfam database search!"
#Output run commands to summary file
echo "hmmscan --cpu 8 --domtblout TrinotatePFAM.out $dbFile $inputsPath > pfam.log" >> "$inputOutFile"
