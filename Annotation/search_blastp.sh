#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N search_blastp_jobOutput
#$ -pe smp 24
#Script to use blastp to translate the nucleotide sequences of a reference genome
# for searching a protein database
#Usage: qsub search_blastp.sh transcriptomeFasta proteinDB
#Usage Ex: qsub search_blastp.sh trimmed_run1E05_assemblyTrinity/clusteredNucleotides_cdhit_0.98 swissprot
#Alternate usage Ex: qsub search_blastp.sh PA42_proteins swissprot
#Alternate usage Ex: qsub search_blastp.sh PA42_cds swissprot
#Alternate usage Ex: qsub search_blastp.sh PA42_transcripts swissprot

#Load necessary modules for ND CRC servers
module load bio
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Determine database selection
if [[ "$2" == "ncbi" ]]; then
	#Retreive ncbi database storage path
	dbFile=$(grep "ncbi:" ../InputData/databasePaths.txt | tr -d " " | sed "s/ncbi://g")
elif [[ "$2" == "uniprot" ]]; then
	#Retreive uniprot database storage path
	dbFile=$(grep "uniprot:" ../InputData/databasePaths.txt | tr -d " " | sed "s/uniprot://g")
elif [[ "$2" == "swissprot" ]]; then
	#Retreive uniprot database storage path
	dbFile=$(grep "trinotateSwissprot:" ../InputData/databasePaths.txt | tr -d " " | sed "s/trinotateSwissprot://g")
else
	#Error message
	echo "Invalid database selection entered (ncbi or uniprot only)... exiting!"
	exit 1
fi
#Determine input query transcriptome for blastp
if [[ "$1" == *assemblyTrinity* ]]; then
	#Retrieve reads input absolute path
	assemblyPath=$(grep "assemblingFree:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingFree://g")
	inputsPath="$assemblyPath"/"$1"/decoded_transdecoder
	#Set outputs absolute path
	outputFolder="$assemblyPath"/"$1"/searched_blastp_"$2"
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
	outputFolder="$assemblyPath"/"$1"/searched_blastp_"$2"
	#Determine input file type
	if [[ "$1" == */clusteredNucleotide* ]]; then
		inputsPath="$inputsPath"/cdhitEst.transdecoder.pep
	elif [[ "$1" == */clusteredProtein* ]]; then
		inputsPath="$inputsPath"/cdhit.transdecoder.pep
	else 
		inputsPath="$inputsPath"/Trinity.fasta.transdecoder.pep
	fi
elif [[ "$1" == PA42_proteins ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "proteinSequences:" ../InputData/databasePaths.txt | tr -d " " | sed "s/proteinSequences://g")
	#Set outputs absolute path
	outputPath=$(dirname "$inputsPath")
	outputFolder="$outputPath"/searched_blastp_"$2"
elif [[ "$1" == PA42_cds ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "codingSequences:" ../InputData/databasePaths.txt | tr -d " " | sed "s/codingSequences://g")
	outputPath=$(dirname "$inputsPath")
	inputsPath="$outputPath"/decoded_transdecoder/PA42.3.0.cds_new.fasta.transdecoder.pep
	#Set outputs absolute path
	outputFolder="$outputPath"/searched_blastp_"$2"	
elif [[ "$1" == PA42_transcripts ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "transcriptSequences:" ../InputData/databasePaths.txt | tr -d " " | sed "s/transcriptSequences://g")
	outputPath=$(dirname "$inputsPath")
	inputsPath="$outputPath"/decoded_transdecoder/PA42.3.0.transcripts_new.fasta.transdecoder.pep
	#Set outputs absolute path
	outputFolder="$outputPath"/searched_blastp_"$2"	
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
inputOutFile="$outputFolder"/searched_blastp_summary.txt
#Use blastp to search a database
# and output with outfmt6 header:
#qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
echo "Beginning blastp database search..."
blastp -query "$inputsPath" -db "$dbFile" -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -num_threads 8 > blastp.outfmt6
echo "Finished blastp database search!"
#Output run commands to summary file
echo "blastp -query $inputsPath -db $dbFile  -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -num_threads 8 > blastp.outfmt6" >> "$inputOutFile"
