#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N search_blastp_jobOutput
#$ -pe smp 8
#Script to use blastp to translate the nucleotide sequences of a reference genome
# for searching a protein database
#Usage: qsub search_blastp.sh transcriptomeFasta proteinDB
#Usage Ex: qsub search_blastp.sh trimmed_run1E05_assemblyTrinity swissprot
#Alternate usage Ex: qsub search_blastp.sh PA42 swissprot

#Load necessary modules for ND CRC servers
module load bio/blast+
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Determine database selection
if [[ "$2" == "ncbi" ]]; then
	#Retreive ncbi database storage path
	dbFile=$(grep "ncbiDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/ncbiDB://g")
elif [[ "$2" == "uniprot" ]]; then
	#Retreive uniprot database storage path
	dbFile=$(grep "uniprotDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/uniprotDB://g")
elif [[ "$2" == "swissprot" ]]; then
	#Retreive uniprot database storage path
	dbFile=$(grep "trinotateSwissprotDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/trinotateSwissprotDB://g")
else
	#Error message
	echo "Invalid database selection entered (ncbi or uniprot only)... exiting!"
	exit 1
fi
#Determine input query transcriptome for blastp
if [[ "$1" == *assembly* ]]; then
	#Retrieve reads input absolute path
	assemblyPath=$(grep "assembling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assembling://g")
	inputsPath="$assemblyPath"/"$1"/decoded_transdecoder
	#Set outputs absolute path
	outputFolder="$assemblyPath"/"$1"/searched_blastp_"$2"
	#Set input transcriptome path
	cd $inputsPath
	inputDB=Trinity.fasta.transdecoder.pep
	inputsPath="$inputsPath"/"$inputDB"
elif [[ "$1" == PA42 ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "transcriptomeDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/transcriptomeDB://g")
	#Set outputs absolute path
	outputPath=$(dirname "$inputsPath")
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
inputOutFile="$outputFolder"/searched_blastp_"$1"_summary.txt
#Use blastp to search a database
# and output with outfmt6 header:
#qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
echo "Beginning blastp database search..."
blastp -query "$inputsPath" -db "$dbFile" -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -num_threads 8 > blastp.outfmt6
echo "Finished blastp database search!"
#Output run commands to summary file
echo "blastp -query $inputsPath -db $dbFile  -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -num_threads 8 > blastp.outfmt6" >> "$inputOutFile"
