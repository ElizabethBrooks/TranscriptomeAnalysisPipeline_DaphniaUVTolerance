#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N decoding_transdecoder_jobOutput
#$ -pe smp 8
#Script to use blastx to translate the nucleotide sequences of a reference genome
# for searching a protein database
#Usage: qsub searchProteins_blastx.sh databaseSelection
#Usage Ex: qsub searchProteins_blastx.sh ncbi

#Load necessary modules for ND CRC servers
module load bio/blast+
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Determine input database for blastx
if [[ "$1" == "ncbi" ]]; then
	#Set slected database to ncbi
	dbPath=$(grep "ncbiDB:" ../InputData/inputPaths.txt | tr -d " " | sed "s/ncbiDB://g")
elif [[ "$1" == "uniprot" ]]; then
	#Set slected database to uniprot
	dbPath=$(grep "uniprotDB:" ../InputData/inputPaths.txt | tr -d " " | sed "s/uniprotDB://g")
else
	#Error message
	echo "Invalid database selection entered (ncbi or uniprot only)... exiting!"
	exit 1
fi
#Retrieve genome reference absolute path for alignment
inputsPath=$(grep "genomeReference:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
#Retrieve outputs absolute path
outputsPath=$(grep "proteinSearch:" ../InputData/outputPaths.txt | tr -d " " | sed "s/proteinSearch://g")
outputFolder="$outputsPath"/searchedProteins_"$1"
mkdir "$outputFolder"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputFolder directory already exsists... please remove before proceeding."
	exit 1
fi
#Move to output folder
cd "$outputFolder"
#Name output file of inputs
inputOutFile="$outputFolder"/"$1"_searchedProteins_summary.txt
#Use blastx to search a protein database
echo "Beginning blastx protein database search..."
blastx -query "$searchFile" -db "$dbPath"  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 8 > blastx.outfmt6
echo "Finished blastx protein database search!"
#Output run commands to summary file
echo "blastx -query" "$searchFile" "-db" "$dbPath"  "-max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 8 >" "blastx.outfmt6" >> "$inputOutFile"
#Copy previous summaries
cp "$inputsPath"/"$1"/*.txt "$outputFolder"