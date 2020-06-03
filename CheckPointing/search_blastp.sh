#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N search_blastp_jobOutput
#$ -pe smp 8
#Script to use blastp to translate the nucleotide sequences of a reference genome
# for searching a protein database
#Usage: qsub search_blastp.sh transcriptomeFasta
#Usage Ex: qsub search_blastp.sh trimmed_run1E05_assemblyTrinity

#Load necessary modules for ND CRC servers
module load bio/blast+
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve genome reference absolute path for querying
dbPath=$(grep "transcriptomeDB:" ../InputData/inputPaths.txt | tr -d " " | sed "s/transcriptomeDB://g")
#Retrieve outputs absolute path
outputsPath=$(grep "transcriptSearch:" ../InputData/outputPaths.txt | tr -d " " | sed "s/transcriptSearch://g")
outputFolder="$outputsPath"/searchedTranscripts_"$1"
#Determine input database for blastp
if [[ "$1" == *assembly* ]]; then
	#Retrieve reads input absolute path
	assemblyPath=$(grep "assembling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assembling://g")
	inputsPath="$assemblyPath"/"$1"/decoded_transdecoder
	inputDB=Trinity.fasta.transdecoder.pep
	#Make blastable DB of transcriptome
	cd $inputsPath
	makeblastdb -in $inputDB -dbtype prot
	inputsPath="$inputsPath"/"$inputDB"
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
inputOutFile="$outputFolder"/searchedTranscripts_"$1"_summary.txt
#Use blastp to search a database
echo "Beginning blastp database search..."
blastp -query "$inputsPath" -db "$dbPath" -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 8 > blastp.outfmt6
echo "Finished blastp database search!"
#Output run commands to summary file
echo "blastp -query $inputsPath -db $dbPath  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 8 > blastp.outfmt6" >> "$inputOutFile"
#Switch query and search paths for reciprocal search
echo "Beginning reciprocal blastp database search..."
blastp -query "$dbPath" -db "$inputsPath" -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 8 > blastp.outfmt6
echo "Finished reciprocal blastp database search!"
#Output run commands to summary file
echo "blastp -query $dbPath -db $inputsPath  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 8 > blastp.outfmt6" >> "$inputOutFile"